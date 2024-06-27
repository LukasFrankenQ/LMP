# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This scripts compiles model outputs from a date, period pair into output files ready to
be used for web visualisation.

**Outputs**

- ``RESULTS/half-hourly/{date}_{period}.json``

"""

import json
import pypsa
import shutil
import logging
import numpy as np
import pandas as pd
import geopandas as gpd
from pathlib import Path
from copy import deepcopy

logger = logging.getLogger(__name__)

from _helpers import configure_logging, to_datetime
from _aggregation_helpers import get_demand


def price_to_zones(n, regions):
    c = n.buses_t.marginal_price.columns
    return n.buses_t.marginal_price.iloc[0].loc[regions.index.intersection(c)]

def load_to_zones(n, regions):
    return (
        n.loads.loc[
            regions.index.intersection(n.loads.index), 'p_set'
            ]
            .rename('load')
    )

def p_nom_to_zones(n, regions):
    return (g := n.generators.groupby('bus')['p_nom'].sum()).loc[regions.index.intersection(g.index)]

def dispatch_to_zones(n, regions):
    con = pd.concat((n.generators, n.generators_t.p.T), axis=1)
    con = con.rename(columns={con.columns[-1]: 'dispatch'})
    con = con.groupby('bus').sum()['dispatch']
    return con.loc[regions.index.intersection(con.index)]

def get_generation_stack(n):
    return (
        pd.concat((
            n.generators[["carrier"]],
            (p := n.generators_t.p.T).rename(columns={p.columns[0]: "dispatch"})
        ), axis=1)
        .groupby("carrier").sum()
        [["dispatch"]]
    )

def get_gen_revenue(n):

    m = n.buses_t.marginal_price.iloc[0]

    gen = pd.concat([
        n.generators[['bus']],
        n.generators_t.p.T.rename(columns={n.generators_t.p.index[0]: 'dispatch'}),
        ], axis=1).groupby('bus').sum()
    
    return gen.multiply(m, axis=0).fillna(0).sum().sum()


def get_consumer_cost(n):

    m = n.buses_t.marginal_price.iloc[0]
    load = n.loads.p_set
    return load.multiply(m, axis=0).fillna(0).sum()


def get_cfd_cost(n, strike_prices):

    cfd_gens = n.generators.loc[
        n.generators.index.intersection(strike_prices.index),
        ['bus', 'carrier']
        ].copy()

    cfd_gens = pd.concat([
        cfd_gens,
        n.generators_t.p[cfd_gens.index].iloc[0].rename('dispatch'),
        cfd_gens['bus'].apply(
            lambda x: n.buses_t.marginal_price.iloc[0].loc[x]
            ).rename('marginal_price'),
        strike_prices.loc[cfd_gens.index, 'strike_price']
    ], axis=1)

    wholesale_cost = (cfd_gens['dispatch'] * cfd_gens['marginal_price']).sum()
    strike_price_cost = (cfd_gens['dispatch'] * cfd_gens['strike_price']).sum()

    cfd_cost = strike_price_cost - wholesale_cost

    return cfd_cost


if __name__ == "__main__":

    configure_logging(snakemake)

    date = snakemake.wildcards.date
    period = int(snakemake.wildcards.period)

    policy_settings = snakemake.params["policy_settings"]
    
    # strike_price = policy_settings["strike_price"]
    strike_prices = pd.read_csv(
        snakemake.input["strike_prices"], index_col=0
    )

    consumer_rent_share = policy_settings["consumer_rent_share"]

    cost_factors = [
        "renewable obligation",
        "feed-in tariffs",
        "old feed-in tariffs",
        "energy company obligation",
        "aahedc",
        "backwardation",
        "cfd",
    ]
    cost_factor_mapper = {}

    index_mapper = {
        'domestic single': slice(0,48),
        'domestic multi': slice(52,100),
        'non-domestic single': slice(104,152),
        'non-domestic multi': slice(156,204),
    }

    single_rate_domestic = pd.read_excel(
        snakemake.input["elexon_demand_profiles"],
        index_col=0,
        header=2).iloc[index_mapper['domestic single']]
    multi_rate_domestic = pd.read_excel(
        snakemake.input["elexon_demand_profiles"],
        index_col=0,
        header=2).iloc[index_mapper['domestic multi']]

    single_rate_nondomestic = pd.read_excel(
        snakemake.input["elexon_demand_profiles"],
        index_col=0,
        header=2).iloc[index_mapper['non-domestic single']]
    multi_rate_nondomestic = pd.read_excel(
        snakemake.input["elexon_demand_profiles"],
        index_col=0,
        header=2).iloc[index_mapper['non-domestic multi']]

    layouts = ['national', 'nodal', 'eso']

    results = {layout: {'geographies': {}} for layout in layouts}

    # treating nodal differently, as it is assumed to be the layout
    # the requires no further balancing
    nodal_generation_stack = (
        get_generation_stack(
            pypsa.Network(snakemake.input["network_nodal"]))
    )

    national_model = pypsa.Network(snakemake.input["network_{}".format('national')])
    national_cfd_cost = get_cfd_cost(national_model, strike_prices)

    carrier_mapper = {
        "PHS": "hydro",
        "hydro-scheme": "hydro",
        "PHS": "hydro",
        "floating wind": "wind",
        "onwind": "wind",
        "offwind": "wind",
        "CCGT": "thermal",
        "CHP": "thermal",
        "biomass": "thermal",
        "gas": "thermal",
        "gas-fired": "thermal",
        "gas turbine": "thermal",
        "coal": "thermal",
        "powerstation": "thermal",
    }

    redispatch_costs = pd.read_csv(
        snakemake.input["redispatch_cost"],
        index_col=0,
        parse_dates=True
        )

    if (ts := pd.Timestamp(date)) in redispatch_costs.index:

        offer_costs = redispatch_costs.loc[pd.Timestamp(date), "offers"]
        bid_costs = redispatch_costs.loc[pd.Timestamp(date), "bids"]
    
    elif ts > redispatch_costs.index.max():
        logger.warning("Date is beyond redispatch costs data. Taking average of last 31 days.")
        offer_costs, bid_costs = redispatch_costs.iloc[-31:].mean().values

    elif ts < redispatch_costs.index.min():
        logger.warning("Date is before redispatch costs data. Taking average of first 31 days.")
        offer_costs, bid_costs = redispatch_costs.iloc[:31].mean().values

    result_store_regional = {}
    result_store_global = {}

    assert layouts[0] == 'national', "Assumes national layout is first in list."
    # for layout in set(layouts):
    for layout in layouts:

        logger.info("Summarising layout {} for {} {}.".format(layout, date, period))

        n = pypsa.Network(snakemake.input["network_{}".format(layout)])

        genstack = get_generation_stack(n)

        diff = (genstack - nodal_generation_stack).iloc[:,0]

        def clean_carriers(d):

            for origin, target in carrier_mapper.items():
                if origin not in d.index:
                    continue

                try:
                    d.loc[target] += d.loc[origin]
                except KeyError:
                    d.loc[target] = d.loc[origin]

                d.drop(index=origin, inplace=True)

            return d

        diff = clean_carriers(diff)

        bid_volume = diff.loc[diff > 0].sum()
        offer_volume = diff.loc[diff < 0].abs().sum()

        bcost = bid_volume * bid_costs + offer_volume * offer_costs

        if layout != 'national':
            bcost = min(bcost, result_store_global['national']['balancing_cost'])

        # bids = bids.loc[bids.index.intersection(redispatch_costs.index)]
        # offers = offers.loc[offers.index.intersection(redispatch_costs.index)]

        congestion_rent = get_gen_revenue(n) - get_consumer_cost(n) # is negative as it reduces consumer payment
        cfd_cost = get_cfd_cost(n, strike_prices) # should be positive (in most cases), increasing consumers payment

        regions_file = snakemake.input["regions_{}".format(layout)]
        regions = gpd.read_file(regions_file).set_index("name")

        regional_results = pd.DataFrame(index=regions.index)

        # object to store wholesale cost, balancing cost, congestion rents and cfd payments
        global_results = {
            "wholesale_cost": get_consumer_cost(n),
            "balancing_cost": bcost,
            "congestion_rent": congestion_rent,
            "cfd_cost": cfd_cost,
        }

        regional_results.loc[:, "wholesale_price"] = price_to_zones(n, regions)
        regional_results.loc[:, "load"] = load_to_zones(n, regions)# .values

        G = regional_results["load"].sum()

        regional_results.loc[:, "post_policy_price"] = (
            (
                regional_results["wholesale_price"] * G
                + global_results["balancing_cost"]
                + global_results["congestion_rent"]
                + global_results["cfd_cost"]
            ) / G
        )

        result_store_regional[layout] = regional_results
        result_store_global[layout] = global_results
        
        # cost savings
        srd = get_demand(single_rate_domestic, date, period)
        mrd = get_demand(multi_rate_domestic, date, period)
        srn = get_demand(single_rate_nondomestic, date, period)
        mrn = get_demand(multi_rate_nondomestic, date, period)

        national_ws = result_store_regional['national']['wholesale_price'].iloc[0]
        national_pb = result_store_regional['national']['post_policy_price'].iloc[0]
        
        for name, demand in zip(
            ['single-rate-domestic', 'multi-rate-domestic', 'single-rate-nondomestic', 'multi-rate-nondomestic'],
            [srd, mrd, srn, mrn]):

            regional_results.loc[:, name + "_wholesale_savings"] = (
                (national_ws - 
                regional_results['wholesale_price']) * demand * 1e-3
                )
            regional_results.loc[:, name + "_total_savings"] = (
                (national_pb - 
                regional_results['post_policy_price']) * demand * 1e-3
                )

        for region in regions.index:

            results[layout]['geographies'][region] = {
                "variables": regional_results.T[region].fillna(0.).astype(np.float32).to_dict()
            }

        global_hold = deepcopy(global_results)

        for cost_factor, value in global_results.items():

            national_value = result_store_global["national"][cost_factor]
            global_hold[cost_factor + '_savings'] = national_value - value
        
        global_results = global_hold

        results[layout]['globals'] = {'variables': global_results}


    with open(snakemake.output[0], "w") as f:
        json.dump({int(to_datetime(date, period).timestamp()): results}, f)

    # clean up resources to reduce disk space
    resource_path = Path(snakemake.input["network_nodal"]).parent
    shutil.rmtree(resource_path)
