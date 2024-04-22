# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This rule collects live and other data to prepare the network for live operation.

**Outputs**

- ``RESOURCES/networks/prepared_live_{data}_{period}.csv``: prepared network

"""

import logging

import pypsa
import pandas as pd

from _helpers import configure_logging, check_network_consistency

logger = logging.getLogger(__name__)

if __name__ == "__main__":

    configure_logging(snakemake)

    date = snakemake.wildcards.date    
    period = snakemake.wildcards.period

    logger.info(f"Preparing live network for {date} settlement period {period}.")
    n = pypsa.Network(snakemake.input["network"])

    isolated_buses = check_network_consistency(n)
    logger.info(f"A total of {len(isolated_buses)} isolated buses:\n" + ",".join(isolated_buses))

    logger.warning("Should be Export Limit for dispatchable generators, but not yet implemented!")
    bmu = pd.read_csv(snakemake.input["elexon_bmus"]).set_index("NationalGridBmUnit")

    load_weights = pd.read_csv(snakemake.input["load_weights"], index_col=0)["load_weight"]

    total_load = bmu.loc[bmu["PN"] > 0, "PN"].sum()

    logger.info("Moving exporting interconnectors from load to generators.")
    interconnectors = n.links[n.links.carrier == "AC"]

    export = bmu.loc[(
        bmu.index.str.startswith('I') &
        (bmu.index.str[3] == '-') &
        bmu.index.str[:3].str.isalpha() &
        (bmu["PN"] < 0)
        ), "PN"]
    
    total_export = abs(export.sum())

    logger.info(f"Total export capacity: {total_export} MW. Removed from load.")

    total_load -= total_export

    loads = load_weights.loc[load_weights > 0].index

    n.madd(
        "Load",
        loads,
        bus=loads,
        p_set=load_weights.loc[loads] * total_load,
    )

    bmu = bmu.loc[n.generators.index]

    mask = n.generators.carrier.str.contains('wind')

    wind_idx = n.generators.index[mask].intersection(bmu.loc[bmu['PN'] > 0.].index)
    others_idx = n.generators.index[~mask].intersection(bmu.loc[bmu['PN'] > 0.].index)

    logger.info("Inserting Export Limit as p_nom for non-wind BMUs.")
    
    # these carriers need cleaning up
    melp_plants = [
        'gas-fired',
        'gas',
        'biomass',
        'battery',
        'CHP',
        'CCGT',
        'coal',
        'gas turbine',
        'powerstation',
        ]

    n.generators.loc[others_idx, 'p_nom'] = bmu.loc[bmu["PN"] > 0].max(axis=1).loc[others_idx]

    disp_idx = n.generators.index[n.generators.carrier.isin(melp_plants)].intersection(bmu.index)
    n.generators.loc[disp_idx, 'p_nom'] = bmu.loc[disp_idx].max(axis=1)

    logger.info("Inserting Physical Notification as p_nom for wind BMUs. To these bid volumes are added.")
    n.generators.loc[wind_idx, 'p_nom'] = bmu.loc[wind_idx, "PN"]

    # adding bid volumes to respective wind BMUs
    balancing_actions = pd.read_csv(snakemake.input["real_balancing_actions"], index_col=0)

    inside = balancing_actions.index.intersection(n.generators.index)
    outside = balancing_actions.index.difference(n.generators.index)

    n.generators.loc[inside, 'p_nom'] += balancing_actions.loc[inside, 'bid volume']
    
    # Simulate export by enforcing negative generation at interconnectors
    missing = export.loc[export.index.difference(n.generators.index)]
    if len(missing):
        logger.info("Dropping missing exporting BMUs {}.".format(",".join(missing.index)))
        export = export.drop(missing.index)    

    n.generators.loc[(g := n.generators.loc[export.index]).loc[g["p_nom"] == 0.].index, "p_nom"] = export.abs()

    pu = export.div(n.generators.loc[export.index, 'p_nom'])

    n.generators.loc[export.index, 'p_max_pu'] = pu
    n.generators.loc[export.index, 'p_min_pu'] = pu

    logger.info("Adjusting marginal costs for dispatchable generators according to wholesale price.")
    real_price = pd.read_csv(snakemake.input["price_stats"], index_col=0).iloc[0,0]

    method_cutoff = snakemake.params["elexon"]["cost_assignment_method_cutoff"]
    
    method = 'adjust_' + ['dispatch_only', 'all'][int(real_price < method_cutoff)]

    cost_estimated_generators = pd.Index(
        pd.read_csv(snakemake.input["cost_estimated_generators"], index_col=0).iloc[:,0].tolist()
        )

    export_volume = n.generators.loc[export.index, 'p_nom'].sum()

    above_price = (
        (m := n.generators.drop(export.index)[['marginal_cost', 'p_nom']].sort_values(by='marginal_cost'))
        .loc[m['p_nom'].cumsum() >= n.loads.p_set.sum() + export_volume]
    )
    price_setter = above_price.index[0]

    if method == 'adjust_dispatch_only':

        factor = real_price / above_price.at[price_setter, 'marginal_cost']
        n.generators.loc[cost_estimated_generators, 'marginal_cost'] *= factor
    
    elif method == 'adjust_all':

        cost_subtraction = n.generators.loc[price_setter, 'marginal_cost'] - real_price
        n.generators['marginal_cost'] -= cost_subtraction


    isolated_buses = check_network_consistency(n)
    logger.info(f"A total of {len(isolated_buses)} isolated buses:\n" + ",".join(isolated_buses))

    n.generators = n.generators.loc[n.generators.p_nom > 0]
    n.export_to_netcdf(snakemake.output["network"])