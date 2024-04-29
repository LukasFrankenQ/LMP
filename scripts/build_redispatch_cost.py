# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
#
# SPDX-License-Identifier: MIT
"""
Iterates over a BMUs and their respective balancing actions to estimate the cost of redispatching
The dataset winter_2024_balancing_actions.csv contains accepted balancing actions between January 2024 and April 2024
that were SO-flagged, i.e. caused by transmission constraints. The script calculates the weighted average of the
costs of the balancing actions for each BMU.

**Inputs**

- ``data/winter_2024_balancing_actions.csv``: transmission constraint related accepted balancing actions
- ``RESOURCES/networks/gen.nc``: pypsa network including modelled BMUs. Used to assign carriers to generators
- ``RESOURCES/bmu_cost_estimates.csv``, estimated costs of electricity generation for each BMU in wholesale

**Outputs**

- ``RESOURCES/redispatch_cost.csv``: 

"""

import logging

logger = logging.getLogger(__name__)

import pypsa
import pandas as pd

from _helpers import configure_logging

idx = pd.IndexSlice


if __name__ == "__main__":

    configure_logging(snakemake)

    balancing_actions = pd.read_csv(
        snakemake.input["balancing_actions"],
        index_col=[0,1],
        parse_dates=True
        )
    n = pypsa.Network(snakemake.input["network"])
    bmus_costs = pd.read_csv(snakemake.input["wholesale_cost_estimates"], index_col=0).iloc[:,0]

    revs = balancing_actions.loc[idx[:, ["offer revenue", "bid revenue"]],:]
    prices = balancing_actions.loc[idx[:, ["offer price", "bid price"]],:]

    results = pd.DataFrame(columns=["offer_cost", "offer_vol", "bid_cost", "bid_vol"])

    for bmu in revs.index.get_level_values(0).unique():

        rev = revs.loc[idx[bmu, :], :]    
        price = prices.loc[idx[bmu, :], :]

        rev = rev.loc[:, ~rev.isna().all(axis=0)]
        price = price.loc[:, rev.columns]

        vol = pd.DataFrame(rev.values / price.values, index=rev.index, columns=rev.columns)
        offer_col, bid_vol = vol.sum(axis=1).values

        vol = vol.div(vol.sum(axis=1), axis=0)

        weighted_avg = (vol.values * price.values).sum(axis=1)
        results.loc[bmu] = pd.Series({
            "offer_cost": weighted_avg[0],
            "offer_vol": offer_col,
            "bid_cost": weighted_avg[1],
            "bid_vol": bid_vol,
        })
    
    ren_carriers = [
        'onwind',
        'offwind',
        'wind',
        'hydro',
        'hydro-scheme',
        'floating wind',
    ]

    offer_results = results.loc[results["offer_vol"] > 0]
    bid_results = results.loc[results["bid_vol"] > 0]

    offinter = offer_results.index.intersection(bmus_costs.index)
    bidinter = bid_results.index.intersection(bmus_costs.index)

    renewable_bidinter = (
        n.generators
        .loc[n.generators.carrier.isin(ren_carriers)]
        .index.intersection(bid_results.index)
    )

    offer_total = pd.concat((
        bmus_costs.rename('wholesale cost').loc[offinter], 
        offer_results.loc[offinter, ['offer_cost', 'offer_vol']],
    ), axis=1)

    dispatch_bid_total = pd.concat((
        bmus_costs.rename('wholesale cost').loc[bidinter], 
        bid_results.loc[bidinter, ['bid_cost', 'bid_vol']],
    ), axis=1)

    renewable_bid_total = pd.concat((
        n.generators.loc[renewable_bidinter, 'marginal_cost'].rename('wholesale cost'), 
        bid_results.loc[renewable_bidinter, ['bid_cost', 'bid_vol']],
    ), axis=1)

    total = pd.concat((offer_total, dispatch_bid_total, renewable_bid_total), axis=0).drop_duplicates()
    total = total.groupby(total.index).first()
    total = (
        pd.concat((
            total,
            n.generators.loc[total.index.intersection(n.generators.index), 'carrier']
            ), axis=1)
    )

    total.groupby('carrier').median().to_csv(snakemake.output[0])

