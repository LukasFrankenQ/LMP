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

import os
import json
import pypsa
import logging
import numpy as np
import pandas as pd
import geopandas as gpd

logger = logging.getLogger(__name__)

from _helpers import configure_logging, to_datetime
    

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


if __name__ == "__main__":

    configure_logging(snakemake)

    date = snakemake.wildcards.date
    period = int(snakemake.wildcards.period)

    layouts = ['national', 'nodal', 'fti', 'eso']
    results = {layout: {'geographies': {}} for layout in layouts}

    # treating nodal differently, as it is assumed to be the layout
    # the requires no further balancing
    nodal_generation_stack = (
        get_generation_stack(
            pypsa.Network(snakemake.input["network_nodal"]))
    )

    balancing_cost = snakemake.params["balancing"]

    logger.warning("Unsolved issue with dispatch and p_nom for nodal layout.")

    assert layouts[0] == 'national', "Assumes national layout is first in list."
    for layout in set(layouts):

        logger.info("Summarising layout {} for {} {}.".format(layout, date, period))

        n = pypsa.Network(snakemake.input["network_{}".format(layout)])

        genstack = get_generation_stack(n)
        diff = (genstack - nodal_generation_stack).abs().iloc[:,0]
        
        bcost = pd.Series(
                    map(
                        lambda carrier: balancing_cost.get(carrier, balancing_cost['default']),
                        diff.index  
                    ),
                    index=diff.index
                ).mul(diff, axis=0).sum()
        bvol = diff.sum()

        regions_file = snakemake.input["regions_{}".format(layout)]
        regions = gpd.read_file(regions_file).set_index("name")

        layout_results = pd.DataFrame(index=regions.index)

        layout_results.loc[:, "marginal_price"] = price_to_zones(n, regions)
        layout_results.loc[:, "load"] = load_to_zones(n, regions)# .values
        layout_results.loc[:, "available_capacity"] = p_nom_to_zones(n, regions)
        layout_results.loc[:, "dispatch"] = dispatch_to_zones(n, regions)# .values

        G = layout_results["load"].sum()

        layout_results.loc[:, "post_balancing_price"] = (
            (layout_results['marginal_price'] * G + bcost) / G
        )

        for region in regions.index:

            results[layout]['geographies'][region] = {
                "variables": layout_results.T[region].fillna(0.).astype(np.float16).to_dict()
            }

    with open(snakemake.output[0], "w") as f:
        json.dump({int(to_datetime(date, period).timestamp()): results}, f)

    os.remove(snakemake.input["regions_nodal"])
