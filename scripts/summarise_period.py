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
import logging
import pandas as pd
import geopandas as gpd

logger = logging.getLogger(__name__)

from _helpers import configure_logging
    

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


if __name__ == "__main__":

    configure_logging(snakemake)

    date = snakemake.wildcards.date
    period = int(snakemake.wildcards.period)

    layouts = ['nodal', 'national', 'fti', 'eso']
    results = {layout: {} for layout in layouts}

    for layout in layouts:

        logger.info("Summarising layout {} for {} {}.".format(layout, date, period))

        n = pypsa.Network(snakemake.input["network_{}".format(layout)])
        regions = gpd.read_file(snakemake.input["regions_{}".format(layout)]).set_index("name")

        layout_results = pd.DataFrame(index=regions.index)

        layout_results.loc[:, "marginal_price"] = price_to_zones(n, regions)
        layout_results.loc[:, "load"] = load_to_zones(n, regions)
        layout_results.loc[:, "available_capacity"] = p_nom_to_zones(n, regions)
        layout_results.loc[:, "dispatch"] = dispatch_to_zones(n, regions)

        for region in regions.index:

            results[layout][region] = {
                "variables": layout_results.T[region].fillna(0.).to_dict()
            }

    with open(snakemake.output[0], "w") as f:
        json.dump(results, f)
