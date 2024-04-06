# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This rule makes a summary of prices for the wildcards passed in the config file

**Outputs**

- ``RESOURCES/plots/prices.pdf``: comparison of wholesale prices
- ``RESOURCES/plots/generation.pdf``: comparison of the generation stack
"""

import logging

import pypsa
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

from itertools import product

logger = logging.getLogger(__name__)

from _helpers import configure_logging, process_scenarios


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

def generation_to_zones(n, regions):
    con = pd.concat((n.generators, n.generators_t.p.T), axis=1)
    con = con.rename(columns={con.columns[-1]: 'dispatch'})
    con = con.groupby('bus').sum()['dispatch']
    return con.loc[regions.index.intersection(con.index)]


if __name__ == "__main__":

    configure_logging(snakemake)

    logger.info("Gathering results for scenarios: %s", snakemake.config["scenario"])

    scenarios = process_scenarios(snakemake.config["scenario"])

    layouts = scenarios["layouts"]

    axs_dict = {
        output:
        plt.subplots(1, len(scenarios["layouts"]), figsize=(15, 3.5*len(scenarios["layouts"])))
        for output in ["prices", "generation"]
        }
    
    for i, layout in enumerate(layouts):

        networks_dict = (
            {
                (date, period): 
                snakemake.params.RESOURCES
                + "/live_data/"
                + f"{date}_{period}/"
                + f"network_s_{layout}_solved.nc"

                for date in scenarios["dates"]
                for period in scenarios["periods"]
            }
        )

        index = pd.MultiIndex.from_tuples(
            networks_dict.keys(),
            names=["date", "period"]
        )

        real_prices = [pd.read_csv(
            snakemake.params.RESOURCES + "/live_data/" + f"{date}_{period}/price_stats.csv", index_col=0)
            .loc["market_index"].iloc[0] for date, period in networks_dict.keys()]
        real_prices = pd.Series(real_prices, index=index)

        results = {
            output: pd.DataFrame(index=index)
            for output in ["prices", "generation"]
        }

        buses = pypsa.Network(networks_dict.values()[0]).buses.index

        for j, fn in networks_dict.items():
            for output in results.keys():

                n = pypsa.Network(fn)
                get_func = globals()[f"{output}_to_zones"]

                results[output].loc[j] = get_func(n, buses)


        def get_date(date, period):
            return pd.Timestamp(date) + pd.Timedelta(hours=0.5) * (period - 1)

        for output, (fig, axs) in axs_dict.items():

            ax = axs[i]

            results[output].index = (
                results[output].index
                .to_frame()
                .apply(lambda row: get_date(row.date, row.period), axis=1)
            )
            real_prices.index = results[output].index

            detailed_plot = len(results.columns) < 9

            if not detailed_plot:
                results[output].plot(ax, legend=False, color='k', alpha=0.4, linewidth=0.5)

            logger.warning("Should be weighted mean for marginal prices!")
            results[output].mean(axis=1).plot(
                ax=ax,
                linewidth=1.5,
                alpha=1,
                color='darkred',
                label="Avg Marginal Price"
                )

            if output == "prices":
                real_prices.plot(
                    ax=ax,
                    linewidth=1.5,
                    alpha=1,
                    color='blue',
                    label="Market Index"
                    )

            ax.legend()
        
    
