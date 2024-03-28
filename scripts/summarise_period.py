# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This rule makes a summary of the period comparing modelled to actual wholesale prices.

**Outputs**

- ``RESOURCES/{date}_{period}/summary.csv``: the summary
- ``RESOURCES/{date}_{period}/maps.csv``: maps showing the prices in different layouts and the market price

"""

import logging

import pypsa
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

from _helpers import configure_logging

if __name__ == "__main__":

    configure_logging(snakemake)

    layouts = ['nodal', 'fti', 'eso']

    stats = pd.read_csv(snakemake.input["price_stats"], index_col=0)

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

    vmins = {
        'price': 0,
    }
    vmaxs = {
        'price': 100,
    }

    for metric in ['dispatch', 'load', 'price', 'p_nom']:    

        get_func = globals()[f"{metric}_to_zones"]

        fig, axs = plt.subplots(1, len(layouts), figsize=(15, 5))

        for layout, ax in zip(layouts, axs):

            n = pypsa.Network(snakemake.input["network_{}".format(layout)])
            regions = gpd.read_file(snakemake.input["regions_{}".format(layout)]).set_index("name")

            regions.loc[:, metric] = get_func(n, regions)            
            regions[metric] = regions[metric].fillna(0)

            regions.plot(
                column=metric,
                legend=True,
                ax=ax,
                vmin=vmins.get(metric, None),
                vmax=vmaxs.get(metric, None),
                # label='Price in £/MWh',
                edgecolor='k',
                linewidth=0.5
                )

            ax.set_xticks([])
            ax.set_yticks([])

        axs[0].set_title('Modelled {}'.format(metric))

        if metric == 'price':
            axs[1].set_title('Wholesale Market Index: {} £/MWh'.format(stats.loc['market_index'].iloc[0]))
        
        plt.savefig(snakemake.output[f"{metric}_map"], bbox_inches='tight')
        plt.show()

    # pd.DataFrame().to_csv(snakemake.output["summary"])