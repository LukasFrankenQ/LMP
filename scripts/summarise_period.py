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

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))

    stats = pd.read_csv(snakemake.input["price_stats"], index_col=0)

    axs[1].set_title('Wholesale Market Index: {}'.format(stats.loc['market_index'].iloc[0]))

    for layout, ax in zip(['nodal', 'fti', 'eso'], axs):

        n = pypsa.Network(snakemake.input["network_{}".format(layout)])
        regions = gpd.read_file(snakemake.input["regions_{}".format(layout)]).set_index("name")

        c = n.buses_t.marginal_price.columns
        regions['price'] = n.buses_t.marginal_price.iloc[0].loc[regions.index.intersection(c)]
        regions.plot(column='price', legend=True, ax=ax, vmin=0, vmax=100)

        ax.set_axis_off()

    plt.savefig(snakemake.output["maps"], bbox_inches='tight')
    pd.DataFrame().to_csv(snakemake.output["summary"])
    plt.show()

