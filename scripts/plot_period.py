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

import json
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable


logger = logging.getLogger(__name__)

from _helpers import configure_logging, to_datetime

if __name__ == "__main__":

    configure_logging(snakemake)

    layouts = ['national', 'nodal', 'fti', 'eso']

    with open(snakemake.input["data"]) as f:    
        data = json.load(f)
    data = data[list(data)[0]]

    vmin = 0
    vmax = 120

    fig, axs = plt.subplots(1, 4, figsize=(12, 5))

    for i, (layout, ax) in enumerate(zip(layouts, axs)):

        method = 'capitalize' if len(layout) > 3 else 'upper'
        ax.set_title(f"{getattr(layout, method)()} layout")

        inner_data = data[layout]['geographies']

        regions = gpd.read_file(snakemake.input["regions_{}".format(layout)]).set_index("name")
        
        get_price = lambda x: inner_data[x]["variables"]["post_balancing_price"]
        regions["price"] = list(map(get_price, regions.index))

        regions.plot(
            column="price",
            legend=False,
            ax=ax,
            vmin=vmin,
            vmax=vmax,
            edgecolor='k',
            linewidth=0.5,
            legend_kwds={'label': "Price [£/MWh]"}
            )
        
        ax1_divider = make_axes_locatable(ax)
        if i == 3:
            from matplotlib import colors

            norm = colors.Normalize(vmin=vmin, vmax=vmax)
            cbar = plt.cm.ScalarMappable(norm=norm, cmap='viridis')

            cax1 = ax1_divider.append_axes("right", size="7%", pad="2%")
            cb1 = fig.colorbar(cbar, cax=cax1)
            cax1.set_ylabel('Price [£/MWh]')#, rotation=270, labelpad=15)


        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(-8, 2)
        ax.set_ylim(49, 61)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)

    axs[0].text(
        -9, 48,
        "Plotting {}, ".format(to_datetime(snakemake.wildcards.date, snakemake.wildcards.period)),
        fontsize=8)

    plt.savefig(snakemake.output["plot"], bbox_inches='tight')
    plt.show()
