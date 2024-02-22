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

from _helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":

    configure_logging(snakemake)

    date = snakemake.wildcards.date    
    period = snakemake.wildcards.period

    logger.info(f"Preparing live network for {date} settlement period {period}.")
    n = pypsa.Network(snakemake.input["network"])

    bmu = pd.read_csv(snakemake.input["bmu_physical_data"]).set_index("NationalGridBmUnit")["LevelTo"]

    load_weights = pd.read_csv(snakemake.input["load_weights"], index_col=0)["load_weight"]

    total_load = bmu.loc[bmu > 0].sum()

    loads = load_weights.loc[load_weights > 0].index
    n.madd(
        "Load",
        loads,
        bus=loads,
        p_set=load_weights.loc[loads] * total_load,
    )

    bmu = bmu.loc[n.generators.index]

    logger.warning("Just deleting negative ones")
    bmu = bmu.loc[bmu > 0]

    n.generators.loc[bmu.index, 'p_nom'] = bmu

    logger.warning("no sensible costs yet!")

    n.export_to_netcdf(snakemake.output["live_network"])
