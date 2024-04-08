# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
Summarized data for individual periods to a coarser time resolution.
Simple averages are calculated for loads, dispatch and capacity. For marginal_price, 
the scripts computes a weighted mean based on load

**Outputs**

- ``RESULTS/[daily/monthly]/[time-key].json

"""

import json
import logging
import numpy as np
import pandas as pd

from _helpers import configure_logging

logger = logging.getLogger(__name__)


def get_stat(fn, layout, variable):

    with open(fn) as f:
        data = json.load(f)

    ts = list(data)[0]
    geos = list(data[ts][layout]['geographies'])
    return pd.Series({geo: data[ts][layout]["geographies"][geo]['variables'][variable] for geo in geos})


if __name__ == "__main__":

    configure_logging(snakemake)

    layout_dicts = {layout: {"geographies": {}} for layout in ['nodal', 'national', 'eso', 'fti']}
    filelist = list(snakemake.input)

    outfile = snakemake.output[0]
    if 'half-hourly' in outfile:
        logger.info("Aggregating to half-hourly resolution by merging into single file.")

        total_dict = {}

        for fn in filelist:
            with open(fn) as f:
                data = json.load(f)

            ts = list(data)[0]
            total_dict[ts] = data[ts]

        with open(outfile, "w") as f:
            json.dump(total_dict, f)

    elif "daily" in outfile:

        for layout in list(layout_dicts):

            prices = [get_stat(fn, layout, 'marginal_price') for fn in filelist]
            prices = pd.concat(prices, axis=1)

            loads = [get_stat(fn, layout, 'load') for fn in filelist]
            loads = pd.concat(loads, axis=1)

            capacity = [get_stat(fn, layout, 'available_capacity') for fn in filelist]
            capacity = pd.concat(capacity, axis=1)

            dispatch = [get_stat(fn, layout, 'dispatch') for fn in filelist]
            dispatch = pd.concat(dispatch, axis=1)

            normalized_loads = loads.sum() / loads.sum().sum() * loads.shape[1]
            avg_prices = prices.mul(normalized_loads).mean(axis=1)

            results = pd.concat([
                avg_prices,
                loads.mean(axis=1),
                capacity.mean(axis=1),
                dispatch.mean(axis=1)
                ], axis=1)

            results.columns = ['marginal_price', 'load', 'available_capacity', 'dispatch']

            for region in results.index:

                layout_dicts[layout]["geographies"][region] = {
                    "variables": results.T[region].fillna(0.).astype(np.float16).to_dict()
                }

        total_seconds = int(pd.Timestamp(snakemake.params.date).timestamp())

        with open(snakemake.output[0], "w") as f:
            json.dump({total_seconds: layout_dicts}, f)
