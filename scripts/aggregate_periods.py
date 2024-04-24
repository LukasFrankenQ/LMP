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
from tqdm import tqdm

from _helpers import configure_logging, get_datelist

logger = logging.getLogger(__name__)


def get_stat(fn, layout, variable):

    with open(fn) as f:
        data = json.load(f)

    ts = list(data)[0]
    geos = list(data[ts][layout]['geographies'])
    return pd.Series({geo: data[ts][layout]["geographies"][geo]['variables'][variable] for geo in geos})


if __name__ == "__main__":

    configure_logging(snakemake)

    infiles = list(snakemake.input)
    outfiles = snakemake.output

    for outfile in outfiles:

        if 'half-hourly' in outfile:

            layout_dicts = {layout: {"geographies": {}} for layout in ['nodal', 'national', 'eso', 'fti']}

            logger.info(f"Aggregating to half-hourly resolution to {outfile}.")

            total_dict = {}

            date = outfile.split('/')[-1].split('.')[0]

            for fn in [file for file in infiles if date in file]:
                with open(fn) as f:
                    data = json.load(f)

                ts = list(data)[0]
                total_dict[ts] = data[ts]

            with open(outfile, "w") as f:
                json.dump(total_dict, f)

        elif "daily" in outfile:

            daily_results = {}

            # print('filelist')
            # print(filelist)

            for date in tqdm(get_datelist(snakemake.params.date)):

                layout_dicts = {layout: {"geographies": {}} for layout in ['nodal', 'national', 'eso', 'fti']}
                date_files = [fn for fn in infiles if date in fn]

                for layout in list(layout_dicts):

                    whole_prices = [get_stat(fn, layout, 'wholesale_price') for fn in date_files]
                    whole_prices = pd.concat(whole_prices, axis=1)

                    pb_prices = [get_stat(fn, layout, 'post_balancing_price') for fn in date_files]
                    pb_prices = pd.concat(pb_prices, axis=1)

                    loads = [get_stat(fn, layout, 'load') for fn in date_files]
                    loads = pd.concat(loads, axis=1)

                    capacity = [get_stat(fn, layout, 'available_capacity') for fn in date_files]
                    capacity = pd.concat(capacity, axis=1)

                    dispatch = [get_stat(fn, layout, 'dispatch') for fn in date_files]
                    dispatch = pd.concat(dispatch, axis=1)

                    normalized_loads = loads.sum() / loads.sum().sum() * loads.shape[1]

                    results = pd.concat([
                        whole_prices.mul(normalized_loads).mean(axis=1),
                        pb_prices.mul(normalized_loads).mean(axis=1),
                        loads.mean(axis=1),
                        capacity.mean(axis=1),
                        dispatch.mean(axis=1)
                        ], axis=1)

                    results.columns = [
                        'wholesale_price',
                        'post_balancing_price',
                        'load',
                        'available_capacity',
                        'dispatch',
                        ]

                    for region in results.index:

                        layout_dicts[layout]["geographies"][region] = {
                            "variables": results.T[region].fillna(0.).astype(np.float16).to_dict()
                        }

                # total_seconds = int(pd.Timestamp(snakemake.params.date).timestamp())
                total_seconds = int(pd.Timestamp(date).timestamp())

                daily_results[total_seconds] = layout_dicts

            # with open(snakemake.output[0], "w") as f:
            with open(outfile, "w") as f:
                json.dump(daily_results, f)
