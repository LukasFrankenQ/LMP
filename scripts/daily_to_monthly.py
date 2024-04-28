# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT

import logging

logger = logging.getLogger(__name__)

import json
import numpy as np
import pandas as pd

from _helpers import configure_logging
from aggregate_periods import aggregate_stats

if __name__ == "__main__":

    configure_logging(snakemake)

    infiles = snakemake.input

    months = [fn.split('/')[-1].split('.')[0] for fn in infiles]
    layouts = ['nodal', 'national', 'eso', 'fti']

    date = snakemake.wildcards.year
    assert len(date.split('-')) == 1, 'Provided date should be year'

    logger.info(f"Aggregating to monthly resolution for year {date}, and months {months}.")

    months = [str(n).zfill(2) for n in range(1, 13)]
    monthly_results = {}

    for month in months:

        source_files = [file for file in infiles if month in file]
        if not source_files:
            logger.warning(f"No files found for {date}-{month}.")
            continue

        layout_dicts = {
            layout: {
                "geographies": {}
                } for layout in layouts
                }

        for layout in layouts:
            results = aggregate_stats(layout, source_files)

            for region in results.index:

                layout_dicts[layout]["geographies"][region] = {
                    "variables": results.T[region].fillna(0.).astype(np.float16).to_dict()
                }

        total_seconds = int(pd.Timestamp(f"{date}-{month}").timestamp())

        monthly_results[total_seconds] = layout_dicts

    with open(snakemake.output["monthly_aggregate"], "w") as f:
        json.dump(monthly_results, f)