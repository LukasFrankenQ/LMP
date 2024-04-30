# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT

import logging

import json
from tqdm import tqdm

from _helpers import configure_logging
from _aggregation_helpers import aggregate_stats


if __name__ == '__main__':

    configure_logging(snakemake)

    logging.info(f'Aggregating {len(snakemake.input)} periods to monthly resolution in {snakemake.output[0]}.')
    monthly = {}

    for infile in tqdm(sorted(snakemake.input)):

        with open(infile) as f:
            daily_data = json.load(f)

        monthly.update(
            {list(daily_data)[0]: aggregate_stats(daily_data)}
            )
    
    with open(snakemake.output[0], 'w') as f:
        json.dump(monthly, f)