# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT

import logging

import json
from tqdm import tqdm

from _helpers import configure_logging


if __name__ == '__main__':

    configure_logging(snakemake)

    assert len(snakemake.input), (
        'Cannot create half-hourly periods without input data.'
        'Consider "hard" aggregation mode instead of "soft".'
    )
    logging.info(f'Gathering {len(snakemake.input)} periods in {snakemake.output[0]}.')

    halfhourly = {}

    for infile in tqdm(sorted(snakemake.input)):

        with open(infile) as f:
            data = json.load(f)
        
        halfhourly.update(data)
    
    with open(snakemake.output[0], 'w') as f:
        json.dump(halfhourly, f)
