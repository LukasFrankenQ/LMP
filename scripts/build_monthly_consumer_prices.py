# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT

import logging

logger = logging.getLogger(__name__)

import json
import numpy as np
import pandas as pd

from _helpers import configure_logging

if __name__ == "__main__":

    configure_logging(snakemake)

    logger.info(f"Processing results to estimate monthly consumer prices.")

    year = snakemake.wildcards.year

    with open(snakemake.input["monthly_results"]) as f:
        results = json.load(f)

    print(snakemake.input)
    for ts, date in results.item():
        print(pd.Timestamp(ts).year, pd.Timestamp(ts).month)

        break

    



