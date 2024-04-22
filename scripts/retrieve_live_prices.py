# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This rule downloads live data on wholesale market prices from the Elexon Insights API.

**Outputs**

- ``RESOURCES/{date}_{period}/price_stats.csv``: collected Physical Data

"""

import logging

import requests
import pandas as pd

from io import StringIO

logger = logging.getLogger(__name__)

from _helpers import configure_logging

template = "https://data.elexon.co.uk/bmrs/api/v1/balancing/pricing/market-index?from={}T00:00Z&to={}T00:00Z&settlementPeriodFrom={}&settlementPeriodTo={}"

if __name__ == "__main__":

    configure_logging(snakemake)

    date = snakemake.wildcards.date
    period = int(snakemake.wildcards.period)

    logger.info(f"Retrieving live wholesale price from Elexon Insights API for {date} settlement period {period}.")

    results = pd.Series()

    url = template.format(date, date, period, period)
    response = requests.get(url)
    df = pd.read_csv(StringIO(response.text))

    for col1, col2 in zip(df.columns[:-1], df.columns[1:]):
        if not ('price' in col1 and 'volume' in col2):
            continue

        results.loc['market_index'] = float(col1.split(':')[-1])
        break
    
    results.to_csv(snakemake.output["price_stats"])