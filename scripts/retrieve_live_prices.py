# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This rule downloads live data on BMU activity from the Elexon Insights API.

**Outputs**

- ``RESOURCES/{date}_{period}/price_stats.csv``: collected Physical Data

"""

import logging
import requests
import pandas as pd

from io import StringIO

logger = logging.getLogger(__name__)

from _helpers import configure_logging


def get_value(df, feature):
    """Helper function to clean received value."""
    for col in df.columns:
        if feature in col:
            feature = col.split(':')[-1]

            try:
                return float(feature)
            except ValueError:
                return float(feature[:-1])


url = "https://data.elexon.co.uk/bmrs/api/v1/balancing/pricing/market-index?from={}T00:00Z&to={}T00:00Z&settlementPeriodFrom={}&settlementPeriodTo={}".format(date, date, period, period)

if __name__ == "__main__":

    configure_logging(snakemake)

    date = snakemake.wildcards.date
    period = int(snakemake.wildcards.period)

    logger.info(f"Retrieving live wholesale price from Elexon Insights API for {date} settlement period {period}.")

    results = pd.Series()

    url = "https://data.elexon.co.uk/bmrs/api/v1/balancing/pricing/market-index?from={}T00:00Z&to={}T00:00Z&settlementPeriodFrom={}&settlementPeriodTo={}".format(date, date, period, period)
    response = requests.get(url)
    df = pd.read_csv(StringIO(response.text))

    results.loc['market_index'] = get_value(df, 'price')
    results.loc['volume'] = get_value(df, 'volume')
    
    results.to_csv(snakemake.output["price_stats"])