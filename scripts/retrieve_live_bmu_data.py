# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This rule downloads live data on BMU activity from the Elexon Insights API.

**Outputs**

- ``RESOURCES/{date}_{period}/bmu_physical.csv``: collected Physical Data

"""

import logging
import requests
import pandas as pd

from io import StringIO

logger = logging.getLogger(__name__)

from _helpers import configure_logging
from _elexon_helpers import process_multiples

url = 'https://data.elexon.co.uk/bmrs/api/v1/balancing/physical/all'

if __name__ == "__main__":

    configure_logging(snakemake)

    date = snakemake.wildcards.date
    period = snakemake.wildcards.period

    logger.info(f"Retrieving Live BMU Data from Elexon Insights API for {date} settlement period {period}.")

    params = {
        'dataset': 'PN',
        'settlementDate': date,
        'settlementPeriod': period,
        'format': 'csv'  # Response data format
    }
        
    response = requests.get(url, params=params)
    df = pd.read_csv(StringIO(response.text))

    df = (
        process_multiples(df)
        .set_index(["SettlementDate", "SettlementPeriod", "NationalGridBmUnit"])
        [["LevelFrom", "LevelTo", "BmUnit"]]
    ).to_csv(snakemake.output["bmu_physical_data"])