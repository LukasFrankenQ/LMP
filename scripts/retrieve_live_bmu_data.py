# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This rule downloads live data on BMU activity from the Elexon Insights API.

**Outputs**

- ``RESOURCES/{date}_{period}/elexon_bmus.csv``: collected Physical Data

"""

import logging
import requests
import pandas as pd

from io import StringIO

logger = logging.getLogger(__name__)

from _helpers import configure_logging
from _elexon_helpers import process_multiples


prep_time = lambda x: str(x).zfill(2)

# url = 'https://data.elexon.co.uk/bmrs/api/v1/balancing/physical/all'
pn_url = "https://data.elexon.co.uk/bmrs/api/v1/datasets/PN?settlementDate={}&settlementPeriod={}&format=csv"
mels_url = "https://data.elexon.co.uk/bmrs/api/v1/datasets/MELS?from={}T{}%3A{}Z&to={}T{}%3A{}Z&format=csv"

if __name__ == "__main__":

    configure_logging(snakemake)

    date = snakemake.wildcards.date
    period = int(snakemake.wildcards.period)

    logger.info(f"Retrieving Live BMU Data from Elexon Insights API for {date} settlement period {period}.")

    response = requests.get(pn_url.format(date, period))
    pn = (
        process_multiples(
            pd.read_csv(
                StringIO(
                    response.text
                    )
                )
            )
        .set_index("NationalGridBmUnit")
        [["LevelTo"]]
        .rename(columns={"LevelTo": "PN"})
    )

    end = pd.Timestamp(date) + period * pd.Timedelta("30min")
    start = end - pd.Timedelta("30min")

    response = requests.get(
        mels_url.format(
            date,
            prep_time(start.hour),
            prep_time(start.minute),
            date,
            prep_time(end.hour),
            prep_time(end.minute)
        )
    )

    mels = (
        process_multiples(
            pd.read_csv(
                StringIO(response.text)
                )
            )
        .set_index("NationalGridBmUnit")
    )

    mels = (
        mels
        .rename(columns={"LevelTo": "MELS"})
        .sort_values(by="SettlementPeriod", ascending=True)
        .reset_index()
        .drop_duplicates(subset='NationalGridBmUnit', keep='first')
        .set_index("NationalGridBmUnit")
        ["MELS"]
    )

    pd.concat([pn, mels], axis=1).to_csv(snakemake.output["elexon_bmus"])