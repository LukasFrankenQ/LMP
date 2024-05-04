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

from _helpers import configure_logging, to_datetime, to_date_period
from _elexon_helpers import process_multiples, robust_request


prep_time = lambda x: str(x).zfill(2)

# url = 'https://data.elexon.co.uk/bmrs/api/v1/balancing/physical/all'
pn_url = "https://data.elexon.co.uk/bmrs/api/v1/datasets/PN?settlementDate={}&settlementPeriod={}&format=csv"
mels_url = "https://data.elexon.co.uk/bmrs/api/v1/datasets/MELS?from={}T{}%3A{}Z&to={}T{}%3A{}Z&format=csv"

if __name__ == "__main__":

    configure_logging(snakemake)

    date = snakemake.wildcards.date
    period = int(snakemake.wildcards.period)

    logger.info(f"Retrieving Live BMU Data from Elexon Insights API for {date} settlement period {period}.")

    max_tries = 3
    try_date, try_period = date, period

    for i in range(max_tries):

        response = robust_request(requests.get, pn_url.format(try_date, try_period))
        df = pd.read_csv(StringIO(response.text))

        if df.empty:
            try_date, try_period = to_date_period(
                to_datetime(try_date, try_period)
                - pd.Timedelta("30min")
            )

            logger.warning(f"Data unavailable; taking PN data for earlier {try_date} period {try_period}.")
            continue
        
        else:

            pn = (
                process_multiples(df)
                .set_index("NationalGridBmUnit")
                [["LevelTo"]]
                .rename(columns={"LevelTo": "PN"})
            )
            break

    # get Export Limit Data
    start = to_datetime(date, period)
    end = start + pd.Timedelta("30min")

    response = robust_request(
        requests.get,
        mels_url.format(
            start.strftime('%Y-%m-%d'),
            prep_time(start.hour),
            prep_time(start.minute),
            end.strftime('%Y-%m-%d'),
            prep_time(end.hour),
            prep_time(end.minute)
        )
    )

    df = pd.read_csv(StringIO(response.text))
    if not df.empty:
        mels = (
            process_multiples(df)
            .set_index("NationalGridBmUnit")
        )
    else:
        pn = pd.DataFrame(columns=["PN"])

    mels = (
        mels
        .rename(columns={"LevelTo": "MELS"})
        .sort_values(by="SettlementPeriod", ascending=True)
        .reset_index()
        .drop_duplicates(subset='NationalGridBmUnit', keep='first')
        .set_index("NationalGridBmUnit")
        ["MELS"]
    )

    output = pd.concat([pn, mels], axis=1).fillna(0.)
    pd.concat([pn, mels], axis=1).to_csv(snakemake.output["elexon_bmus"])