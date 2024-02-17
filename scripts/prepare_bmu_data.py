# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This rule downloads a current list of Balancing Mechanism Units fromPhysical Data
Elexon Insights, and matches them with datasources to determine their location.

**Outputs**

- ``RESOURCES/bmunits_loc.csv``: collected data

"""

import logging

import requests
import pandas as pd

from io import StringIO
from copy import deepcopy

from _helpers import configure_logging
from _elexon_helpers import process_multiples

logger = logging.getLogger(__name__)

url = 'https://data.elexon.co.uk/bmrs/api/v1/balancing/physical/all'

if __name__ == "__main__":

    configure_logging(snakemake)

    logger.info(f"Retrieving list of all available BMUs from Elexon Insights from \n{url}")

    global query_params
    query_params = {
        'dataset': 'PN',  # Dataset to query
        'settlementDate': "2024-01-01",  # Settlement date in yyyy-MM-dd format
        'settlementPeriod': 1,  # Settlement period, an integer from 1-48
        'format': 'csv'  # Response data format
    }

    def get_bm_data(date, period):

        if not isinstance(date, str):
            date = date.strftime("%Y-%m-%d")
        
        params = deepcopy(query_params)
        params.update({"settlementDate": date, "settlementPeriod": period})

        response = requests.get(url, params=params)
        return pd.read_csv(StringIO(response.text))

    retrieval_date = snakemake.params["elexon"]["bmu_date"]

    df = get_bm_data(retrieval_date, 1)

    df = (
        process_multiples(df)
        .set_index(["SettlementDate", "SettlementPeriod", "NationalGridBmUnit"])
        [["LevelFrom", "LevelTo", "BmUnit"]]
    )

    print(df.head())
    print(df.tail())

    all_units = pd.DataFrame(index=df.index.get_level_values(2).unique())
    logger.info(f"Retrieved {len(all_units)} BMUs from Elexon Insights.")

    logger.info("Adding locations provided by wikidata.")

    wiki_df = pd.read_csv(snakemake.input["wiki_data"])

    print("wiki_data")
    print(wiki_df.head())

    merged = pd.concat((
        all_units,
        wiki_df[["bmrs_id", "lat", "lon", "capacity"]].groupby("bmrs_id").mean(),
    ), axis=1).fillna(0).loc[all_units.index]

    found_idx = merged.loc[merged["lat"] != 0].index

    logger.info(f"Found {len(found_idx)/len(merged)*100:.2f}% of BMUs in wikidata.")

    print("merged")
    print(merged)