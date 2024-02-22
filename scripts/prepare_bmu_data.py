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
import numpy as np
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

    all_units = pd.DataFrame(index=df.index.get_level_values(2).unique())
    logger.info(f"Retrieved {len(all_units)} BMUs from Elexon Insights.")

    logger.info("Adding locations provided by wikidata.")

    wiki_df = pd.read_csv(snakemake.input["wiki_data"])

    located_df = pd.concat((
        all_units,
        wiki_df[["bmrs_id", "lat", "lon", "capacity"]].groupby("bmrs_id").mean(),
    ), axis=1).fillna(0).loc[all_units.index]

    found_idx = located_df.loc[located_df["lat"] != 0].index
    logger.info(f"Found {len(found_idx)/len(located_df)*100:.2f}% of BMUs in wikidata.")

    logger.info("Adding locations provided by github.com/OSUKED/Power-Station-Dictionary.")

    osuked_ids = pd.read_csv(snakemake.input["osuked_ids"])
    osuked_plant_locations = pd.read_csv(snakemake.input["osuked_plant_locations"])

    for name in all_units.index:

        row = osuked_ids.loc[osuked_ids['ngc_bmu_id'].fillna('').str.contains(name)]

        if row.empty:
            continue

        bmu_id = row.index[0]

        try:
            located_df.loc[name, 'lat'] = osuked_plant_locations.loc[row.index[0], 'latitude'] or np.nan
            located_df.loc[name, 'lon'] = osuked_plant_locations.loc[row.index[0], 'longitude'] or np.nan
        except KeyError:
            continue

    found_idx = located_df.loc[located_df["lat"] != 0].index
    logger.info(f"Found share {len(found_idx)/len(located_df)*100:.2f}% after wikidata and OSUKED.")
    
    logger.info("Adding manual locations.")

    manual_bmus = pd.read_csv(snakemake.input["manual_bmus"]).set_index("nationalGridBmUnit")

    located_df.update(manual_bmus[["lon", "lat"]])

    found_idx = located_df.loc[located_df["lat"] != 0].index
    logger.info(f"Found share {len(found_idx)/len(located_df)*100:.2f}% after wikidata, OSUKED and manual data.")

    located_df.to_csv(snakemake.output["bmunits_loc"])
    logger.info("Saved BMU data to " + snakemake.output["bmunits_loc"])