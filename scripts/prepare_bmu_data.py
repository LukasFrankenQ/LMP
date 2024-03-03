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

import yaml
import requests
import numpy as np
import pandas as pd

from io import StringIO
from copy import deepcopy

from _helpers import configure_logging
from _elexon_helpers import process_multiples


power_hierarchy = {
    "wind": ["offwind", "onwind"],
    "gas": ["CCGT", "CHP"],
    "interconnector": ["HVDC submarine", ""],
    "hydro": ["hydro-scheme", "PHS", "dam"],
    "thermal": ["coal", "biomass"]
    }

def intersection(a, b):
    return list(set(a) & set(b))

def assign_carrier(cands):
    if "battery" in cands:
        return "battery"
    
    for outer, inner in power_hierarchy.items():
        if (i := intersection(cands, inner)):
            return i[0]

        if outer in cands:
            return outer
        
    return "powerstation"


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

    print("wiki df")
    print(wiki_df)
    print("--------------------------------------")

    carrier_mapper = yaml.safe_load(open(snakemake.input["carrier_mapper"]))
    print('carrier')
    print(carrier_mapper)

    wiki_df["carrier"] = wiki_df["instance"].apply(lambda entry: carrier_mapper[entry])

    print("--------------------------------------")
    print("wiki df")
    print(wiki_df)
    print("--------------------------------------")

    located_df = pd.concat((
        all_units,
        wiki_df[["bmrs_id", "lat", "lon", "capacity"]].groupby("bmrs_id").mean(),
    ), axis=1).fillna(0).loc[all_units.index]

    print("located_df")
    print(located_df)

    carriers = wiki_df["instance"].apply(lambda entry: carrier_mapper[entry])
    carriers.index = wiki_df["bmrs_id"]

    print("carrier")
    print(carriers)
    
    def get_carrier(entry):

        try:
            value = carriers.loc[entry]
        except KeyError:
            return np.nan

        if isinstance(value, pd.Series):
            return assign_carrier(value.tolist())

        elif isinstance(value, str):
            return value

        return "other"

    located_df["carrier"] = list(map(get_carrier, located_df.index))

    print("located_df")
    print(located_df.head())
    print(located_df.carrier.isna().mean())
    print(located_df["carrier"].tolist())
    print("bmrs_id" in located_df["carrier"].tolist())
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++===")

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

    print("=======================================")
    print('manual bmus')
    print(manual_bmus.head())
    print("bmrs_id" in manual_bmus["bmu_type"].values.tolist())
    print("++++++++++++++++++++++++++++++++++++++++++++++++++++++===")

    located_df.update(manual_bmus.rename(columns={"bmu_type": "carrier"})[["lon", "lat", "carrier"]])

    print("=======================================")
    print('final located')
    print(located_df.head(10))

    found_idx = located_df.loc[located_df["lat"] != 0].index
    logger.info(f"Found share {len(found_idx)/len(located_df)*100:.2f}% after wikidata, OSUKED and manual data.")

    located_df.to_csv(snakemake.output["bmunits_loc"])
    logger.info("Saved BMU data to " + snakemake.output["bmunits_loc"])