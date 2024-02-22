# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
Uses FES 2021 data to assign load shares to different regions

**Outputs**

- ``RESOURCES/load_weights.csv``: share of load attributatable to each region

"""

import logging

from _helpers import configure_logging

import numpy as np
import pandas as pd
import geopandas as gpd

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    configure_logging(snakemake)

    scenario = snakemake.params["load"]["load_distribution_scenario"]    
    year = int(snakemake.params["load"]["load_distribution_year"])

    logger.info(f"Retrieving load distribution data for FES {year} scenario {scenario}.")

    assert year == 2021, "Only FES 2021 data is currently supported, change in cfg: 'load_distribution_year'"

    load_data = pd.read_csv(snakemake.input[f"fes_{year}_{scenario}_demandpeaks"])
    load_data = (
        load_data[["Grid Supply Points", f"2023 (MW)"]]
        .rename(columns={
            "Grid Supply Points": "GSP",
            "2023 (MW)": "load",})
        )
    weights = (l := load_data.set_index("GSP")).div(l.sum()).reset_index()
    
    bus_regions = gpd.read_file(snakemake.input["regions_onshore"])

    gsp_regions = (
        gpd.read_file(snakemake.input["gsp_regions"])
        .rename(
            columns={
                "RegionName": "region_name",
                "RegionID": "region_id",
            })
        .set_index("region_id")
        )

    gsp_lookup = pd.read_csv(snakemake.input["gsp_regions_lookup"]).set_index("region_id")

    def get_region_id(name):
        try:
            return int(gsp_lookup[gsp_lookup["gsp_name"] == name].index[0])
        except IndexError:
            return np.nan

    weights["region_id"] = weights["GSP"].apply(lambda x: get_region_id(x))
    weights["num_regions"] = weights["GSP"].apply(lambda x: len(x.split(";")))

    special = weights[weights["num_regions"] > 1]

    gsp_regions[["weight"]] = np.nan

    for region_id in gsp_regions.index:
        try:

            weights_idx = weights.loc[weights["region_id"] == region_id].index[0]
            gsp_regions.at[region_id, "weight"] = weights.loc[weights_idx, "load"]

        except IndexError:
            region_name = gsp_lookup.loc[region_id, "gsp_name"]
            if isinstance(region_name, pd.Series):
                region_name = region_name.values[0]

            hit = weights.loc[weights["GSP"].str.contains(region_name)].index

            if hit.empty:
                continue

            n_regions = weights.loc[hit[0], "num_regions"]
            gsp_regions.at[region_id, "weight"] = weights.loc[hit[0], "load"] / n_regions

    load_weights = pd.DataFrame(index=bus_regions.index, columns=["load_weight"])

    hold = gsp_regions.copy().to_crs("EPSG:4326")

    for i, row in bus_regions.iterrows():

        overlap = hold.geometry.intersection(row.geometry).area
        if overlap.sum() == 0.:
            continue

        overlap = overlap / overlap.sum()

        load_weights.at[i, "load_weight"] = overlap.mul(hold["weight"]).sum()

    load_weights = load_weights.fillna(0)
    load_weights = load_weights.div(load_weights.sum())

    (
        pd.concat((bus_regions, load_weights), axis=1)
        [["name", "load_weight"]]
        .set_index("name")
        .fillna(0.)
        .to_csv(snakemake.output["load_weights"])
    )

    logger.info(f"Saved load weights to {snakemake.output['load_weights']}.")







