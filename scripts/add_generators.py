# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
Assigns located generators to buses based on regions

**Outputs**

- ``RESOURCES/networks/gen.nc``: Network with generators assigned to buses

"""

import logging

import pypsa
import pandas as pd
import geopandas as gpd

from _helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    configure_logging(snakemake)

    logger.info("Assigning generators to buses based on regions.")

    n = pypsa.Network(snakemake.input["base_network"])

    bmus = pd.read_csv(snakemake.input["bmunits_loc"])
    bmus = bmus.loc[bmus.lat != 0.]

    bmus = (
        gpd.GeoDataFrame(bmus, geometry=gpd.points_from_xy(bmus.lon, bmus.lat))
        [["geometry", "NationalGridBmUnit", "capacity"]]
        ).set_crs(epsg=4326)

    onshore = gpd.read_file(snakemake.input["regions_onshore"]).set_index("name")
    offshore = gpd.read_file(snakemake.input["regions_offshore"]).set_index("name")

    non_assigned = 0

    for _, bmu in bmus.iterrows():

        hit = onshore.loc[onshore.geometry.contains(bmu.geometry)].index
        if len(hit) == 1:
            n.add("Generator", bmu["NationalGridBmUnit"], bus=hit[0], p_nom=bmu["capacity"])
            continue

        hit = offshore.loc[offshore.geometry.contains(bmu.geometry)].index
        if len(hit) == 1:
            n.add("Generator", bmu["NationalGridBmUnit"], bus=hit[0], p_nom=bmu["capacity"])
            continue

        non_assigned += 1

    logger.info(f"Assigned {len(n.generators)} generators to buses based on regions.")
    logger.info(f"Failed to assign {non_assigned} generators to buses based on regions.")
    logger.info(f"Saving network to {snakemake.output['gen_network']}.")

    n.export_to_netcdf(snakemake.output["gen_network"])