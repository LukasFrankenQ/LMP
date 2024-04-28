# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
Procomputes the overlap between regions used for modelling and the DNO regions to which
inhomogeneous costs are assigned in postprocessing. 

**Outputs**
- ``RESOURCES/overlap_[layout].csv

"""

import logging

logger = logging.getLogger(__name__)

import pypsa
import numpy as np
import pandas as pd
import geopandas as gpd

from _helpers import configure_logging


if __name__ == "__main__":
    configure_logging(snakemake)

    dno_regions = gpd.read_file(snakemake.input["dno_regions"]).set_index('LongName')["geometry"]

    for layout in ['eso', 'fti', 'nodal']:
        regions = gpd.read_file(snakemake.input[f"{layout}_regions"]).set_index('name')["geometry"]
        load_weights = (l := pypsa.Network(snakemake.input[f"{layout}_network"]).loads.p_set) / l.sum()

        overlap = pd.DataFrame(columns=regions.index, index=dno_regions.index, data=np.nan)

        for region, geom in regions.items():

            overlap[region] = dno_regions.intersection(regions.loc[region].buffer(0.01)).area
            # overlap[region] *= overlap[region].mul(load_weights) # this goes purely by area, should also go by load/population?
            overlap[region] /= overlap[region].sum()

        overlap = overlap.fillna(0.)

        overlap.to_csv(snakemake.output[layout])
        logger.info(f"Computed overlap between DNO regions and {layout} layout.")

