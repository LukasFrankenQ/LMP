# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This script assigns BMUs to their respective strike prices using two datasets.

Dataset on UK strike prices is from the UK government and is available at:
https://cfd.lowcarboncontracts.uk/

And the respective mapping to BMUs:
https://dp.lowcarboncontracts.uk/dataset/cfd-to-bm-unit-mapping

**Outputs**

- ``RESOURCES/bmu_strike_prices.csv``
"""

import logging

import pandas as pd

logger = logging.getLogger(__name__)

from _helpers import configure_logging

if __name__ == "__main__":

    configure_logging(snakemake)

    logger.info("Mapping BMUs to strike prices.")

    # d = pd.read_csv(path / 'cfd' / 'CfD Register 6-Jun-2024.csv')
    reg = pd.read_excel(
        snakemake.input['cfd_register'], index_col=0
        ).rename(columns={
            'Current strike price (field_cfd_current_strikeprice)': 'strike_price',
            'Technology type (field_cfd_technology_type)': 'carrier',
            'Unique Identifier (field_cfd_unique_id)': 'CFD_Id',
            })

    reg = (
        reg[['carrier', 'strike_price', 'CFD_Id']]
        .dropna()
        .set_index('CFD_Id')
    )

    mapper = pd.read_csv(
        snakemake.input['cfd_bmu_mapping'],
        index_col=1
        )

    mapper['strike_price'] = mapper['CFD_Id'].apply(lambda x: reg.loc[x, 'strike_price'])
    mapper['carrier'] = mapper['CFD_Id'].apply(lambda x: reg.loc[x, 'carrier'])
    mapper = mapper[['CFD_Id', 'strike_price', 'carrier']]

    mapper.index = mapper.index.map(lambda x: x.split('_')[-1])
    mapper = mapper[~mapper.index.duplicated(keep='first')]

    mapper.to_csv(snakemake.output["strike_prices"])
