# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This rule processes wholesale and balancing modelling output and adds network
and policy costs, in total adding up to the end consumer price.

**Inputs**

- ``RESULTS/periods/{date}_{period}.json``: modelling output including wholesale and balancing prices
- ``data/price_postprocessing/Annex_3.xlsx``: network cost allowances
- ``data/price_postprocessing/Annex_4.xlsx``: policy cost allowances

**Outputs**

- ``RESULTS/periods_final/{date}_{period}.json``: modelling output with added variable 'consumer prices'

"""

import logging

import numpy as np
import pandas as pd
import geopandas as gpd

from _helpers import configure_logging

logger = logging.getLogger(__name__)

idx = pd.IndexSlice


def clean_date(s):
    s = s.replace('–','-')
    return pd.Timestamp(s.replace(' ','').split('-')[0])


def get_network_sheet(fn, sheet_name):
    """Processes Annex 3 to DataFrame."""

    df = pd.read_excel(
        fn,
        sheet_name=sheet_name,
        header=1,
        index_col=[1,2,3,5]
    )

    df.columns = df.iloc[5]
    df.index.names = ['metering', 'consumption', 'unit', 'region']
    df.columns.name = 'quarter'

    df = df.iloc[8:,4:]
    df = df.loc[:, ~df.columns.isna()]
    df.columns = df.columns.map(clean_date)

    return df


def get_ro(fn, date):
    """Obtains Renewable Obligation (RO) from OFGEM price cap Annex 4. Unit: £/MWh."""

    series = pd.read_excel(
        fn,
        sheet_name='3e ECO',
        header=8,
    ).iloc[-1, 8:].dropna()
    series.name = 'RO cost estimate'

    series.index = series.index.map(clean_date)
    return series.loc[date].values[0]


def get_fit(fn, date):
    """Obtains New Methodology Feed-In Tariffs (FIT) from OFGEM price cap Annex 4. Unit: £/MWh supplied."""

    series = pd.read_excel(
        fn,
        sheet_name='3i New FIT methodology',
        header=219,
        index_col=1,
    ).iloc[-1, 2:]

    series = series.loc[~series.index.str.contains('Unnamed')].fillna(0.)
    series.index = series.index.map(clean_date)

    return series.loc[date].values[0]


def get_old_fit(fn, date):
    """Obtains (Outdated) Feed-In Tariffs (FIT) from OFGEM price cap Annex 4. Unit: £/MWh supplied."""

    series = pd.read_excel(
        fn,
        sheet_name='3d FIT',
        header=7,
        index_col=1,
    ).iloc[-1, 6:]

    series = series.loc[~series.index.str.contains('Unnamed')].fillna(0.)
    series.index = series.index.map(clean_date)
    
    try:
        return series.loc[date].values[0]
    except KeyError:
        return 0.


def get_eco(fn, date):
    """Obtains Energy Company Obligation (ECO) from OFGEM price cap Annex 4. Unit: £/MWh supplied."""

    series = pd.read_excel(
        fn,
        sheet_name='3e ECO',
        header=7,
        index_col=1,
    ).iloc[-1,6:].dropna()

    series.index = series.index.map(clean_date)

    return series.loc[date].values[0]


def get_whd(fn, date):
    """Obtains Warm Home Discount (WHD) from OFGEM price cap Annex 4. Unit: £/customer"""

    series = pd.read_excel(
        fn,
        sheet_name='3f WHD',
        header=7,
        index_col=1,
    ).iloc[-1,6:].dropna()

    series.index = series.index.map(clean_date)

    return series.loc[date].values[0]


def get_aahedc(fn, date):
    """Obtains Assistance for Areas with High Electricity Distribution
    Costs (AAHEDC) from OFGEM price cap Annex 4. Unit: £/MWh supplied"""

    series = pd.read_excel(
        fn,
        sheet_name='3g AAHEDC',
        header=7,
        index_col=1,
    ).iloc[-1,6:].dropna()

    series.index = series.index.map(clean_date)

    return series.loc[date].values[0]


def get_policy_summary(fn, date, mode='multi'):
    """Obtains a summary of policy cost allowance from OFGEM price cap Annex 4. Unit: £/MWh supplied"""

    assert mode in ['multi', 'single']
    nice_modes = {
        'single': 'Electricity - Single-Rate Metering Arrangement',
        'multi': 'Electricity - Multi-Register Metering Arrangement',
    }

    df = pd.read_excel(
        fn,
        sheet_name='1a Policy Cost Allowance',
        header=11,
        index_col=[1,2],
    ).iloc[2:30, 4:]

    df = df.loc[idx[nice_modes[mode]], ~df.columns.str.contains('Unnamed')]
    df.columns = df.columns.map(clean_date)

    return df[pd.Timestamp(date)]


if __name__ == "__main__":

    configure_logging(snakemake)    

    policy_file = snakemake.input["policy_allowances"]
    network_file = snakemake.input["network_allowances"]

    regions = gpd.read_file(snakemake.input["dno_regions"]).set_index('LongName')

    price_config = snakemake.params["consumer_price"]

    print('imported regions')
    print(regions.head())

    date, period = pd.Timestamp(snakemake.wildcards["date"]), int(snakemake.wildcards["period"])

    nice_modes = {
        'single': 'Electricity - Single-Rate Metering Arrangement',
        'multi': 'Electricity - Multi-Register Metering Arrangement',
    }
    assert price_config['metering'] in nice_modes, f"Mode {price_config['mode']} not in {nice_modes.keys()}"

    logger.info("Retrieving network cost allowances: TNUoS, DUoS, BSUoS.")
    tnuos, duos, bsuos = (
        get_network_sheet(network_file, sheet_name)
        for sheet_name in ['2a TNUoS', '2b DUoS', '2c BSUoS']
        )

    def add_network_allowance(regions, pdf, colname):

        pdf.index = pdf.index.get_level_values('region')
        pdf = pdf.rename(columns={pdf.columns[0]: colname})
        pdf = pd.concat([regions, pdf], axis=1)
        pdf[colname] = pdf[colname].astype(float)

        return gpd.GeoDataFrame(pdf)

    for quant in ['tnuos', 'duos', 'bsuos']:

        df = globals()[quant].copy()

        modes = {
            "multi": "Multi-Register Metering Arrangement",
            "single": "Single-Rate Metering Arrangement",
        }
        consumptions = {
            "multi": "m (4,200 kWh)",
            "single": "m (3,100 kWh)",
        }

        df = df.loc[:, :date].iloc[:,[-1]]

        try:
            standing_df = df.copy().loc[
                idx[modes[price_config["metering"]], "Nil"],
                ]
        except KeyError:
            # No standing charge for BSUoS, effectively replacing lack of values with 0
            standing_df = (
                tnuos.copy().loc[
                    idx[modes[price_config["metering"]], "Nil"], :date
                ]
                .iloc[:,[-1]]
            )

            standing_df.values[:] = 0.

        regions = add_network_allowance(regions, standing_df, quant + ' standing charge')

        consumption_df = df.copy().loc[
            idx[modes[price_config["metering"]], consumptions[price_config["metering"]]],
            ]
        
        regions = add_network_allowance(regions, consumption_df, quant + consumptions[price_config["metering"]])
        
    print('final regions with all network costs')
    print(regions.head())

