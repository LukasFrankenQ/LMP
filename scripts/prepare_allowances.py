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
        sheet_name='3b RO',
        header=8,
    ).iloc[-1, 8:].dropna()

    series.name = 'RO cost estimate'

    series.index = series.index.map(clean_date)
    return series.loc[pd.Timestamp(date)]


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
    return series.loc[pd.Timestamp(date)]


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
        return series.loc[date]
    except KeyError:
        return 0.


def get_eco(fn, date):
    """Obtains Energy Company Obligation (ECO) from OFGEM price cap Annex 4. Unit: £/MWh supplied."""

    series = pd.read_excel(
        fn,
        sheet_name='3e ECO',
        header=8,
        index_col=1,
    ).iloc[-1,6:].dropna()

    series.index = series.index.map(clean_date)
    return series.loc[pd.Timestamp(date)]


def get_whd(fn, date):
    """Obtains Warm Home Discount (WHD) from OFGEM price cap Annex 4. Unit: £/customer"""

    series = pd.read_excel(
        fn,
        sheet_name='3f WHD',
        header=8,
        index_col=1,
    ).iloc[-1,6:].dropna()

    series.index = series.index.map(clean_date)

    return series.loc[pd.Timestamp(date)]


def get_aahedc(fn, date, mode):
    """Obtains Assistance for Areas with High Electricity Distribution
    Costs (AAHEDC) from OFGEM price cap Annex 4. Unit: £/MWh supplied"""

    series = pd.read_excel(
        fn,
        sheet_name='3g AAHEDC',
        header=8,
        index_col=1,
    ).iloc[-1,6:].dropna()

    series.index = series.index.map(clean_date)
    flat_value = series.loc[date].values[0]

    losses = pd.read_excel(
        fn,
        sheet_name='3h Losses',
        header=11,
        index_col=[1,3],
    ).iloc[2:,6:].dropna(axis=1)

    losses.columns = losses.columns.map(clean_date)
    
    keys = {'multi': 'Multi-Register', 'single': 'Single Rate'}
    aa = (losses.loc[idx[keys[mode], :], pd.Timestamp(date)] * flat_value)

    aa.index = aa.index.get_level_values(1)

    return aa


def get_demand(fn, mode):
    """Obtains demand from OFGEM price cap Annex 4. Unit: MWh supplied"""

    df = pd.read_excel(
        fn,
        sheet_name='3a Demand',
    )

    if mode == 'single':
        return df.iloc[7, 2]
    elif mode == 'multi':
        return df.iloc[8, 2]


def get_weights(fn, mode):
    """Obtains share of electricity usages winter vs summer"""

    df = pd.read_excel(
        fn,
        sheet_name='3a Demand',
    )

    if mode == 'single':
        return {
            "summer": df.iloc[15, 2],
            "winter": df.iloc[15, 3],
        }
    elif mode == 'multi':
        return {
            "summer": df.iloc[16, 2],
            "winter": df.iloc[16, 3],
        }


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

    date = snakemake.wildcards["quarter"]

    logger.info(f"Preparing allowances for {date}.")

    nice_modes = {
        'single': 'Electricity - Single-Rate Metering Arrangement',
        'multi': 'Electricity - Multi-Register Metering Arrangement',
    }
    assert price_config['metering'] in nice_modes, f"Mode {price_config['mode']} not in {nice_modes.keys()}"
    
    logger.info("Gathering policy cost allowances: RO, FIT, ECO, WHD, AAHEDC.")

    regions.loc[:, 'renewable obligation'] = get_ro(policy_file, date)

    regions.loc[:, 'feed-in tariffs'] = get_fit(policy_file, date)

    regions.loc[:, 'old feed-in tariffs'] = get_old_fit(policy_file, date)

    regions.loc[:, 'energy company obligation'] = get_eco(policy_file, date)

    regions.loc[:, 'warm home discount'] = get_whd(policy_file, date)

    regions.loc[:, 'aahedc'] = get_aahedc(policy_file, date, price_config['metering'])

    regions.loc[:, 'total policy price check'] = get_policy_summary(policy_file, date, mode=price_config['metering'])

    per_delivered = [
        "renewable obligation",
        "feed-in tariffs",
        "old feed-in tariffs",
        "energy company obligation",
        "aahedc"
        ]

    demand = get_demand(policy_file, price_config['metering'])

    regions.loc[:, 'computed total'] = (
        regions[per_delivered].sum(axis=1) * demand 
        + regions['warm home discount']
    )

    logger.info("Gathering network cost allowances: TNUoS, DUoS, BSUoS.")

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
        
        regions = add_network_allowance(
            regions,
            consumption_df,
            quant + " " + consumptions[price_config["metering"]])
        

    logger.info("Gathering wholesale allowance: Direct Fuel, Backwardation, CFD.")

    def quarter_to_date(q):
        return q.split(' ')[1] + "-" + str(3 * (int(q.split(' ')[0][-1]) - 1) + 1).zfill(2)

    wholesale_slide = pd.read_excel(
        snakemake.input["wholesale_allowances"],
        sheet_name='1a Wholesale allowance',
        header=7,
        index_col=[1,2,3],
    )


    def get_direct_fuel(sheet, date, mode):
        """Takes Annex 2, 1a Wholesale Allowance, and returns direct fuel cost for requested date and mode."""

        modes = {
            "single": "Single-Rate Metering Arrangement",
            "multi": "Multi-Register Metering Arrangement",
        }

        df = sheet.copy().loc[idx["Electricity", modes[mode]]].iloc[:14,22:]

        df = df.loc[:, ~df.columns.str.contains("Unnamed")]
        df.columns = df.columns.map(quarter_to_date).map(pd.Timestamp)

        return df[pd.Timestamp(date)]


    def get_backwardation(sheet, date, mode):
        """Takes Annex 2, 1a Wholesale Allowance, and returns backwardation cost for requested date and mode."""

        modes = {
            "single": "Single-Rate Metering Arrangement",
            "multi": "Multi-Register Metering Arrangement",
        }

        df = sheet.copy().loc[idx["Electricity", modes[mode]]].iloc[14:28,22:]

        df = df.loc[:, ~df.columns.str.contains("Unnamed")]
        df.columns = df.columns.map(quarter_to_date).map(pd.Timestamp)

        return df[pd.Timestamp(date)]


    def get_cfd(sheet, date, mode):
        """Takes Annex 2, 1a Wholesale Allowance, and returns CFD cost for requested date and mode."""

        modes = {
            "single": "Single-Rate Metering Arrangement",
            "multi": "Multi-Register Metering Arrangement",
        }

        df = sheet.copy().loc[idx["Electricity", modes[mode]]].iloc[28:42,22:]

        df = df.loc[:, ~df.columns.str.contains("Unnamed")]
        df.columns = df.columns.map(quarter_to_date).map(pd.Timestamp)

        return df[pd.Timestamp(date)]
    

    regions.loc[:, 'direct fuel'] = get_direct_fuel(wholesale_slide, date, price_config['metering'])
    regions.loc[:, 'backwardation'] = get_backwardation(wholesale_slide, date, price_config['metering'])
    regions.loc[:, 'cfd'] = get_cfd(wholesale_slide, date, price_config['metering'])

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    pdf = regions.drop(columns=['geometry', 'total policy price check', 'computed total'])
    pos = pdf.loc[:, pdf.mean() > 0.]
    neg = pdf.loc[:, pdf.mean() <= 0.]

    pos.plot.bar(stacked=True, ax=ax, legend=True)
    neg.plot.bar(stacked=True, ax=ax, legend=True)

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(
        handles[::-1],
        labels[::-1],
        title='Cost Factor',
        bbox_to_anchor=(0.9, 1.1),
        ncols=5,
    )
    ax.set_ylabel('Cost per year per person (£)')
    ax.set_xlabel('Region')
    ax.grid(axis='y', linestyle='--', alpha=0.5)
    ax.xaxis.set_tick_params(rotation=45)

    plt.tight_layout()
    plt.show()

    demand = get_demand(policy_file, price_config['metering'])

    linear_cols = [
        "renewable obligation",
        "feed-in tariffs",
        "old feed-in tariffs",
        "energy company obligation",
        "aahedc",
        "direct fuel",
        "backwardation",
        "cfd",
        ] + [col for col in regions.columns if ' m ' in col]

    standing_cols = [
        "warm home discount"
    ] + [col for col in regions.columns if 'standing charge' in col]

    print('demand')
    print(demand)
    print('linear')
    print(linear_cols)
    print((regions[linear_cols].head() / demand).iloc[:,:len(linear_cols)//2])
    print((regions[linear_cols].head() / demand).iloc[:,len(linear_cols)//2:])

    print('standing')
    print(standing_cols)
    print(regions[standing_cols].head())

    regions[standing_cols].to_csv(snakemake.output["zeroth_order"])
    (regions[linear_cols] / demand).to_csv(snakemake.output["first_order"])
