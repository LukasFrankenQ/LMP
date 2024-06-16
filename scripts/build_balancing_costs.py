# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
From accepted bids and offers that are soFlagged computes the daily weighted average

! THIS FUNCTION IS CURRENTLY NOT EMBEDDED IN THE WORKFLOW AS THE INPUT FILES ARE TOO LARGE !

**Inputs**

- ``bmrs-data-all/brms_detsysprices_weekly``

**Outputs**

- ``RESOURCES/daily_balancing_costs.py``

"""

import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from pathlib import Path

from _helpers import to_datetime


def get_period_prices(df, mode):

    result = pd.Series(index=df.index.unique())

    for date, period in tqdm(df.index.unique()):

        result.loc[idx[date, period]] = np.average(
            df.loc[idx[date, period], mode+'price'],
            weights=df.loc[idx[date, period], mode+'volume']
            )
    
    return result


path = Path.cwd().parent / 'data' / 'bmrs-data-all' / 'bmrs_detsysprices_weekly'

if __name__ == '__main__':

    idx = pd.IndexSlice

    print("Calculating daily balancing costs between Jan 2022 and May 2024.")

    all_offers = list()
    all_bids = list()

    for i, fn in enumerate(os.listdir(path)):

        if not ('2022' in fn or '2023' in fn or '2024' in fn):
            continue

        cols = [
            'volume',
            'settlementDate',
            'settlementPeriod',
            'soFlag',
            'storFlag',
            'bidPrice',
            'bidVolume',
            'offerVolume',
            'offerPrice',
            'acceptanceId',
        ]
        df = pd.read_parquet(path / fn)# [cols]
        df = df.replace(to_replace=['F'], value=False).replace(to_replace=['T'], value=True)
        df = df.loc[df[['soFlag', 'storFlag']].any(axis=1)]
        df = df.loc[df['acceptanceId'] != 'NULL']

        df.index = df.apply(
            lambda x: to_datetime(x['settlementDate'], x['settlementPeriod']),
            axis=1
            )

        bids = df.loc[df.bidVolume.astype(float) < 0]
        offers = df.loc[df.offerVolume.astype(float) > 0]

        all_offers.append(offers[cols])
        all_bids.append(bids[cols])

    o = pd.concat(all_offers)
    o = o.loc[o['soFlag']]

    o = o.set_index(['settlementDate', 'settlementPeriod'])[['offerVolume', 'offerPrice']]

    o['offerprice'] = o['offerPrice'].astype(float)
    o['offervolume'] = o['offerVolume'].astype(float)
    o.loc[~o['offerPrice'].isin([999., 9999.])]

    o = get_period_prices(o, 'offer')

    b = pd.concat(all_bids)

    b = b.loc[b['soFlag']]

    b = b.set_index(['settlementDate', 'settlementPeriod'])[['bidVolume', 'bidPrice']]

    b['bidprice'] = b['bidPrice'].astype(float)
    b['bidprice'] = b['bidprice'].abs()
    b.loc[~b['bidPrice'].isin([999., 9999.])]

    b['bidvolume'] = b['bidVolume'].astype(float)
    b['bidvolume'] = b['bidvolume'].abs()

    b = get_period_prices(b, 'bid')

    o.index = list(map(to_datetime, o.index.get_level_values(0), o.index.get_level_values(1)))
    b.index = list(map(to_datetime, b.index.get_level_values(0), b.index.get_level_values(1)))

    freq = 'd'
    o = o.resample(freq).mean().interpolate(method='linear').rename('offers')
    b = b.resample(freq).mean().interpolate(method='linear').rename('bids')

    pd.concat((o, b), axis=1).to_csv(
        Path.cwd().parent /
        'data' /
        f'{freq}_balancing_cost_2022_2023_2024.csv'
        )
