# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This rule gathers are accepted bids and offers from the balancing mechanism for
one settlement period. Only actions tagged by the So-Flag are considered - the flag
indicates the action being caused by transmission constraints

**Outputs**

- ``RESOURCES/{date}_{period}/real_balancing_actions.csv``: An overview of BMunits partaking in said balancing actions

"""

import logging

from _helpers import configure_logging
from _elexon_helpers import robust_request

import requests
import pandas as pd
from io import StringIO
from pandas import json_normalize

logger = logging.getLogger(__name__)

if __name__ == "__main__":

    configure_logging(snakemake)

    date = snakemake.wildcards.date    
    period = int(snakemake.wildcards.period)

    logger.info(f"Retrieving balancing actions for {date} settlement period {period}.")

    accepts_url = (
        'https://data.elexon.co.uk/bmrs/api/v1/balancing/acceptances/' +
        'all?settlementDate={}&settlementPeriod={}&format=csv'.format(date, period)
        )

    response = robust_request(requests.get, accepts_url)
    accepts = pd.read_csv(StringIO(response.text))
    logger.info(f"Retrieved {len(accepts)} balancing actions.")

    bidsoffers_url = (
    "https://data.elexon.co.uk/bmrs/api/v1/balancing/settlement/"
    "indicative/volumes/all/{}/{}/{}?{}&format=json"
    )
    unit_insert = "bmUnit={}"
    accepted_units = accepts["NationalGridBmUnit"].unique().tolist()
    
    print('accepted units')
    print(list(accepted_units))
    print('N accepted units')
    print(len(list(accepted_units)))

    so_units = accepts.loc[accepts.SoFlag]
    print('so accepted units')
    print(so_units[['NationalGridBmUnit']])
    print('N so accepted units')
    print(len(so_units))

    print('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
    

    # print('first ours then real')

    # accepted_units = ["PEMB-41", "GRAI-6"]

    insertion = "&".join(unit_insert.format(unit) for unit in accepted_units)

    # print(bidsoffers_url.format('offer', date, period, insertion))
    # url = 'https://data.elexon.co.uk/bmrs/api/v1/balancing/settlement/indicative/volumes/all/offer/2023-11-10/46?bmUnit=PEMB-41&bmUnit=GRAI-6&format=json'
    # print(url)

    # print(requests.get(url).json()['data'])

    bids_offers_url = "https://data.elexon.co.uk/bmrs/api/v1/balancing/settlement/indicative/volumes/all/{}/{}/{}?format=json"


    offers = json_normalize(
        # requests.get("https://data.elexon.co.uk/bmrs/api/v1/balancing/settlement/indicative/volumes/all/offer/2024-02-01/2?format=json").json()['data']
        robust_request(requests.get, bids_offers_url.format('offer', date, period)).json()['data']
    )
    bids = json_normalize(
        # requests.get("https://data.elexon.co.uk/bmrs/api/v1/balancing/settlement/indicative/volumes/all/bid/2024-02-01/2?format=json").json()['data']
        robust_request(requests.get, bids_offers_url.format('bid', date, period)).json()['data']
    )

    # print('all in there')
    # col = 'nationalGridBmUnit'
    # total = set(offers[col].tolist() + bids["nationalGridBmUnit"].tolist())

    # print(pd.Series(accepted_units).isin(total).all())

    print('bids')
    print(bids)
    print('offers')
    print(offers)

    # import sys
    # sys.exit()

    def get_trades(mode, date, period, insertion):
        print(f"Getting trades {mode}")
        print(bidsoffers_url.format(mode, date, period, insertion))

        response = requests.get(bidsoffers_url.format(mode, date, period, insertion))    
        print('response')
        print(response)
        print(response.json())
        return json_normalize(response.json()["data"])
    
    # bids = get_trades("bid", date, period, insertion)
    # offers = get_trades("offer", date, period, insertion)

    # print('bids')
    # print(bids)
    # print('offers')
    # print(offers)

    # data on all bid offers: Used to calculate the cost of in 
    trades_url = (
        "https://data.elexon.co.uk/bmrs/api/v1/balancing/bid-offer/" +
        "all?settlementDate={}&settlementPeriod={}&format=csv".format(date, period)
    )

    response = requests.get(trades_url)
    trades = pd.read_csv(StringIO(response.text))

    print(offers)

    def get_balancing_summary(bm):
        # accepts.loc[accepts["SoFlag"]]["NationalGridBmUnit"].value_counts()

        # print(f"====================== {bm} ======================")
        # drops = bids.columns[~bids.columns.str.contains("Volume")]
        # print(bids.loc[bids['nationalGridBmUnit'] == bm].drop(columns=drops))
        # print('------------------')
        # print(trades.loc[trades['NationalGridBmUnit'] == bm])
        # print("=============================================")

        if bm == 'DUNGW-1':
            print('================================')
            print(offers.loc[offers['nationalGridBmUnit'] == bm, ["totalVolumeAccepted"]])
            print('--------------------------------')
            print(trades.loc[trades['NationalGridBmUnit'] == bm, 'Offer'])
            print('--------------------------------')
            print(bids.loc[bids['nationalGridBmUnit'] == bm, ["totalVolumeAccepted"]])
            print('--------------------------------')
            print(trades.loc[trades['NationalGridBmUnit'] == bm, 'Bid'])
            print('================================')


        bm_summary = pd.Series({
            "offer volume": offers.loc[offers['nationalGridBmUnit'] == bm, ["totalVolumeAccepted"]].max().max(),
            "offer price": trades.loc[trades['NationalGridBmUnit'] == bm, 'Offer'].max(),
            "bid volume": bids.loc[bids['nationalGridBmUnit'] == bm, ["totalVolumeAccepted"]].abs().max().max(),
            "bid price": trades.loc[trades['NationalGridBmUnit'] == bm, 'Bid'].abs().max(),
            },
        name=bm)

        bm_summary['offer revenue'] = bm_summary['offer volume'] * bm_summary['offer price']
        bm_summary['bid revenue'] = bm_summary['bid volume'] * bm_summary['bid price']

        return bm_summary

    # The So Flag tracks which balancing actions are taken due to transmission constraints
    # Only these are considered here
    so_units = accepts.loc[accepts.SoFlag]["NationalGridBmUnit"].unique()


    print('offers')
    print(offers)
    print('bids')
    print(bids)

    if len(so_units) > 0:
        # for bm in so_units:
            # print('---------------------------------------------')
            # print(bm)
            # print(get_balancing_summary(bm))

        so_actions = pd.concat((get_balancing_summary(bm) for bm in so_units), axis=1).T
        print('final so actions')
        print(so_actions)
    else:
        # this is not nice...
        so_actions = pd.DataFrame(
            columns=[
                "offer volume",
                "offer price",
                "offer revenue",
                "bid volume",
                "bid price",
                "bid revenue"
                ]
            )

    so_actions.to_csv(snakemake.output["real_balancing_actions"])
