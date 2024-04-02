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

    response = requests.get(accepts_url)
    accepts = pd.read_csv(StringIO(response.text))
    logger.info(f"Retrieved accted {len(accepts)} balancing actions.")

    bidsoffers_url = (
    "https://data.elexon.co.uk/bmrs/api/v1/balancing/settlement/"
    "indicative/volumes/all/{}/{}/{}?{}&format=json"
    )
    unit_insert = "bmUnit={}"
    accepted_units = accepts["NationalGridBmUnit"].unique().tolist()

    insertion = "&".join(unit_insert.format(unit) for unit in accepted_units)

    def get_trades(mode, date, period, insertion):
        response = requests.get(bidsoffers_url.format(mode, date, period, insertion))    
        return json_normalize(response.json()["data"])
    
    bids = get_trades("bid", date, period, insertion)
    offers = get_trades("offer", date, period, insertion)

    # data on all bid offers: Used to calculate the cost of in 
    trades_url = (
        "https://data.elexon.co.uk/bmrs/api/v1/balancing/bid-offer/" +
        "all?settlementDate={}&settlementPeriod={}&format=csv".format(date, period)
    )

    response = requests.get(trades_url)
    trades = pd.read_csv(StringIO(response.text))

    def get_balancing_summary(bm):
        accepts.loc[accepts["SoFlag"]]["NationalGridBmUnit"].value_counts()

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

    if len(so_units) > 0:
        so_actions = pd.concat((get_balancing_summary(bm) for bm in so_units), axis=1).T
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
