# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This rule gathers from national grid ESO per main boundary the flow limit for a given
settlement period.

**Outputs**

- ``RESOURCES + "live_data/{date}_{period}/constraint_flows.csv"``
"""

import logging

from _helpers import configure_logging, to_datetime, to_date_period
from _elexon_helpers import robust_request

import requests
import pandas as pd
from urllib import parse
from pathlib import Path

logger = logging.getLogger(__name__)


def retrieve_constraints(date, period):

    end = to_datetime(date, period) - pd.Timedelta(microseconds=1)
    start = end - pd.Timedelta(minutes=30) + pd.Timedelta(microseconds=1)

    sql_query = (
        '''SELECT COUNT(*) OVER () AS _count, * FROM "38a18ec1-9e40-465d-93fb-301e80fd1352"'''+
        ''' WHERE "Date (GMT/BST)" >= '{}' '''.format(start) +
        '''AND "Date (GMT/BST)" <= '{}' '''.format(end) +
        '''ORDER BY "_id" ASC LIMIT 1000'''
    )

    params = {'sql': sql_query}

    response = robust_request(
        requests.get,
        'https://api.nationalgrideso.com/api/3/action/datastore_search_sql',
        params=parse.urlencode(params)
        )

    data = response.json()["result"]

    df = (
        pd.DataFrame(data["records"])
        .set_index("Constraint Group")
        [["_count", "Limit (MW)", "Flow (MW)", "Date (GMT/BST)"]]
        .rename(columns={
            "Limit (MW)": "limit",
            "Flow (MW)": "flow",
            "Date (GMT/BST)": "date",
            })
    )

    return df


if __name__ == "__main__":

    configure_logging(snakemake)

    date = snakemake.wildcards.date
    period = int(snakemake.wildcards.period)

    logger.info(f"Retrieving constraint flows for {date} settlement period {period}.")

    max_tries = 400

    try_date, try_period = date, period
    df = pd.DataFrame()

    def is_valid(df):
        """Makes sure that obtained flows are neither all zero nor empty"""
        return not (df.empty or df.limit.sum() == 0)

    for i in range(max_tries):

        try:
            df = retrieve_constraints(try_date, try_period)
            method = "Live"

        except (requests.exceptions.RequestException, KeyError) as e:
            print(f'Detected error, trying previous period: {e}.')

            if (fn := Path(snakemake.params.RESOURCES +
                f"/live_data/{try_date}_{try_period}/constraint_flows.csv")
                ).exists():

                df = pd.read_csv(fn, index_col=0)
                method = "Past"

        if not is_valid(df):

            try_date, try_period = to_date_period(
                to_datetime(try_date, try_period)
                - pd.Timedelta("30min")
            )
            continue 

        else:
            break

    if i == 0:
        method = "Live"

    logger.info(f"{method} constraint data for {date}, {period}, taken from {try_date}, {try_period}.")
    df.to_csv(snakemake.output["constraint_flows"])