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

from _helpers import configure_logging

import requests
import pandas as pd

logger = logging.getLogger(__name__)

if __name__ == "__main__":

    date = snakemake.wildcards.date
    period = int(snakemake.wildcards.period)

    logger.info(f"Retrieving constraint flows for {date} settlement period {period}.")

    start = pd.Timestamp(date) + pd.Timedelta(minutes=(period - 1) * 30)
    end = pd.Timestamp(date) + pd.Timedelta(minutes=period * 30) - pd.Timedelta(microseconds=1)

    sql_query = (
        '''SELECT COUNT(*) OVER () AS _count, * FROM "38a18ec1-9e40-465d-93fb-301e80fd1352"'''+
        ''' WHERE "Date (GMT/BST)" >= '{}' '''.format(start) +
        '''AND "Date (GMT/BST)" <= '{}' '''.format(end) +
        '''ORDER BY "_id" ASC LIMIT 1000'''
    )

    params = {'sql': sql_query}

    try:
        resposne = requests.get('https://api.nationalgrideso.com/api/3/action/datastore_search_sql', params = parse.urlencode(params))
        data = resposne.json()["result"]

        df = (
            pd.DataFrame(data["records"])
            .set_index("Constraint Group")
            [["_count", "Limit (MW)", "Flow (MW)", "Date (GMT/BST)"]]
            .rename(columns={
                "Limit (MW)": "limit",
                "Flow (MW)": "flow",
                "Date (GMT/BST)": "date",
                "_count": "period",
                })
        )
    except requests.exceptions.RequestException as e:
        print(e.response.text)
    
    df.to_csv(snakemake.output["constraint_flows"])