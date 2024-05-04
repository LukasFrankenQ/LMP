# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT


import time
import numpy as np
import pandas as pd
from requests.exceptions import HTTPError


def robust_request(getfunc, *args, max_retries=100, wait_time=5, **kwargs):

    for _ in range(max_retries):
        try:
            response = getfunc(*args, **kwargs)
            return response
        except HTTPError as http_err:
            print(f"HTTP error occurred: {http_err}")
        except Exception as err:
            print(f"Other error occurred: {err}")
        time.sleep(wait_time)

    raise Exception(f"Failed to get response after {max_retries} retries")


def process_multiples(df):
    """
    Units that were adjusted multiple times within a single settlement period
    are processed such that their entry is compressed into a single line.
    The 'LevelFrom' column is taken from the first entry, 'LevelTo' from the last 
    """
    
    new_rows = []

    for bmunit in (removals := (m := df["BmUnit"].value_counts()).loc[m > 1].index):
        
        new_row = (b := df.loc[df["BmUnit"] == bmunit]).iloc[-1].copy()
        new_row.loc["LevelFrom"] = b["LevelFrom"].iat[0]

        new_rows.append(new_row)
    
    df = df.loc[~df["BmUnit"].isin(removals)]

    return pd.concat([df, pd.concat(new_rows, axis=1).T], axis=0, ignore_index=True)


def infer_index(name, df, x):

    has_data = df[x] != 0.

    try:
        if "MOYL1" in name:
            return df.loc[(has_data) & (df.index.str.contains("MOYL1"))].index[0]
        
        else:
            return df.loc[(has_data) & (df.index.str.startswith(name[:3]))].index[0]
    except IndexError:
        return np.nan


def fill_interconnector_locations(df):
    """Fills in missing locations of interconnector BMUs based on existing data.
    Uses that many BMUs of are referring to the same interconnector."""

    x = df.columns.intersection(["x", "lon"])[0]
    y = df.columns.intersection(["y", "lat"])[0]

    def infer_location(name, df):
        if (df.loc[name, x] != 0.) or not name.startswith("I"):

            return df.loc[name, [x, y]].values
        
        elif ("MOYL1" in name) or (name[3] == "-"):

            other = infer_index(name, df, x)

        else:
            return df.loc[name, [x, y]].values
        
        if isinstance(other, str):
            return df.loc[other, [x, y]].values
        else:
            return [0., 0.]

    df.loc[:, x], df.loc[:, y] = zip(
        *df.reset_index().iloc[:,0].apply(lambda name: infer_location(name, df))
        )

    return df


def find_other_interconnectors(new, df):
    if isinstance(new, pd.Series):
        new = new.to_frame()
    
    others = new.reset_index().iloc[:,0].apply(lambda name: infer_index(name, df, "bus"))

    return others