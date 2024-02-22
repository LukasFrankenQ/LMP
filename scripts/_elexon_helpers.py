# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT


import pandas as pd

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