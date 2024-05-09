# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT

import sys
import json
import collections
import pandas as pd
from pathlib import Path
from copy import deepcopy
from itertools import product

sys.path.append(str(Path.cwd() / 'scripts'))
from _aggregation_helpers import (
    get_nested_value,
    set_nested_value,
    layouts,
    method_mapper,
    scale_stats,
    aggregate_stats,
)

def to_last_month(dt):

    if isinstance(dt, str):
        dt = pd.Timestamp.fromtimestamp(int(dt))
    
    return pd.Timestamp(year=dt.year, month=dt.month, day=1, hour=1)


def update_monthly(now: dict, monthly: dict) -> dict:
    """
    Updates 'monthly' data to by absorbing 'now' data into the
    latest month, and cutting back the values in the oldest month.  
    """

    monthly = deepcopy(monthly)
    now = deepcopy(now)

    # ensure monthly is sorted
    monthly = collections.OrderedDict(sorted(monthly.items()))

    # indicates that new month data is already present; updates it
    if to_last_month(list(monthly)[-1]) == to_last_month(list(now)[0]):

        monthly[list(monthly)[-1]] = aggregate_stats({
            list(monthly)[-1]: monthly[list(monthly)[-1]],
            list(now)[0]: now[list(now)[0]]
            })
    # indicates that month data is not present; adds it
    else:
        now_ts = list(now)[0]
        now_dt = pd.Timestamp.fromtimestamp(int(now_ts))

        if now_dt.day == 1:
            mstart = pd.Timestamp(
                year=now_dt.year,
                month=now_dt.month,
                day=1,
                hour=1
                )
        else:
            mstart = (now_dt - pd.offsets.MonthBegin(1)).floor('d')

        now[str(int(mstart.timestamp()))] = now.pop(now_ts)

    now_ts = list(now)[0]
    now_dt = pd.Timestamp.fromtimestamp(int(now_ts))

    # potentially removes oldest data, updates first timestamp
    monthly_max_length = pd.Timedelta('365 days')

    # remove months that are too far into the past
    while (
        (
            (now_dt - monthly_max_length).month !=
            (mdt := pd.Timestamp.fromtimestamp(int(list(monthly)[0]))).month
            )
        & 
        (
            now_dt - mdt > monthly_max_length
            )
        ):
        del monthly[list(monthly)[0]]

    monthly_ts = list(monthly)[0]
    monthly_dt = pd.Timestamp.fromtimestamp(int(monthly_ts))

    delta = now_dt - monthly_dt

    if delta > monthly_max_length:

        mend = monthly_dt + pd.offsets.MonthBegin(1)
        remaining = mend - monthly_dt

        too_much = delta - monthly_max_length
        reduction_factor = 1 - too_much / remaining

        if reduction_factor == 0:
            del monthly[monthly_ts]

        else:
            new_ts = str(int((now_dt - monthly_max_length).timestamp()))

            monthly[new_ts] = scale_stats(
                monthly.pop(monthly_ts),
                reduction_factor
                )

    monthly = collections.OrderedDict(sorted(monthly.items()))
    return monthly