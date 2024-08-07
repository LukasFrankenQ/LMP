# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT

import sys
import json
import collections
import numpy as np
import pandas as pd
from pathlib import Path
from copy import deepcopy
from itertools import product

sys.path.append(str(Path.cwd() / 'scripts'))
from _aggregation_helpers import (
    set_nested_value,
    get_nested_value,
    layouts,
    flexible_aggregate,
    flexible_scale,
    aggregate_variable,
    remove_leaves,
    _get_key_chains,
)


def prepare_system_total(
        data,
        now_timestamp,
        cfd_consumer_benefit_share=1.
        ):

    agg = flexible_aggregate(data)

    if isinstance(now_timestamp, pd.Timestamp):
        now_timestamp = str(int(now_timestamp.timestamp()))

    total_savings = {
        list(data)[-1]: {'last_update': now_timestamp}
    }

    def _get_total_savings(d):
        return (
            d['balancing_cost_savings'] +
            d['congestion_rent_savings'] * cfd_consumer_benefit_share +
            d['wholesale_cost_savings'] +
            d['cfd_cost_savings']
        )

    for l in ['nodal', 'eso']:
            total_savings[list(total_savings)[-1]][f'{l}_total_savings'] = (
            _get_total_savings(agg[l]['globals']['variables'])
        )

    return total_savings


def prepare_household_total(
    regional_total,
    now_datetime,
    ):

    yearly_demands = {
        "single-rate-domestic": 4.211,
        "multi-rate-domestic": 5.495,
        "single-rate-nondomestic": 16.667,
        "multi-rate-nondomestic": 25.6824,
    }

    remove_leaves(
        regional_total,
        [
            # 'load',
            'wholesale_price',
            'post_policy_price',
        ]
    )

    start_datetime = pd.Timestamp.fromtimestamp(int(list(regional_total)[0]))
    
    # if savings are accumulated over more than a year, they are scaled down to one year
    scaling_factor = min(365 / (now_datetime - start_datetime).days, 1.)
    print('scaling factor to obtain one-yearly numbers: ', scaling_factor)

    agg = flexible_aggregate(regional_total)

    hh = {
        str(int((now_datetime - pd.Timedelta(days=365)).timestamp()))
        : agg
    }

    hh = flexible_scale(hh, scaling_factor)

    # to £/MWh
    for k in _get_key_chains(hh):

        for d in yearly_demands.keys():
            if d in k[-1]:
                set_nested_value(hh, k, get_nested_value(hh, k) / yearly_demands[d])

    hh['last_update'] = str(int(now_datetime.timestamp()))

    return hh


def to_last_month(dt):

    if isinstance(dt, str):
        dt = pd.Timestamp.fromtimestamp(int(dt))

    return pd.Timestamp(year=dt.year, month=dt.month, day=1, hour=1)


def get_variables(data):

    def deeper(key, item):
        if key == 'variables':
            return list(item)
        else:
            return deeper(list(item)[0], item[list(item)[0]])

    return deeper(list(data)[0], data[list(data)[0]])


def easy_aggregate(data):
    """Just sums over all timesteps"""

    data = collections.OrderedDict(sorted(data.items()))

    target_data = deepcopy(data[list(data)[0]])

    keys_template = "{layout},geographies,{region},variables,{variable}"

    region_mapper = {
        l: list(list(data.values())[0][l]["geographies"])
        for l in layouts
        }
    
    variables = get_variables(data)

    for layout, variable in product(layouts, variables):
        for region in region_mapper[layout]:

            keychain = (
                keys_template
                .format(layout=layout, region=region, variable=variable)
                .split(',')
            )

            set_nested_value(
                target_data,
                keychain,
                aggregate_variable(
                    data,
                    keychain,
                    "sum",
                    None,
                    )
                )

    return {list(data)[0]: target_data}


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

        monthly[list(monthly)[-1]] = flexible_aggregate({
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

        monthly[str(int(mstart.timestamp()))] = now[now_ts]

    '''
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
    '''

    monthly = collections.OrderedDict(sorted(monthly.items()))
    return monthly


def update_daily(daily, new, date):

    daily = deepcopy(daily)
    new = deepcopy(new)

    day_td = str(int(pd.Timestamp(date).timestamp()))

    if not day_td in daily:
        daily[day_td] = new[list(new)[0]]
        return daily

    to_agg = {day_td: daily[day_td]}
    to_agg.update(new)

    daily.update({day_td: flexible_aggregate(to_agg)})

    return daily


hh_keepers = ['post_balancing_price', 'wholesale_price', 'load']    
def half_hourly_func(d, key):
    """
    Modify the given dictionary 'd' by converting the value associated with 'key' to np.float16.
    If 'key' is not in the list of hh_keepers, remove it from the dictionary.
    """
    if key in hh_keepers:
        pass
    else:
        del d[key]


with open(
    Path(__file__).parent.parent.parent /
    'data' /
    'demand_totals.json', 'r'
    ) as f:
    demand_totals = json.load(f)
key_mapper = {}


def daily_func(d, old_key):
    """
    Modify the given dictionary 'd' by replacing the 'old_key' with a new key based
    on the key_mapper dictionary.
    If the new key exists in demand_totals, calculate the new value based on the old
    value and demand_totals[new_key]. Finally, remove the 'old_key' from the dictionary.
    (unless it is 'load')
    """
    if not old_key in key_mapper:

        if 'post_balancing_price' == old_key:
            key_mapper[old_key] = 'post_balancing_price'
        else:
            key_mapper[old_key] = old_key.split('_')[0].replace('-', '_')
    
    new_key = key_mapper[old_key]

    if new_key in demand_totals:
        d[new_key] = d[old_key] / demand_totals[new_key] * 1e3

    if old_key != new_key:
        del d[old_key]


def summary_func(d, old_key):
    """
    Modify the given dictionary 'd' by replacing the 'old_key' with a new key based
    on the key_mapper dictionary.
    If the new key exists in demand_totals, calculate the new value based on the old
    value and demand_totals[new_key]. Finally, remove the 'old_key' from the dictionary.
    (unless it is 'load')
    """
    if not old_key in key_mapper:
        key_mapper[old_key] = old_key.split('_')[0].replace('-', '_')
    
    new_key = key_mapper[old_key]

    if new_key in demand_totals:
        d[new_key] = d[old_key] / demand_totals[new_key] * 1e3

    if old_key != new_key and new_key in demand_totals:
        del d[old_key]


def prepare_frontend_dict(data, prep_func):
    """
    Apply the given 'prep_func' to each key-value pair in the 'data' dictionary.
    The 'prep_func' should modify the dictionary in place.
    """

    hold = deepcopy(data)

    def traverse_n_apply(d, target):
        for key in d.keys():
            if isinstance(d[key], dict):
                traverse_n_apply(d[key], target[key])
            else:
                prep_func(target, key)

    traverse_n_apply(data, hold)

    return hold


region_assignment = {
    "GB0 Z1_1": "GB0 Z1_3",
    "GB0 Z4": "GB0 Z3",
    "GB0 Z1_2": "GB0 Z1_3",
}

def fix_zonal_remote_regions(data):

    def traverse_n_apply(d, target):
        for key in d.keys():
            if key == 'globals':
                continue

            if not key == 'geographies':
                traverse_n_apply(d[key], target[key])
            else:
                for remote_region, source_region in region_assignment.items():
                    target[key][remote_region] = target[key][source_region]
                    target[key][remote_region]['variables']['load'] = 0.

    for key in data.keys():

        hold = deepcopy(data[key]['eso'])

        traverse_n_apply(data[key]['eso'], hold)
        data[key]['eso'] = hold


def prepare_constituency_total(hh_total, const_mapper_file):
    '''Aggregates household data to constituency level'''

    total_nodal = hh_total[list(hh_total)[0]]['nodal']['geographies']
    total_nodal = pd.DataFrame(
        {key: item['variables'] for key, item in total_nodal.items()}
        ).T.reset_index()
    total_nodal['index'] = total_nodal['index'].astype(str)

    total_zonal = hh_total[list(hh_total)[0]]['eso']['geographies']
    total_zonal = pd.DataFrame(
        {key: item['variables'] for key, item in total_zonal.items()}
        ).T
    
    const = pd.read_csv(const_mapper_file)
    const['nodal_region'] = const['nodal_region'].astype(str)

    const = (
        const.merge(
            total_zonal[['single-rate-domestic_total_savings']].reset_index(),
            left_on='zonal_region',
            right_on='index',
            how='left'
        )
        .drop(columns='index') 
        .rename(columns={
            'single-rate-domestic_total_savings': 'single_rate_domestic_zonal'
            })
    )

    const = (
            const.merge(
            total_nodal[['index', 'single-rate-domestic_total_savings']],
            left_on='nodal_region',
            right_on='index',
            how='left'
        )
        .drop(columns='index')
        .rename(columns={
            'single-rate-domestic_total_savings': 'single_rate_domestic_nodal'
            })
    )

    const_results = {}
    for l in ['nodal', 'zonal']:

        # assumes household savings are already in £/MWh
        const[f'savings_{l}'] = const[f'single_rate_domestic_{l}'].mul(const['demand'])
        const_results[f'const_savings_{l}'] = const.groupby('constituency_code')[f'savings_{l}'].sum()

    f = pd.concat((
        const_results['const_savings_nodal'],
        const_results['const_savings_zonal'],
        ), axis=1)
        
    fdict = dict(f.T)

    for key, item in fdict.items():
        fdict[key] = {'variables': item.to_dict()}

    fdict = {list(hh_total)[0]: {'constituencies': {'geographies': fdict}}}

    return fdict
