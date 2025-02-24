# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
Provides helper functions for aggregation of front-end friendly result data

"""

import pandas as pd
from copy import deepcopy
from functools import reduce
from itertools import product


method_mapper = {
    "wholesale_price": "mean",
    "post_policy_price": "mean",
    "load": "sum",
    # "dispatch": "mean",
    # "available_capacity": "mean",
    "single-rate-domestic_wholesale_savings": "sum",
    "multi-rate-domestic_wholesale_savings": "sum",
    "single-rate-nondomestic_wholesale_savings": "sum",
    "multi-rate-nondomestic_wholesale_savings": "sum",
    "single-rate-domestic_total_savings": "sum",
    "multi-rate-domestic_total_savings": "sum",
    "single-rate-nondomestic_total_savings": "sum",
    "multi-rate-nondomestic_total_savings": "sum",
    "single_rate_domestic": "sum",
    "multi_rate_domestic": "sum",
    "single_rate_nondomestic": "sum",
    "multi_rate_nondomestic": "sum",
    "wholesale_cost": "sum",
    "balancing_cost": "sum",
    "congestion_rent": "sum",
    "cfd_cost": "sum",
    "wholesale_cost_savings": "sum",
    "balancing_cost_savings": "sum",
    "congestion_rent_savings": "sum",
    "cfd_cost_savings": "sum",
}

# layouts = ['national', 'eso', 'fti', 'nodal']
layouts = ['national', 'eso', 'nodal']


def get_demand(demand, date, period):
    """Retrieve the demand for a specific date and period from the demand data."""

    day_mapper = {
        i: 'Wd' for i in range(5)
    }
    day_mapper.update({
        5: 'Sat', 6: 'Sun'
    })

    '''Seasons according to OFGEM (passive aggressive tone accidental)'''
    season_mapper = pd.concat((
        pd.Series('Smr', pd.date_range('2022-05-03', '2022-07-16', freq='D')),
        pd.Series('Hsr', pd.date_range('2022-07-17', '2022-08-24', freq='D')),
        pd.Series('Aut', pd.date_range('2022-08-25', '2022-10-27', freq='D')),
        pd.Series('Wtr', pd.date_range('2022-10-28', '2023-03-31', freq='D')),
        pd.Series('Spr', pd.date_range('2023-04-01', '2023-05-02', freq='D')),
        pd.Series('Smr', pd.date_range('2023-05-03', '2023-07-16', freq='D')),
        pd.Series('Hsr', pd.date_range('2023-07-17', '2023-08-24', freq='D')),
        pd.Series('Aut', pd.date_range('2023-08-25', '2023-10-27', freq='D')),
        pd.Series('Wtr', pd.date_range('2023-10-28', '2024-03-31', freq='D')),
        pd.Series('Spr', pd.date_range('2024-04-01', '2024-05-02', freq='D')),
        pd.Series('Smr', pd.date_range('2024-05-03', '2024-07-16', freq='D')),
        pd.Series('Hsr', pd.date_range('2024-07-17', '2024-08-24', freq='D')),
        pd.Series('Aut', pd.date_range('2024-08-25', '2024-10-27', freq='D')),
        pd.Series('Wtr', pd.date_range('2024-10-28', '2025-03-31', freq='D')),
        pd.Series('Spr', pd.date_range('2025-04-01', '2025-05-02', freq='D')),
        pd.Series('Smr', pd.date_range('2025-05-03', '2025-07-16', freq='D')),
        pd.Series('Hsr', pd.date_range('2025-07-17', '2025-08-24', freq='D')),
    ))

    if isinstance(date, str):
        date = pd.Timestamp(date)
    if isinstance(period, str):
        period = int(period)

    col = f"{season_mapper.loc[date]} {day_mapper[date.weekday()]}"
    return demand[col].iloc[period-1]


def get_timestep_weights(data):
    """Averaging over time should be weighted by load and not by time.
    This function calculates the weights for each timestep based on the load."""    

    loads = pd.Series([
        item['national']["geographies"]["GB"]["variables"]["load"] for item in data.values()
    ], index=data.keys()
    )
    loads /= loads.sum()

    return loads


def get_nested_value(origin_dict, keys):
    """Retrieve a value from a nested dictionary using a list of keys."""
    current_level = origin_dict
    for key in keys:
        current_level = current_level[key]
    return current_level


def set_nested_value(target_dict, keys, value):
    """Set a value in a nested dictionary using a list of keys. If a key
    does not exist, it is created."""
    current_level = target_dict
    for key in keys[:-1]:

        if key not in current_level:
            current_level[key] = {}

        current_level = current_level[key]
    current_level[keys[-1]] = value


def aggregate_variable(data, keys, method, weights=None):
    """Group the data by a certain variable and aggregate using a certain method."""

    if weights is None:
        weights = pd.Series(1, index=data.keys())
    elif method == 'mean' and not weights is None:
        weights = weights.copy() / weights.sum() * len(weights)

    return getattr(
        pd.Series([
            w * get_nested_value(item, keys)
              for w, item in zip(weights, data.values())
        ]), method)()


def aggregate_stats(origin_data):
    """Aggregates the statistics from the origin data.
    Returns a (key, item) pair where the key is the first timestep in origin data.
    """

    region_mapper = {
        l: list(list(origin_data.values())[0][l]["geographies"])
        for l in layouts
        }

    weights = get_timestep_weights(origin_data)
    weights_mapper = {
        "wholesale_price": weights,
        "post_balancing_price": weights,
        # "load": None,
        # "dispatch": None,
        # "available_capacity": None,
    }

    target_data = {}

    keys_template = "{layout},geographies,{region},variables,{variable}"

    for layout, (variable, method) in product(layouts, method_mapper.items()):
        for region in region_mapper[layout]:

            keychain = (
                keys_template
                .format(layout=layout, region=region, variable=variable)
                .split(',')
            )

            try:
                set_nested_value(
                    target_data,
                    keychain,
                    aggregate_variable(
                        origin_data,
                        keychain,
                        method,
                        weights_mapper.get(variable, None)
                        )
                    )
            except KeyError:
                pass
    
    return target_data


def get_variables(data):
    return list(data['national']['geographies']['GB']['variables'])


def scale_stats(data, factor, inplace=False):
    """scales variables in data by factor. Only does so for variables that
    have aggregation method 'sum' in the method_mapper.
    Tailored for dicts as they are passed from .github/scripts/_live_helpers.py"""

    region_mapper = {
        l: list(data[l]["geographies"])
        for l in layouts
        }
    data_vars = get_variables(data)

    keys_template = "{layout},geographies,{region},variables,{variable}"

    for layout, (variable, method) in product(layouts, method_mapper.items()):

        if method != 'sum':
            continue
    
        if variable not in data_vars:
            continue

        for region in region_mapper[layout]:

            keychain = (
                keys_template
                .format(layout=layout, region=region, variable=variable)
                .split(',')
            )

            set_nested_value(
                data,
                keychain,
                get_nested_value(data, keychain) * factor
                )
    
    if not inplace:
        return data


def _get_key_chains(nested_dict, current_chain=None):
    if current_chain is None:
        current_chain = []

    if isinstance(nested_dict, dict):
        key_chains = []
        for key, value in nested_dict.items():
            key_chains.extend(_get_key_chains(value, current_chain + [key]))
        return key_chains
    else:
        return [current_chain]


def _add_to_nested_dict(d, keys, value):
    reduce(lambda d, key: d.setdefault(key, {}), keys[:-1], d)[keys[-1]] = value


def flexible_aggregate(data):
    
    keychains = _get_key_chains(data[list(data)[0]])

    agg = {}

    for keychain in keychains:
        method = method_mapper.get(keychain[-1], None)
    
        assert method is not None, f"Method for {keychain[-1]} not defined in method_mapper."

        _add_to_nested_dict(
            agg,
            keychain,
            aggregate_variable(data, keychain, method)
        )

    return agg


def flexible_scale(data, factor):
    """Scales all leaves in a nested dictionary by a factor."""

    key_chains = _get_key_chains(data)
    hold = deepcopy(data)

    for k in key_chains:
        set_nested_value(hold, k, get_nested_value(hold, k) * factor)

    return hold


def _remove_key_from_nested_dict(data, key_chain):
    if not key_chain:
        return

    current = data
    for key in key_chain[:-1]:
        if key in current:
            current = current[key]
        else:
            print(f"Key {key} not found in the dictionary")
            return
    
    final_key = key_chain[-1]
    if final_key in current:
        del current[final_key]


def remove_leaves(d, leavers):
    '''Removes from d all leaves that end in any of the strings in leavers'''

    _key_chain = _get_key_chains(d)

    for k in _key_chain:
        if k[-1] in leavers:

            _remove_key_from_nested_dict(d, k)
