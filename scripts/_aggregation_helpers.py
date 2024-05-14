# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
Provides helper functions for aggregation of front-end friendly result data

"""

import pandas as pd
from itertools import product


method_mapper = {
    "wholesale_price": "mean",
    "post_balancing_price": "mean",
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
}

layouts = ['national', 'eso', 'fti', 'nodal']


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


def scale_stats(data, factor, inplace=False):
    """scales variables in data by factor. Only does so for variables that
    have aggregation method 'sum' in the method_mapper.
    Tailored for dicts as they are passed from .github/scripts/_live_helpers.py"""

    region_mapper = {
        l: list(data[l]["geographies"])
        for l in layouts
        }

    keys_template = "{layout},geographies,{region},variables,{variable}"

    for layout, (variable, method) in product(layouts, method_mapper.items()):

        if method != 'sum':
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