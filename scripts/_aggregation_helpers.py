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


def group_quantity(data, keys, method, weights=None):
    """Group the data by a certain variable and aggregate using a certain method."""

    if method == 'mean' and weights is None:
        weights = pd.Series(1/len(data), index=data.keys())
    elif method == 'sum' and weights is None:
        weights = pd.Series(1, index=data.keys())
    elif method == 'mean' and not weights is None:
        weights = weights.copy() / weights.sum() * len(weights)

    return getattr(
        pd.Series([
            w * get_nested_value(item, keys)
              for w, item in zip(weights, data.values())
        ]), method)()


def aggregate_stats(
        origin_data,
        layouts
        ):
    """Aggregates the statistics from the origin data.
    returns a key, item pair where the key is the first timestep in origin data.
    """

    layouts = ['national', 'eso', 'fti', 'nodal']
    region_mapper = {
        l: list(list(origin_data.values())[0][l]["geographies"])
        for l in layouts
        }

    method_mapper = {
        "wholesale_price": "mean",
        "post_balancing_price": "mean",
        "load": "sum",
        "dispatch": "mean",
        "available_capacity": "mean",
    }

    weights = get_timestep_weights(origin_data)
    weights_mapper = {
        "wholesale_price": weights,
        "post_balancing_price": weights,
        "load": None,
        "dispatch": None,
        "available_capacity": None,
    }

    target_data = {}

    keys_template = "{layout},geographies,{region},variables,{variable}"

    for layout, (variable, method) in product(layouts, method_mapper.items()):
        for region in region_mapper[layout]:

            keychain = keys_template.format(layout=layout, region=region, variable=variable).split(',')

            set_nested_value(
                target_data,
                keychain,
                group_quantity(
                    origin_data,
                    keychain,
                    method,
                    weights_mapper[variable]
                    )
                )
    
    return target_data