# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT

import logging

logger = logging.getLogger(__name__)

import json
import numpy as np
import pandas as pd
from tqdm import tqdm
from itertools import product

from _helpers import configure_logging, to_quarter
from prepare_allowances import get_demand, get_weights

if __name__ == "__main__":

    configure_logging(snakemake)

    logger.info(f"Processing results to estimate monthly consumer prices.")

    year = snakemake.wildcards.year

    with open(snakemake.input["monthly_results"]) as f:
        results = json.load(f)

    policy_allowance_file = snakemake.input["policy_allowance"]
    dno_load_weights = pd.read_csv(snakemake.input["dno_load_weights"], index_col=0)["load_weight"]
                    

    print('loaded dno_load_weights')
    print(dno_load_weights)

    for ts, data in results.items():

        time = pd.Timestamp.fromtimestamp(int(ts))
        quarter = to_quarter(time)

        # single_standing = pd.read_csv(snakemake.input["single_st"], index_col=0)
        def get_allowance(files, quarter, order, mode):
            
            def check_fn(fn, quarter, order, mode):
                return (quarter in fn) and (order in fn) and (mode in fn)

            fn = [fn for fn in files if check_fn(fn, quarter, order, mode)][0]
            return pd.read_csv(fn, index_col=0) 

        single_standing = get_allowance(snakemake.input, quarter, 'zeroth', 'single')
        single_linear = get_allowance(snakemake.input, quarter, 'first', 'single')
        multi_standing = get_allowance(snakemake.input, quarter, 'zeroth', 'multi')
        multi_linear = get_allowance(snakemake.input, quarter, 'first', 'multi')

        zeroth_cols = {
            "policy": ["warm home discount"],
            "network": [col for col in multi_standing.columns if 'standing' in col],
            }
        linear_cols = {
            "policy": [
                "renewable obligation",
                "feed-in tariffs",
                "old feed-in tariffs",
                "energy company obligation",
                "aahedc",
                ],
            "network": [col for col in single_linear.columns if ' m ' in col],
        }

        single_demand = get_demand(policy_allowance_file, 'single')
        multi_demand = get_demand(policy_allowance_file, 'multi')

        single_weight = get_weights(policy_allowance_file, 'single')
        multi_weight = get_weights(policy_allowance_file, 'multi')
        dt_quarter = pd.Timestamp(quarter)

        # reflects that demand is considered yearly while we apply costs monthly
        if dt_quarter.month == 1 or dt_quarter.month == 10:
            single_demand = single_demand * single_weight["winter"] / 6.
            multi_demand = multi_demand * multi_weight["winter"] / 6.
        else:
            single_demand = single_demand * single_weight["summer"] / 6.
            multi_demand = multi_demand * multi_weight["summer"] / 6.

        for layout, layout_data in data.items():
            logger.info(f"Processing layout {layout} for month {time.strftime('%Y-%m')}.")

            if layout == "national":
                overlap = pd.DataFrame(1., index=single_standing.index, columns=["GB"])
            else:
                overlap = pd.read_csv(snakemake.input[f"overlap_{layout}"], index_col=0)

            layout_data = layout_data["geographies"] 

            for region, region_data in tqdm(layout_data.items()):

                for mode, cost_factor in product(["single", "multi"], ["policy", "network"]):

                    weighting = (w := overlap[region] * dno_load_weights) / w.sum()
                    weighted_mean = (
                        globals()[mode+'_standing'][zeroth_cols[cost_factor]]
                        .sum(axis=1)
                        .mul(weighting, axis=0)
                        .sum()
                    )

                    standing_total = weighted_mean / 12

                    try:
                        df = globals()[mode+'_linear'][linear_cols[cost_factor]]
                    except KeyError:
                        # very ad-hoc fix, should be made more robust
                        cols = [col.replace('3,100', '4,200') for col in linear_cols[cost_factor]]
                        df = globals()[mode+'_linear'][cols]

                    weighted_mean = (
                        df
                        .sum(axis=1)
                        .mul(weighting, axis=0)
                        .sum()
                    )

                    linear_total = weighted_mean * globals()[mode+'_demand']
                    layout_data[region]["variables"][f"{mode}-rate {cost_factor} cost"] = standing_total + linear_total

    with open(snakemake.output["monthly_results"], "w") as f:
        json.dump(results, f)
