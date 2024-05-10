# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT

import logging

import json
import pandas as pd
from tqdm import tqdm
from itertools import product

from _helpers import configure_logging
from _aggregation_helpers import get_demand

logger = logging.getLogger(__name__)


if __name__ == '__main__':

    configure_logging(snakemake)

    logger.info('Building demand totals.')

    index_mapper = {
        'domestic single': slice(0,48),
        'domestic multi': slice(52,100),
        'non-domestic single': slice(104,152),
        'non-domestic multi': slice(156,204),
    }

    single_rate_domestic = pd.read_excel(
        snakemake.input["elexon_demand_profiles"],
        index_col=0,
        header=2).iloc[index_mapper['domestic single']]
    multi_rate_domestic = pd.read_excel(
        snakemake.input["elexon_demand_profiles"],
        index_col=0,
        header=2).iloc[index_mapper['domestic multi']]

    single_rate_nondomestic = pd.read_excel(
        snakemake.input["elexon_demand_profiles"],
        index_col=0,
        header=2).iloc[index_mapper['non-domestic single']]
    multi_rate_nondomestic = pd.read_excel(
        snakemake.input["elexon_demand_profiles"],
        index_col=0,
        header=2).iloc[index_mapper['non-domestic multi']]
    
    demands = {
        'single_rate_domestic': [],
        'multi_rate_domestic': [],
        'single_rate_nondomestic': [],
        'multi_rate_nondomestic': [],
    }

    for date, period in tqdm(product(pd.date_range('2023', '2024', freq='d'), range(1, 49))):

        date = date.strftime('%Y-%m-%d')

        for key in demands.keys():
            demands[key].append(get_demand(locals()[key], date, period))


    for key, demand in demands.items():
        demands[key] = sum(demand)

    with open(snakemake.output[0], 'w') as f:
        json.dump(demands, f)

