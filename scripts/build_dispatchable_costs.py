# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This rules uses 6-months of dispatch data of BMUs to estimate the marginal cost of electricity generation
by correlating it with the wholesale price during the same time window.

**Outputs**

- ``RESOURCES/estimated_bmu_prices.csv``: estimated BMU prices

"""

import logging

logger = logging.getLogger(__name__)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
from sklearn.preprocessing import MinMaxScaler

from _helpers import configure_logging


def get_extremes(dispatch, mode='on', neighbour_filter=10):
    """
    Returns indices of timeseries of dispatch in which a generator is switched on or off.
    Works for high-cost peak-demand meeting generators 
    
    Parameters
    ----------
    dispatch : pd.Series
        Timeseries of dispatch
    mode : str, optional
        'on' or 'off', by default 'on'; on to return times of on-switching else off-switching
    """

    if isinstance(dispatch, pd.DataFrame):
        assert len(dispatch.columns) == 0, 'unclear how to handle multiple columns in dispatch dataframe'
        dispatch = dispatch.iloc[:,0]

    if dispatch.max() > 1.:
        dispatch = pd.Series(
            MinMaxScaler().fit_transform(dispatch.values.reshape(-1, 1)).flatten(),
            index=dispatch.index
        )
    
    rounded = dispatch.fillna(0.).round()

    deriv_window_size = 10
    deriv = (
        rounded
        .rolling(deriv_window_size, center=True)
        .apply(lambda x: x[deriv_window_size//2:].mean() - x[:deriv_window_size//2].mean())
    )

    if mode == 'on':
        extreme_func = np.less_equal
    elif mode == 'off':
        extreme_func = np.greater_equal
    else:
        raise ValueError('mode must be either "on" or "off"')

    extremes = np.array(argrelextrema(deriv.values, extreme_func, order=5)[0])

    extremes = extremes[deriv.iloc[extremes] != 0]

    # minor cleanup; remove neighbouring values
    mask = np.abs(np.roll(extremes, 1) - extremes) < neighbour_filter
    mask = mask + np.roll(mask, -1)

    extremes = extremes[~mask]

    return extremes


def get_bmu_price_distribution(dispatch, prices, visval=False, **kwargs):
    """
    For the dispatch of a single BMU, returns the average price of electricity
    during the time it is switched on. 

    Parameters
    ----------
    dispatch : pd.Series
        Timeseries of dispatch
    prices : pd.Series
        wholesale market price for electricity
    visval : bool, optional
        if True plot the dispatch and the inferred switch-on and switch-off times

    """

    full_index = pd.date_range(prices.index[0], prices.index[-1], freq='30min')

    prices = prices.copy().reindex(full_index).interpolate()
    dispatch = dispatch.copy().reindex(full_index).interpolate()

    dispatch = dispatch.loc[~dispatch.isna()]
    prices = prices.loc[dispatch.index]

    switchons = get_extremes(dispatch, mode='on', **kwargs)
    switchoffs = get_extremes(dispatch, mode='off', **kwargs)

    switchons = pd.Index([dispatch.index[n] for n in switchons])
    switchoffs = pd.Index([dispatch.index[n] for n in switchoffs])

    if visval:

        start = dispatch.index[0]
        if len(dispatch) > 1e3:
            end = dispatch.index[min(len(dispatch)-1, 1000)]

        fig, ax = plt.subplots(1, 1, figsize=(10, 3.5))
        ax.plot(dispatch, color='k')
        for on in switchons:
            ax.axvline(on, color='r', linestyle='--')
        for off in switchoffs:
            ax.axvline(off, color='b', linestyle='--')
        ax.set_title(dispatch.name)
        
        ax.set_xlim(start, end)

        ax.set_ylabel('Dispatch [MW]')
        ax.set_xlabel('Time')
        ax.xaxis.set_tick_params(rotation=45)

        plt.show()

    avg_prices = list()

    for on in switchons:

        try:
            off = switchoffs[switchoffs > on][0]
        except IndexError:
            continue

        avg_prices.append(prices.loc[on:off].mean())
    
    if len(avg_prices) == 0:
        return np.nan
    else:
        return avg_prices


if __name__ == "__main__":

    configure_logging(snakemake)

    logger.info("Estimating BMU prices from dispatch data of the last half year.")

    dispatch = pd.read_csv(snakemake.input["bmus_dispatch"], index_col=0, parse_dates=True)
    prices = pd.read_csv(snakemake.input["market_prices"], index_col=0, parse_dates=True).iloc[:,0]

    bmu_data = pd.read_csv(snakemake.input["bmus_data"], index_col=0)

    default_cost = 68 # £/MWh

    dispatchers = (
        (b := (
            pd.concat((
                bmu_data.loc[bmu_data.index.intersection(dispatch.columns)],
                dispatch.sum()
            ), axis=1)
            .rename(columns={0: 'generation'})
            .sort_values(by='generation', ascending=False)
        ))
        .loc[
            (b["carrier"].isin(['gas', 'CCGT', 'biomass', 'coal', 'CHP', 'gas turbine', 'powerstation']))
            & (b["generation"] > 0.)
        ]
    )

    total_generation = dispatchers.generation.sum()
    logger.info(f'total generation : {np.around(total_generation*1e-6, decimals=2)} TWh')

    avg_prices = {}

    from tqdm import tqdm
    for unit in tqdm(dispatchers.index):

        generator_prices = get_bmu_price_distribution(
                dispatch[unit],
                prices,
                visval=False,
                neighbour_filter=5)
        
        avg_prices[unit] = pd.Series(generator_prices).median()

    avg_prices = pd.Series(avg_prices)

    found_share = dispatch[avg_prices.index].sum().sum() / total_generation
    logger.info(f"Found share of generation: {found_share*100:.8f}%")

    avg_prices = avg_prices.fillna(default_cost)

    avg_prices.to_csv(snakemake.output["bmu_cost_estimates"])