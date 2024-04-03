# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
Simple method to solve

**Inputs**

- ``RESOURCES/{date}_{period}/network_s_{layout}.nc``: prepared network

**Outputs**

- ``RESOURCES/{date}_{period}/network_s_{layout}_solved.nc``: solved network

"""

import logging

logger = logging.getLogger(__name__)

import pypsa
import pandas as pd
import matplotlib.pyplot as plt

from _helpers import configure_logging, check_network_consistency

if __name__ == "__main__":

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input["network"])

    # remove generator if p_nom is 0
    # n.generators.drop(n.generators.loc[n.generators.p_nom == 0].index, inplace=True)

    print('Network to solver -')
    # Get all bus ids
    bus_ids = n.buses.index

    """
    # Check for each component if there are any attached to each bus
    mask = list()
    for bus_id in bus_ids:
        has_components = (
            not n.generators[n.generators.bus == bus_id].empty or
            not n.lines[n.lines.bus0 == bus_id].empty or
            not n.lines[n.lines.bus1 == bus_id].empty or
            not n.links[(n.links.bus0 == bus_id) | (n.links.bus1 == bus_id)].empty
        )

        has_load = not n.loads[n.loads.bus == bus_id].empty

        if has_load and not has_components:
            print(f"Bus {bus_id} has a load but no attached components.")
            mask.append(True)
        else:
            mask.append(False)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    ax.scatter(n.buses['x'], n.buses['y'], label='all')
    ax.scatter(n.buses.loc[mask, 'x'], n.buses.loc[mask, 'y'], label='mask')

    ax.legend()
    plt.show()

    print(n.buses)
    print(n.generators)
    print(n.loads)
    print(n.lines)
    print(n.links)
    """
    print('pypsa consistency check:')
    print(n.consistency_check())

    isolated_buses = check_network_consistency(n)
    logger.info(f"A total of {len(isolated_buses)} isolated buses:\n" + ",".join(isolated_buses))

    def interp_color(c1, c2, t):
        r1, g1, b1 = c1
        r2, g2, b2 = c2

        r = r1 * (1 - t) + r2 * t
        g = g1 * (1 - t) + g2 * t
        b = b1 * (1 - t) + b2 * t

        return (r/255., g/255., b/255.)

    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))#, subplot_kw={"projection": ccrs.PlateCarree()})

    def line_to_plot(bus0, bus1, quant):

        color = interp_color((255, 0, 0), (0, 255, 0), quant)
        ax.plot(
            [n.buses.loc[bus0, "x"], n.buses.loc[bus1, "x"]],
            [n.buses.loc[bus0, "y"], n.buses.loc[bus1, "y"]],
            color=color
            )

    ax.scatter(
        n.buses.loc[:, "x"],
        n.buses.loc[:, "y"],
        color="blue",
        label="all buses",
        alpha=0.5
        )
    ax.scatter(
        n.buses.loc[isolated_buses, "x"]
        , n.buses.loc[isolated_buses, "y"],
        color="red",
        label="isolated buses",
        alpha=0.5
        )

    point = 'r'
    for i, (_, line) in enumerate(n.lines.iterrows()):
        line_to_plot(
            line.bus0,
            line.bus1,
            # getattr(line, point) / getattr(n.lines, point).max()
            i / len(n.lines)
            )

    ax.legend()
    plt.show()
    print(n.lines.tail(30))

    factor = snakemake.params["solving"]["p_nom_multiplier"]

    bad = n.lines.loc[n.lines.x == 0.0].index
    good = n.lines.loc[n.lines.x != 0.0].index

    pd.set_option('display.max_columns', 500)
    print(n.lines.loc[n.lines.x == 0.0])
    print('avg values')
    print(n.lines.loc[n.lines.x != 0.0, 'x'].mean())

    # for quant in ["x", "r", "b", "x_pu", "r_pu", "b_pu", "r_pu_eff", "x_pu_eff"]:
    for quant in ["x", "r", "x_pu", "r_pu", "r_pu_eff", "x_pu_eff"]:
        n.lines.loc[n.lines[quant] == 0.0, quant] = n.lines.loc[n.lines[quant] != 0.0, quant].mean()

    n.lines.loc[n.lines["b"] == 0.0, "b"] = 0.0001
    n.lines.loc[n.lines["b_pu"] == 0.0, "b_pu"] = 16.
    
    # n.lines.loc[n.lines["i_nom"].isna()] = 2.58

    print('============= after =================')
    print(n.lines.columns)
    print('good\n: ', n.lines.loc[good].head())
    print('good\n: ', n.lines.loc[good].shape)
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('bad\n: ', n.lines.loc[bad].head())
    print('bad\n: ', n.lines.loc[bad].shape)


    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    for c in n.iterate_components(n.passive_branch_components):
        for attr in ["r", "x"]:
            print(attr)
            bad = c.df[attr] == 0.
            if bad.any():
                print(100 * bad.sum() / len(bad), "% of", attr, "are zero.")

    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')


    n.generators.loc[:, "p_nom"] *= (
        n.loads.p_set.sum() / 
        n.generators.p_nom.sum() 
        * factor
    )

    # n.plot()
    # plt.show()
    # Assuming 'n' is your Network object

    logger.warning("Solver configuration not yet taken from gurobi!")

    n.optimize(solver_name="gurobi")
    logger.warning("Solver does not yet check for infeasibility!")

    n.export_to_netcdf(snakemake.output["network"])
