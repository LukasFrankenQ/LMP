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

from _helpers import configure_logging

if __name__ == "__main__":

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input["network"])

    # remove generator if p_nom is 0
    # n.generators.drop(n.generators.loc[n.generators.p_nom == 0].index, inplace=True)

    print('Network to solver -')
    # Get all bus ids
    bus_ids = n.buses.index

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
    n.consistency_check()
    # n.export_to_csv_folder('csv')

    factor = snakemake.params["solving"]["p_nom_multiplier"]

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
