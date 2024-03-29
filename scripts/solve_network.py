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
    n.generators.drop(n.generators.loc[n.generators.p_nom == 0].index, inplace=True)

    print(n.buses)
    print(n.generators)
    print(n.loads)
    print(n.lines)
    print(n.links)
    n.consistency_check()
    n.export_to_csv_folder('csv')

    factor = snakemake.params["solving"]["p_nom_multiplier"]

    n.generators.loc[:, "p_nom"] *= (
        n.loads.p_set.sum() / 
        n.generators.p_nom.sum() 
        * factor
    )

    n.plot()
    plt.show()

    logger.warning("Solver configuration not yet taken from gurobi!")

    n.optimize(solver_name="gurobi")
    logger.warning("Solver does not yet check for infeasibility!")

    n.export_to_netcdf(snakemake.output["network"])
