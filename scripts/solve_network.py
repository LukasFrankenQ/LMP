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
    nodal_n = pypsa.Network(snakemake.input["network"])

    constraints = pd.read_csv(snakemake.input["network_constraints"], index_col=0)
    
    print(snakemake.wildcards.date)
    print(snakemake.wildcards.period)
    print(constraints)
    constraints = constraints["limit"]

    # print(constraints)
    # print(constraints.shape)

    layout = snakemake.wildcards.layout
    boundaries = snakemake.params["boundaries"][layout]

    # print(n.lines)
    # print(type(n.lines.index[0]))
    # print(n.lines.index[0])

    no_data = pd.Index(boundaries).difference(constraints.index)
    logger.info(f"No day-head constraint flow data for boundaries {no_data}. Filling in from nodal layout")
    
    def get_nodal_constraints(n, boundary):
        return n.lines.loc[pd.Index(boundaries[boundary], dtype=str), "s_nom"].sum()

    no_data = pd.Series(no_data, no_data).apply(lambda x: get_nodal_constraints(nodal_n, x))

    # print("before constraints", constraints)
    constraints = pd.concat([constraints, no_data])
    # print("after constraints", constraints)

    for boundary, lines in boundaries.items():
        # print(boundary, lines)

        if not len(lines):
            logger.info(f"Boundary {boundary} not representable in layout {layout} - no constraints enforced.")
            continue

        lines = pd.Series(lines).astype(str).tolist()

        value = constraints.loc[boundary]

        logger.info(f"Setting flow limits for boundary {boundary} to {int(value)} MW.")
        n.lines.loc[lines, "s_nom"] *= value / n.lines.loc[lines, "s_nom"].sum()

        # print(n.lines.loc[lines, "s_nom"].sum())

    # remove generator if p_nom is 0
    # n.generators.drop(n.generators.loc[n.generators.p_nom == 0].index, inplace=True)

    bus_ids = n.buses.index

    isolated_buses = check_network_consistency(n)
    logger.info(f"A total of {len(isolated_buses)} isolated buses:\n" + ",".join(isolated_buses))

    factor = snakemake.params["solving"]["p_nom_multiplier"]

    n.generators.loc[:, "p_nom"] *= (
        n.loads.p_set.sum() / 
        n.generators.p_nom.sum() 
        * factor
    )

    logger.warning("Solver configuration not yet taken from gurobi!")

    n.optimize(solver_name="gurobi")
    logger.warning("Solver does not yet check for infeasibility!")

    n.export_to_netcdf(snakemake.output["network"])
