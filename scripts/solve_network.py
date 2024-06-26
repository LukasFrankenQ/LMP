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
import numpy as np
import pandas as pd

from _helpers import configure_logging, check_network_consistency

if __name__ == "__main__":

    configure_logging(snakemake)

    n = pypsa.Network(snakemake.input["network"])
    nodal_n = pypsa.Network(snakemake.input["nodal_network"])

    constraints = pd.read_csv(snakemake.input["network_constraints"], index_col=0)

    constraints = constraints["limit"]

    layout = snakemake.wildcards.layout
    boundaries = snakemake.params["boundaries"][layout]
    boundaries_nodal = snakemake.params["boundaries"]['nodal']

    no_data = pd.Index(boundaries).difference(constraints.index)

    logger.info(f"No day-head constraint flow data for boundaries {','.join(no_data)}. Filling in from nodal layout.")
    
    def get_nodal_constraints(n, boundary):
        return n.lines.loc[pd.Index(boundaries_nodal[boundary], dtype=str), "s_nom"].sum()

    no_data = pd.Series(no_data, no_data).apply(lambda x: get_nodal_constraints(nodal_n, x))

    constraints = pd.concat([constraints, no_data])

    for boundary, lines in boundaries.items():

        if not len(lines):
            logger.info(f"Boundary {boundary} not representable in layout {layout} - no constraints enforced.")
            continue

        lines = pd.Series(lines).astype(str).tolist()
        current_capacity = n.lines.loc[lines, "s_nom"].sum()

        if boundary in ['SCOTEX', 'SSHARN']:
            # Western HVDC counts towards these boundaries, scaling accordingly
            current_capacity = current_capacity + n.links.loc['8009', 'p_nom']

        constraint_factor = np.around(
            constraints.loc[boundary] / current_capacity,
            decimals=3)

        verdict = "Accepted" if constraint_factor < 1 else "Skipped"
        logger.info(f"Reduction factor {constraint_factor} for boundary {boundary}: {verdict}.")

        if verdict == "Skipped":
            continue

        n.lines.loc[lines, "s_nom"] *= constraint_factor

        if boundary in ['SCOTEX', 'SSHARN']:
            # also reduce HVDC capacity
            n.links.loc['8009', 'p_nom'] *= constraint_factor
        
    isolated_buses = check_network_consistency(n)
    logger.info(f"A total of {len(isolated_buses)} isolated buses:\n" + ",".join(isolated_buses))

    factor = snakemake.params["solving"]["p_nom_multiplier"]

    logger.warning("Solver configuration not yet taken from gurobi!")

    # n.optimize(solver_name="gurobi")
    n.optimize(solver_name="highs")
    logger.warning("Solver does not yet check for infeasibility!")

    n.export_to_netcdf(snakemake.output["network"])