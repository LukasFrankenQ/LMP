# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This rule collects live and other data to prepare the network for live operation.

**Outputs**

- ``RESOURCES/networks/prepared_live_{data}_{period}.csv``: prepared network

"""

import logging

import pypsa
import pandas as pd

from _helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":

    configure_logging(snakemake)

    date = snakemake.wildcards.date    
    period = snakemake.wildcards.period

    logger.info(f"Preparing live network for {date} settlement period {period}.")
    n = pypsa.Network(snakemake.input["network"])

    import pypsa
    import networkx as nx

    # Assuming 'network' is your PyPSA network
    graph = n.graph()

    # Find the connected components
    connected_components = list(nx.connected_components(graph))

    # The largest component is usually the main system
    main_system = max(connected_components, key=len)

    # print("Main system:", main_system)
    print("Number of buses in the main system:", len(main_system))
    # print("Buses:", network.buses.index)

    # Find the buses not in the main system
    isolated_buses = [bus for bus in n.buses.index if bus not in main_system]

    print("Isolated buses:", isolated_buses)
    print("Number of isolated buses:", len(isolated_buses))
    
    logger.warning("Should be Export Limit for dispatchable generators, but not yet implemented!")
    bmu = pd.read_csv(snakemake.input["elexon_bmus"]).set_index("NationalGridBmUnit")

    load_weights = pd.read_csv(snakemake.input["load_weights"], index_col=0)["load_weight"]

    total_load = bmu.loc[bmu["PN"] > 0, "PN"].sum()

    loads = load_weights.loc[load_weights > 0].index
    n.madd(
        "Load",
        loads,
        bus=loads,
        p_set=load_weights.loc[loads] * total_load,
    )

    bmu = bmu.loc[n.generators.index]

    logger.warning("Just deleting negative ones")
    bmu = bmu.loc[bmu["PN"] > 0].max(axis=1)

    n.generators.loc[bmu.index, 'p_nom'] = bmu

    logger.warning("no sensible costs yet!")

    print('Network after live network -')
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

    n.export_to_netcdf(snakemake.output["network"])