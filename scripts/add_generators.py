# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
Assigns located generators to buses based on regions

**Outputs**

- ``RESOURCES/networks/gen.nc``: Network with generators assigned to buses

"""

import logging

import yaml
import pypsa
import pandas as pd
import geopandas as gpd

from _helpers import configure_logging
from cluster_network import make_busmap

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    configure_logging(snakemake)

    logger.info("Assigning generators to buses based on regions.")

    n = pypsa.Network(snakemake.input["base_network"])

    bmus = pd.read_csv(snakemake.input["bmunits_loc"])
    bmus = bmus.loc[bmus.lat != 0.]

    bmus = (
        gpd.GeoDataFrame(bmus, geometry=gpd.points_from_xy(bmus.lon, bmus.lat))
        [["geometry", "NationalGridBmUnit", "capacity", "carrier"]]
        ).set_crs(epsg=4326)

    onshore = gpd.read_file(snakemake.input["regions_onshore"]).set_index("name")
    offshore = gpd.read_file(snakemake.input["regions_offshore"]).set_index("name")

    custom_busmap = make_busmap(n, onshore)

    '''
    non = custom_busmap.loc[custom_busmap.isna()].index

    logger.warning(f"Excluding {len(non)} buses from clustering with an ad-hoc method.")
    for c in n.iterate_components(n.one_port_components):
        (c := c.df).drop(c.loc[c.bus.isin(non)].index, inplace=True)

    for c in n.iterate_components(n.branch_components):
        (c := c.df).drop(c.loc[(c.bus0.isin(non)) | (c.bus1.isin(non))].index, inplace=True)

    n.buses.drop(non, inplace=True)
    '''

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

    non_assigned = 0

    for _, bmu in bmus.iterrows():
        
        hit = onshore.loc[onshore.geometry.contains(bmu.geometry)].index
        if len(hit) == 1:
            n.add("Generator", bmu["NationalGridBmUnit"], bus=hit[0])
            continue

        hit = offshore.loc[offshore.geometry.contains(bmu.geometry)].index
        if len(hit) == 1:
            n.add("Generator", bmu["NationalGridBmUnit"], bus=hit[0])
            continue

        non_assigned += 1

    logger.info(f"Assigned {len(n.generators)} generators to buses based on regions.")
    logger.info(f"Failed to assign {non_assigned} generators to buses based on regions.")
    logger.info(f"Saving network to {snakemake.output['gen_network']}.")

    n.generators.loc[:, 'carrier'] = bmus.set_index("NationalGridBmUnit").loc[n.generators.index, 'carrier']
    logger.info(f"Added carriers to generators. Share Unknown: {n.generators.carrier.isna().sum()/len(n.generators)}.")

    logger.info(f"Adding marginal costs to generators from '{snakemake.input['carrier_costs']}'.")
    costs = yaml.safe_load(open(snakemake.input["carrier_costs"]))
    n.generators.loc[:, "marginal_cost"] = n.generators.carrier.apply(lambda carrier: costs.get(carrier, 100))

    n.export_to_netcdf(snakemake.output["gen_network"])