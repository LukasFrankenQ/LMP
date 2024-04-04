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

from _helpers import configure_logging, check_network_consistency

logger = logging.getLogger(__name__)

if __name__ == "__main__":

    configure_logging(snakemake)

    date = snakemake.wildcards.date    
    period = snakemake.wildcards.period

    logger.info(f"Preparing live network for {date} settlement period {period}.")
    n = pypsa.Network(snakemake.input["network"])

    isolated_buses = check_network_consistency(n)
    logger.info(f"A total of {len(isolated_buses)} isolated buses:\n" + ",".join(isolated_buses))

    logger.warning("Should be Export Limit for dispatchable generators, but not yet implemented!")
    bmu = pd.read_csv(snakemake.input["elexon_bmus"]).set_index("NationalGridBmUnit")

    load_weights = pd.read_csv(snakemake.input["load_weights"], index_col=0)["load_weight"]

    total_load = bmu.loc[bmu["PN"] > 0, "PN"].sum()

    logger.info("Moving exporting interconnectors from load to generators.")
    interconnectors = n.links[n.links.carrier == "AC"]

    export = bmu.loc[(
        bmu.index.str.startswith('I') &
        (bmu.index.str[3] == '-') &
        bmu.index.str[:3].str.isalpha() &
        (bmu["PN"] < 0)
        ), "PN"]
    print('exporters')
    print(export)
    
    total_export = abs(export.sum())

    logger.info(f"Total export capacity: {total_export} MW. Removed from load.")

    print('before total load')
    print(total_load)
    total_load -= total_export
    print('after total load')
    print(total_load)

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

    print('current state of exporters')
    print(n.generators.loc[export.index, ['p_nom', 'marginal_cost', 'p_max_pu', 'p_min_pu']])

    pu = export.div(n.generators.loc[export.index, 'p_nom'])
    print('export')
    print(export)

    print('pu')
    print(pu)

    n.generators.loc[export.index, 'p_max_pu'] = pu
    n.generators.loc[export.index, 'p_min_pu'] = pu

    print('after state of exporters')
    print(n.generators.loc[export.index, ['p_nom', 'marginal_cost', 'p_max_pu', 'p_min_pu']])

    logger.warning("no sensible costs yet!")

    print('Network after live network -')
    # Get all bus ids
    bus_ids = n.buses.index

    # Check for each component if there are any attached to each bus
    '''
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
    '''

    isolated_buses = check_network_consistency(n)
    logger.info(f"A total of {len(isolated_buses)} isolated buses:\n" + ",".join(isolated_buses))

    n.export_to_netcdf(snakemake.output["network"])