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
from _elexon_helpers import find_other_interconnectors

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
    
    total_export = abs(export.sum())

    logger.info(f"Total export capacity: {total_export} MW. Removed from load.")

    total_load -= total_export

    loads = load_weights.loc[load_weights > 0].index

    n.madd(
        "Load",
        loads,
        bus=loads,
        p_set=load_weights.loc[loads] * total_load,
    )

    bmu = bmu.loc[n.generators.index]

    logger.warning("Just deleting BMUS with negative Physical Notification for now.")
    bmu = bmu.loc[bmu["PN"] > 0].max(axis=1)

    n.generators.loc[bmu.index, 'p_nom'] = bmu

    missing = export.loc[export.index.difference(n.generators.index)]
    if len(missing):
        logger.info("Dropping missing exporting BMUs {}.".format(",".join(missing.index)))
        export = export.drop(missing.index)    

    n.generators.loc[(g := n.generators.loc[export.index]).loc[g["p_nom"] == 0.].index, "p_nom"] = export.abs()

    pu = export.div(n.generators.loc[export.index, 'p_nom'])

    n.generators.loc[export.index, 'p_max_pu'] = pu
    n.generators.loc[export.index, 'p_min_pu'] = pu

    logger.warning("no sensible costs yet!")
    bus_ids = n.buses.index

    isolated_buses = check_network_consistency(n)
    logger.info(f"A total of {len(isolated_buses)} isolated buses:\n" + ",".join(isolated_buses))

    n.generators = n.generators.loc[n.generators.p_nom > 0]
    n.export_to_netcdf(snakemake.output["network"])