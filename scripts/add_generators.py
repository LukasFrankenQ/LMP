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
from scipy import stats
from sklearn.preprocessing import MinMaxScaler

from _helpers import configure_logging, check_network_consistency
from _elexon_helpers import fill_interconnector_locations
from cluster_network import make_busmap

logger = logging.getLogger(__name__)


if __name__ == "__main__":
    configure_logging(snakemake)

    logger.info("Assigning generators to buses based on regions.")

    n = pypsa.Network(snakemake.input["base_network"])
    nds = snakemake.params["electricity"]["network_dataset"]

    onshore = gpd.read_file(snakemake.input["regions_onshore"]).set_index("name")
    offshore = gpd.read_file(snakemake.input["regions_offshore"]).set_index("name")

    if nds == "ENTSO-E":

        # gb_shape = gpd.read_file(snakemake.input["gb_shape"]).geometry[0]
        # 
        # buses_geom = gpd.GeoDataFrame(
        #    index=n.buses.index,
        #     geometry=gpd.points_from_xy(
        #         n.buses["x"].values, n.buses["y"].values
        #         )
        #     ).set_crs(epsg=4326)
        # outside_gb = buses_geom.loc[~buses_geom.within(gb_shape)].index

        # import matplotlib.pyplot as plt
        # fig, ax = plt.subplots()
        # buses_geom.plot(ax=ax, color="black", label="All")
        # buses_geom.loc[outside_gb].plot(ax=ax, color="red", label="Outside GB")
        # ax.legend()
        # plt.show()

        custom_busmap = make_busmap(n, onshore)

        outsiders = custom_busmap.loc[custom_busmap.isna()].index

        logger.warning(f"Excluding {len(outsiders)} buses from network with an ad-hoc method.")
        for c in n.iterate_components(n.one_port_components):
            (c := c.df).drop(c.loc[c.bus.isin(outsiders)].index, inplace=True)

        for c in n.iterate_components(n.branch_components):
            (c := c.df).drop(c.loc[(c.bus0.isin(outsiders)) | (c.bus1.isin(outsiders))].index, inplace=True)

        n.buses.drop(outsiders, inplace=True)

    else:
        raise ValueError(f"Currently only ENTSO-E network dataset is supported, found {nds}.")

    bmus = pd.read_csv(snakemake.input["bmunits_loc"])

    bmus = fill_interconnector_locations(bmus.set_index("NationalGridBmUnit")).reset_index()

    bmus = bmus.loc[bmus.lat != 0.]

    bmus = (
        gpd.GeoDataFrame(bmus, geometry=gpd.points_from_xy(bmus.lon, bmus.lat))
        [["geometry", "NationalGridBmUnit", "capacity", "carrier"]]
        ).set_crs(epsg=4326)

    isolated_buses = check_network_consistency(n)
    logger.info(f"A total of {len(isolated_buses)} isolated buses:\n" + ",".join(isolated_buses))

    non_assigned = 0

    for _, bmu in bmus.iterrows():
        
        hit = onshore.loc[onshore.geometry.contains(bmu.geometry)].index
        if len(hit) >= 1:
            n.add("Generator", bmu["NationalGridBmUnit"], bus=hit[0])
            continue

        hit = offshore.loc[offshore.geometry.contains(bmu.geometry)].index
        if len(hit) >= 1:
            n.add("Generator", bmu["NationalGridBmUnit"], bus=hit[0])
            continue

        non_assigned += 1

    logger.info(f"Assigned {len(n.generators)} generators to buses based on regions.")
    logger.info(f"Failed to assign {non_assigned} generators to buses based on regions.")
    logger.info(f"Saving network to {snakemake.output['gen_network']}.")

    n.generators.loc[:, 'carrier'] = bmus.set_index("NationalGridBmUnit").loc[n.generators.index, 'carrier']
    n.generators.loc[
        (n.generators.index.str.startswith("I")) &
        (n.generators.index.str[3] == "-"),
        "carrier"] = "interconnector"

    logger.info(f"Added carriers to generators. Share Unknown: {n.generators.carrier.isna().sum()/len(n.generators)}.")

    logger.info(f"Adding marginal costs to generators from '{snakemake.input['carrier_costs']}'.")
    costs = yaml.safe_load(open(snakemake.input["carrier_costs"]))
    
    n.generators.loc[:, "marginal_cost"] = n.generators.carrier.apply(lambda carrier: costs.get(carrier, 100))

    logger.info("Inserting cost estimates for dispatchable generators.")

    bmu_costs = pd.read_csv(snakemake.input["bmu_cost_estimates"], index_col=0)
    inter = n.generators.index.intersection(bmu_costs.index)
    inter = bmu_costs.loc[inter].loc[bmu_costs.loc[inter, bmu_costs.columns[0]] != 68.].index

    inserted_costs = bmu_costs.loc[inter].values

    def change_cost_spread(costs, factor):
        costs = costs.copy()
        min_cost = costs.min()

        return (costs - min_cost) * factor + min_cost    


    logger.info(f"Applying spread factor: {(f:= snakemake.params['elexon']['cost_spread'])}.")
    inserted_costs = change_cost_spread(
        inserted_costs, f
        )

    n.generators.loc[inter, ["marginal_cost"]] = inserted_costs

    def sample_more(data, n):
        '''Fits distribution to data, samples n values from it.'''

        scaler = MinMaxScaler()
        data = (
            scaler.fit_transform(
                data.copy()
                .values
                .reshape(-1, 1)
                )
            .flatten()
        )

        a, b, loc, scale = stats.beta.fit(data)
        new = stats.beta.rvs(a, b, loc=loc, scale=scale, size=n)

        return scaler.inverse_transform(new.reshape(-1, 1)).flatten()
    

    logger.warning("Assuming default marginal cost for generators without cost estimates as 68 £/MWh.")
    replacers = n.generators.loc[n.generators.marginal_cost == 68.].index

    n.generators.loc[replacers, "marginal_cost"] = sample_more(
        n.generators.loc[inter, ["marginal_cost"]],
        len(replacers)
        )

    # these costs will later be tuned according to wholesale prices to match the real market
    (
        pd.Series(
            replacers.tolist() + inter.tolist(),
            name='cost_esimated_generators'
            )
            .to_csv(snakemake.output['cost_estimated_generators'])
    )

    n.export_to_netcdf(snakemake.output["gen_network"])