import logging
from functools import reduce

import numpy as np
import pandas as pd
import pypsa
import scipy as sp

from _helpers import (
    configure_logging,
    update_p_nom_max,
    set_scenario_config,
    load_costs,
    check_network_consistency,
)
from cluster_network import cluster_regions, clustering_for_n_clusters
from pypsa.clustering.spatial import (
    aggregateoneport,
    busmap_by_stubs,
    get_clustering_from_busmap,
)
from pypsa.io import import_components_from_dataframe, import_series_from_dataframe
from scipy.sparse.csgraph import connected_components, dijkstra

logger = logging.getLogger(__name__)


"""
def simplify_network_to_380(n):
    Fix all lines to a voltage level of 380 kV and remove all transformers.

    The function preserves the transmission capacity for each line while
    updating its voltage level, line type and number of parallel bundles
    (num_parallel).

    Transformers are removed and connected components are moved from
    their starting bus to their ending bus. The corresponding starting
    buses are removed as well.

    target_value = 400.

    logger.info(f"Mapping all network lines onto a single {int(target_value)}kV layer")

    n.buses["v_nom"] = target_value
    print(n.lines.v_nom.value_counts())

    (linetype,) = n.lines.loc[n.lines.v_nom == target_value, "type"].unique()

    linetype = "Al/St 240/40 4-bundle 380.0"
    n.lines["type"] = linetype
    n.lines["v_nom"] = int(target_value)
    print('linetype: ', linetype)
    print(n.line_types.i_nom)
    print(n.line_types.i_nom[linetype])
    n.lines["i_nom"] = n.line_types.i_nom[linetype]
    n.lines["num_parallel"] = n.lines.eval("s_nom / (sqrt(3) * v_nom * i_nom)")

    trafo_map = pd.Series(n.transformers.bus1.values, n.transformers.bus0.values)
    print(n.transformers)
    print('trafo map')
    print(trafo_map)
    print(trafo_map.shape)
    print(trafo_map.index.duplicated(keep="first"))

    trafo_map = trafo_map[~trafo_map.index.duplicated(keep="first")]
    print('after')
    print(trafo_map)
    several_trafo_b = trafo_map.isin(trafo_map.index)
    print('several trafor b')
    print(several_trafo_b)
    print(several_trafo_b.loc[several_trafo_b])
    print('total')
    print(several_trafo_b.sum())
    trafo_map[several_trafo_b] = trafo_map[several_trafo_b].map(trafo_map)
    missing_buses_i = n.buses.index.difference(trafo_map.index)
    print('missing buses: ', missing_buses_i)
    missing = pd.Series(missing_buses_i, missing_buses_i)
    print('missing: ', missing)
    trafo_map = pd.concat([trafo_map, missing])

    print(n.lines.v_nom.value_counts())

    for name, trafo in n.transformers.iterrows():
        n.add("Line", name, **dict(trafo))
        n.transformers.drop(name, inplace=True)

        # print(name)
        # print(trafo)
        # print(dict(trafo))

        # import sys
        # sys.exit()

    for c in n.one_port_components | n.branch_components:

        df = n.df(c)
        for col in df.columns:
            if col.startswith("bus"):
                df[col] = df[col].map(trafo_map)

        print(f'--------------------- after {c} loop -----------------------------------')
        isolated_buses = check_network_consistency(n)
        logger.info(f"A total of {len(isolated_buses)} isolated buses")

    print('--------------------- beofre transformer removal -----------------------------------')
    isolated_buses = check_network_consistency(n)
    logger.info(f"A total of {len(isolated_buses)} isolated buses:\n" + ",".join(isolated_buses))

    n.mremove("Transformer", n.transformers.index)

    print('--------------------- after transformer removal -----------------------------------')
    isolated_buses = check_network_consistency(n)
    logger.info(f"A total of {len(isolated_buses)} isolated buses:\n" + ",".join(isolated_buses))

    n.mremove("Bus", n.buses.index.difference(trafo_map))

    print('--------------------- after buses removal -----------------------------------')
    isolated_buses = check_network_consistency(n)
    logger.info(f"A total of {len(isolated_buses)} isolated buses:\n" + ",".join(isolated_buses))

    return n, pd.Series(n.buses.index, n.buses.index)
    # return n, trafo_map
"""

def simplify_network_to_380(n):
    """
    Fix all lines to a voltage level of 380 kV and remove all transformers.

    The function preserves the transmission capacity for each line while
    updating its voltage level, line type and number of parallel bundles
    (num_parallel).

    Transformers are removed and connected components are moved from
    their starting bus to their ending bus. The corresponding starting
    buses are removed as well.
    """
    target_value = 380.
    logger.warning(f"Setting lines voltage to {target_value} kV and removing transformers.")

    logger.info(f"Mapping all network lines onto a single {target_value} kV layer")

    n.buses["v_nom"] = target_value

    (linetype_380,) = n.lines.loc[n.lines.v_nom == target_value, "type"].unique() #
    n.lines["type"] = linetype_380 #
    n.lines["v_nom"] = target_value
    n.lines["i_nom"] = n.line_types.i_nom[linetype_380] #
    n.lines["num_parallel"] = n.lines.eval("s_nom / (sqrt(3) * v_nom * i_nom)") #

    trafo_map = pd.Series(n.transformers.bus1.values, n.transformers.bus0.values)
    trafo_map = trafo_map[~trafo_map.index.duplicated(keep="first")]
    several_trafo_b = trafo_map.isin(trafo_map.index)
    trafo_map[several_trafo_b] = trafo_map[several_trafo_b].map(trafo_map)
    missing_buses_i = n.buses.index.difference(trafo_map.index)
    missing = pd.Series(missing_buses_i, missing_buses_i)
    trafo_map = pd.concat([trafo_map, missing])

    for c in n.one_port_components | n.branch_components:
        df = n.df(c)
        for col in df.columns:
            if col.startswith("bus"):
                df[col] = df[col].map(trafo_map)

    n.mremove("Transformer", n.transformers.index)
    n.mremove("Bus", n.buses.index.difference(trafo_map))

    return n, trafo_map



def _prepare_connection_costs_per_link(n, costs, renewable_carriers, length_factor):
    if n.links.empty:
        return {}

    return {
        tech: (
            n.links.length
            * length_factor
            # * (
            #     n.links.underwater_fraction
            #     * costs.at[tech + "-connection-submarine", "capital_cost"]
            #     + (1.0 - n.links.underwater_fraction)
            #     * costs.at[tech + "-connection-underground", "capital_cost"]
            # )
        )
        for tech in renewable_carriers
        if tech.startswith("offwind")
    }

def _compute_connection_costs_to_bus(
    n,
    busmap,
    costs,
    renewable_carriers,
    length_factor,
    connection_costs_per_link=None,
    buses=None,
):
    if connection_costs_per_link is None:
        connection_costs_per_link = _prepare_connection_costs_per_link(
            n, costs, renewable_carriers, length_factor
        )

    if buses is None:
        buses = busmap.index[busmap.index != busmap.values]

    connection_costs_to_bus = pd.DataFrame(index=buses)

    for tech in connection_costs_per_link:
        adj = n.adjacency_matrix(
            weights=pd.concat(
                dict(
                    Link=connection_costs_per_link[tech].reindex(n.links.index),
                    Line=pd.Series(0.0, n.lines.index),
                )
            )
        )

        costs_between_buses = dijkstra(
            adj, directed=False, indices=n.buses.index.get_indexer(buses)
        )
        connection_costs_to_bus[tech] = costs_between_buses[
            np.arange(len(buses)), n.buses.index.get_indexer(busmap.loc[buses])
        ]

    return connection_costs_to_bus


def _adjust_capital_costs_using_connection_costs(n, connection_costs_to_bus, output):
    connection_costs = {}
    for tech in connection_costs_to_bus:
        tech_b = n.generators.carrier == tech
        costs = (
            n.generators.loc[tech_b, "bus"]
            .map(connection_costs_to_bus[tech])
            .loc[lambda s: s > 0]
        )
        if not costs.empty:
            n.generators.loc[costs.index, "capital_cost"] += costs
            logger.info(
                "Displacing {} generator(s) and adding connection costs to capital_costs: {} ".format(
                    tech,
                    ", ".join(
                        "{:.0f} Eur/MW/a for `{}`".format(d, b)
                        for b, d in costs.items()
                    ),
                )
            )
            connection_costs[tech] = costs
    pd.DataFrame(connection_costs).to_csv(output.connection_costs)


def _aggregate_and_move_components(
    n,
    busmap,
    connection_costs_to_bus,
    output,
    aggregate_one_ports={"Load", "StorageUnit"},
    aggregation_strategies=dict(),
    exclude_carriers=None,
):
    def replace_components(n, c, df, pnl):
        n.mremove(c, n.df(c).index)

        import_components_from_dataframe(n, df, c)
        for attr, df in pnl.items():
            if not df.empty:
                import_series_from_dataframe(n, df, c, attr)

    _adjust_capital_costs_using_connection_costs(n, connection_costs_to_bus, output)

    generator_strategies = aggregation_strategies["generators"]

    carriers = set(n.generators.carrier) - set(exclude_carriers)
    generators, generators_pnl = aggregateoneport(
        n,
        busmap,
        "Generator",
        carriers=carriers,
        custom_strategies=generator_strategies,
    )

    replace_components(n, "Generator", generators, generators_pnl)

    for one_port in aggregate_one_ports:
        df, pnl = aggregateoneport(n, busmap, component=one_port)
        replace_components(n, one_port, df, pnl)

    buses_to_del = n.buses.index.difference(busmap)
    n.mremove("Bus", buses_to_del)
    for c in n.branch_components:
        df = n.df(c)
        n.mremove(c, df.index[df.bus0.isin(buses_to_del) | df.bus1.isin(buses_to_del)])


def simplify_links(
    n,
    costs,
    renewables,
    length_factor,
    p_max_pu,
    exclude_carriers,
    output,
    aggregation_strategies=dict(),
):
    ## Complex multi-node links are folded into end-points
    logger.info("Simplifying connected link components")

    if n.links.empty:
        return n, n.buses.index.to_series()

    # Determine connected link components, ignore all links but DC
    adjacency_matrix = n.adjacency_matrix(
        branch_components=["Link"],
        weights=dict(Link=(n.links.carrier == "DC").astype(float)),
    )

    _, labels = connected_components(adjacency_matrix, directed=False)
    labels = pd.Series(labels, n.buses.index)

    G = n.graph()

    def split_links(nodes):
        nodes = frozenset(nodes)

        seen = set()
        supernodes = {m for m in nodes if len(G.adj[m]) > 2 or (set(G.adj[m]) - nodes)}

        for u in supernodes:
            for m, ls in G.adj[u].items():
                if m not in nodes or m in seen:
                    continue

                buses = [u, m]
                links = [list(ls)]  # [name for name in ls]]

                while m not in (supernodes | seen):
                    seen.add(m)
                    for m2, ls in G.adj[m].items():
                        if m2 in seen or m2 == u:
                            continue
                        buses.append(m2)
                        links.append(list(ls))  # [name for name in ls])
                        break
                    else:
                        # stub
                        break
                    m = m2
                if m != u:
                    yield pd.Index((u, m)), buses, links
            seen.add(u)

    busmap = n.buses.index.to_series()

    connection_costs_per_link = _prepare_connection_costs_per_link(
        n, costs, renewables, length_factor
    )
    connection_costs_to_bus = pd.DataFrame(
        0.0, index=n.buses.index, columns=list(connection_costs_per_link)
    )

    for lbl in labels.value_counts().loc[lambda s: s > 2].index:
        for b, buses, links in split_links(labels.index[labels == lbl]):
            if len(buses) <= 2:
                continue

            logger.debug("nodes = {}".format(labels.index[labels == lbl]))
            logger.debug("b = {}\nbuses = {}\nlinks = {}".format(b, buses, links))

            m = sp.spatial.distance_matrix(
                n.buses.loc[b, ["x", "y"]], n.buses.loc[buses[1:-1], ["x", "y"]]
            )
            busmap.loc[buses] = b[np.r_[0, m.argmin(axis=0), 1]]
            connection_costs_to_bus.loc[buses] += _compute_connection_costs_to_bus(
                n,
                busmap,
                costs,
                renewables,
                length_factor,
                connection_costs_per_link,
                buses,
            )

            all_links = [i for _, i in sum(links, [])]

            lengths = n.links.loc[all_links, "length"]
            name = lengths.idxmax() + "+{}".format(len(links) - 1)
            params = dict(
                carrier="DC",
                bus0=b[0],
                bus1=b[1],
                length=sum(
                    n.links.loc[[i for _, i in l], "length"].mean() for l in links
                ),
                p_nom=min(n.links.loc[[i for _, i in l], "p_nom"].sum() for l in links),
                # underwater_fraction=sum(
                #     lengths
                #     / lengths.sum()
                #     * n.links.loc[all_links, "underwater_fraction"]
                # ),
                p_max_pu=p_max_pu,
                p_min_pu=-p_max_pu,
                underground=False,
                under_construction=False,
            )

            logger.info(
                "Joining the links {} connecting the buses {} to simple link {}".format(
                    ", ".join(all_links), ", ".join(buses), name
                )
            )

            n.mremove("Link", all_links)

            static_attrs = n.components["Link"]["attrs"].loc[lambda df: df.static]
            for attr, default in static_attrs.default.items():
                params.setdefault(attr, default)
            n.links.loc[name] = pd.Series(params)

            # n.add("Link", **params)

    logger.debug("Collecting all components using the busmap")

    _aggregate_and_move_components(
        n,
        busmap,
        connection_costs_to_bus,
        output,
        aggregation_strategies=aggregation_strategies,
        exclude_carriers=exclude_carriers,
    )
    return n, busmap


def remove_stubs(
    n,
    costs,
    renewable_carriers,
    length_factor,
    simplify_network,
    output,
    aggregation_strategies=dict(),
):
    logger.info("Removing stubs")

    across_borders = simplify_network["remove_stubs_across_borders"]
    matching_attrs = [] if across_borders else ["country"]
    busmap = busmap_by_stubs(n, matching_attrs)

    connection_costs_to_bus = _compute_connection_costs_to_bus(
        n, busmap, costs, renewable_carriers, length_factor
    )

    _aggregate_and_move_components(
        n,
        busmap,
        connection_costs_to_bus,
        output,
        aggregation_strategies=aggregation_strategies,
        exclude_carriers=simplify_network["exclude_carriers"],
    )

    return n, busmap


def aggregate_to_substations(n, aggregation_strategies=dict(), buses_i=None):
    # can be used to aggregate a selection of buses to electrically closest neighbors
    # if no buses are given, nodes that are no substations or without offshore connection are aggregated

    if buses_i is None:
        logger.info(
            "Aggregating buses that are no substations or have no valid offshore connection"
        )
        buses_i = list(set(n.buses.index) - set(n.generators.bus) - set(n.loads.bus))

    weight = pd.concat(
        {
            "Line": n.lines.length / n.lines.s_nom.clip(1e-3),
            "Link": n.links.length / n.links.p_nom.clip(1e-3),
        }
    )

    adj = n.adjacency_matrix(branch_components=["Line", "Link"], weights=weight)

    bus_indexer = n.buses.index.get_indexer(buses_i)
    dist = pd.DataFrame(
        dijkstra(adj, directed=False, indices=bus_indexer), buses_i, n.buses.index
    )

    dist[buses_i] = (
        np.inf
    )  # bus in buses_i should not be assigned to different bus in buses_i

    for c in n.buses.country.unique():
        incountry_b = n.buses.country == c
        dist.loc[incountry_b, ~incountry_b] = np.inf

    busmap = n.buses.index.to_series()
    busmap.loc[buses_i] = dist.idxmin(1)

    line_strategies = aggregation_strategies.get("lines", dict())
    generator_strategies = aggregation_strategies.get("generators", dict())
    one_port_strategies = aggregation_strategies.get("one_ports", dict())

    clustering = get_clustering_from_busmap(
        n,
        busmap,
        aggregate_generators_weighted=True,
        aggregate_generators_carriers=None,
        aggregate_one_ports=["Load", "StorageUnit"],
        line_length_factor=1.0,
        line_strategies=line_strategies,
        generator_strategies=generator_strategies,
        one_port_strategies=one_port_strategies,
        scale_link_capital_costs=False,
    )
    return clustering.network, busmap


def cluster(
    n,
    n_clusters,
    focus_weights,
    solver_name,
    algorithm="hac",
    feature=None,
    aggregation_strategies=dict(),
):
    logger.info(f"Clustering to {n_clusters} buses")

    clustering = clustering_for_n_clusters(
        n,
        n_clusters,
        custom_busmap=False,
        aggregation_strategies=aggregation_strategies,
        solver_name=solver_name,
        algorithm=algorithm,
        feature=feature,
        focus_weights=focus_weights,
    )

    return clustering.network, clustering.busmap


if __name__ == "__main__":

    configure_logging(snakemake)
    set_scenario_config(snakemake)

    params = snakemake.params
    solver_name = snakemake.config["solving"]["solver"]["name"]

    n = pypsa.Network(snakemake.input.network)

    Nyears = n.snapshot_weightings.objective.sum() / 8760

    # remove integer outputs for compatibility with PyPSA v0.26.0
    n.generators.drop("n_mod", axis=1, inplace=True, errors="ignore")

    n, trafo_map = simplify_network_to_380(n)

    technology_costs = load_costs(
        snakemake.input.tech_costs,
        params.costs,
        params.max_hours,
        Nyears,
    )

    n, simplify_links_map = simplify_links(
        n,
        technology_costs,
        params.renewable_carriers,
        params.length_factor,
        params.p_max_pu,
        params.simplify_network["exclude_carriers"],
        snakemake.output,
        params.aggregation_strategies,
    )

    # busmaps = [trafo_map]
    busmaps = [trafo_map, simplify_links_map]

    if params.simplify_network["remove_stubs"]:
        n, stub_map = remove_stubs(
            n,
            technology_costs,
            params.renewable_carriers,
            params.length_factor,
            params.simplify_network,
            snakemake.output,
            aggregation_strategies=params.aggregation_strategies,
        )
        busmaps.append(stub_map)

    if params.simplify_network["to_substations"]:
        n, substation_map = aggregate_to_substations(n, params.aggregation_strategies)
        busmaps.append(substation_map)

    # treatment of outliers (nodes without a profile for considered carrier):
    # all nodes that have no profile of the given carrier are being aggregated to closest neighbor
    if params.simplify_network["algorithm"] == "hac":
        carriers = params.simplify_network["feature"].split("-")[0].split("+")
        for carrier in carriers:
            buses_i = list(
                set(n.buses.index) - set(n.generators.query("carrier == @carrier").bus)
            )
            logger.info(
                f"clustering preparation (hac): aggregating {len(buses_i)} buses of type {carrier}."
            )
            n, busmap_hac = aggregate_to_substations(
                n, params.aggregation_strategies, buses_i
            )
            busmaps.append(busmap_hac)


    if "":
        n, cluster_map = cluster(
            n,
            int(snakemake.wildcards.simpl),
            params.focus_weights,
            solver_name,
            params.simplify_network["algorithm"],
            params.simplify_network["feature"],
            params.aggregation_strategies,
        )
        busmaps.append(cluster_map)

    # some entries in n.buses are not updated in previous functions, therefore can be wrong. as they are not needed
    # and are lost when clustering (for example with the simpl wildcard), we remove them for consistency:
    remove = [
        "symbol",
        "tags",
        "under_construction",
        "onshore_bus",
        "substation_lv",
        "substation_off",
        "geometry",
        "underground",
    ]
    n.buses.drop(remove, axis=1, inplace=True, errors="ignore")

    n.lines.drop(remove, axis=1, errors="ignore", inplace=True)

    update_p_nom_max(n)

    n.meta = dict(snakemake.config, **dict(wildcards=dict(snakemake.wildcards)))
    n.export_to_netcdf(snakemake.output.network)

    busmap_s = reduce(lambda x, y: x.map(y), busmaps[1:], busmaps[0])
    busmap_s.to_csv(snakemake.output.busmap)

    cluster_regions(busmaps, snakemake.input, snakemake.output)