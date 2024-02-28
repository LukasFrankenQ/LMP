import logging
import linopy
import numpy as np
import pandas as pd
import geopandas as gpd

from functools import reduce

from pypsa.clustering.spatial import get_clustering_from_busmap
from _helpers import configure_logging

logger = logging.getLogger(__name__)
configure_logging(snakemake)


def normed(x):
    return (x / x.sum()).fillna(0.0)


def get_feature_for_hac(n, buses_i=None, feature=None):
    if buses_i is None:
        buses_i = n.buses.index

    if feature is None:
        feature = "solar+onwind-time"

    carriers = feature.split("-")[0].split("+")
    if "offwind" in carriers:
        carriers.remove("offwind")
        carriers = np.append(
            carriers, n.generators.carrier.filter(like="offwind").unique()
        )

    if feature.split("-")[1] == "cap":
        feature_data = pd.DataFrame(index=buses_i, columns=carriers)
        for carrier in carriers:
            gen_i = n.generators.query("carrier == @carrier").index
            attach = (
                n.generators_t.p_max_pu[gen_i]
                .mean()
                .rename(index=n.generators.loc[gen_i].bus)
            )
            feature_data[carrier] = attach

    if feature.split("-")[1] == "time":
        feature_data = pd.DataFrame(columns=buses_i)
        for carrier in carriers:
            gen_i = n.generators.query("carrier == @carrier").index
            attach = n.generators_t.p_max_pu[gen_i].rename(
                columns=n.generators.loc[gen_i].bus
            )
            feature_data = pd.concat([feature_data, attach], axis=0)[buses_i]

        feature_data = feature_data.T
        # timestamp raises error in sklearn >= v1.2:
        feature_data.columns = feature_data.columns.astype(str)

    feature_data = feature_data.fillna(0)

    return feature_data


def distribute_clusters(n, n_clusters, focus_weights=None, solver_name="scip"):
    """
    Determine the number of clusters per country.
    """
    L = (
        n.loads_t.p_set.mean()
        .groupby(n.loads.bus)
        .sum()
        .groupby([n.buses.country, n.buses.sub_network])
        .sum()
        .pipe(normed)
    )

    N = n.buses.groupby(["country", "sub_network"]).size()

    assert (
        n_clusters >= len(N) and n_clusters <= N.sum()
    ), f"Number of clusters must be {len(N)} <= n_clusters <= {N.sum()} for this selection of countries."

    if isinstance(focus_weights, dict):
        total_focus = sum(list(focus_weights.values()))

        assert (
            total_focus <= 1.0
        ), "The sum of focus weights must be less than or equal to 1."

        for country, weight in focus_weights.items():
            L[country] = weight / len(L[country])

        remainder = [
            c not in focus_weights.keys() for c in L.index.get_level_values("country")
        ]
        L[remainder] = L.loc[remainder].pipe(normed) * (1 - total_focus)

        logger.warning("Using custom focus weights for determining number of clusters.")

    assert np.isclose(
        L.sum(), 1.0, rtol=1e-3
    ), f"Country weights L must sum up to 1.0 when distributing clusters. Is {L.sum()}."

    m = linopy.Model()
    clusters = m.add_variables(
        lower=1, upper=N, coords=[L.index], name="n", integer=True
    )
    m.add_constraints(clusters.sum() == n_clusters, name="tot")
    # leave out constant in objective (L * n_clusters) ** 2
    m.objective = (clusters * clusters - 2 * clusters * L * n_clusters).sum()
    if solver_name == "gurobi":
        logging.getLogger("gurobipy").propagate = False
    elif solver_name not in ["scip", "cplex"]:
        logger.info(
            f"The configured solver `{solver_name}` does not support quadratic objectives. Falling back to `scip`."
        )
        solver_name = "scip"
    m.solve(solver_name=solver_name)
    return m.solution["n"].to_series().astype(int)


def busmap_for_n_clusters(
    n,
    n_clusters,
    solver_name,
    focus_weights=None,
    algorithm="kmeans",
    feature=None,
    **algorithm_kwds,
):
    if algorithm == "kmeans":
        algorithm_kwds.setdefault("n_init", 1000)
        algorithm_kwds.setdefault("max_iter", 30000)
        algorithm_kwds.setdefault("tol", 1e-6)
        algorithm_kwds.setdefault("random_state", 0)

    def fix_country_assignment_for_hac(n):
        from scipy.sparse import csgraph

        # overwrite country of nodes that are disconnected from their country-topology
        for country in n.buses.country.unique():
            m = n[n.buses.country == country].copy()

            _, labels = csgraph.connected_components(
                m.adjacency_matrix(), directed=False
            )

            component = pd.Series(labels, index=m.buses.index)
            component_sizes = component.value_counts()

            if len(component_sizes) > 1:
                disconnected_bus = component[
                    component == component_sizes.index[-1]
                ].index[0]

                neighbor_bus = n.lines.query(
                    "bus0 == @disconnected_bus or bus1 == @disconnected_bus"
                ).iloc[0][["bus0", "bus1"]]
                new_country = list(set(n.buses.loc[neighbor_bus].country) - {country})[
                    0
                ]

                logger.info(
                    f"overwriting country `{country}` of bus `{disconnected_bus}` "
                    f"to new country `{new_country}`, because it is disconnected "
                    "from its initial inter-country transmission grid."
                )
                n.buses.at[disconnected_bus, "country"] = new_country
        return n


    if algorithm == "hac":
        feature = get_feature_for_hac(n, buses_i=n.buses.index, feature=feature)
        n = fix_country_assignment_for_hac(n)

    if (algorithm != "hac") and (feature is not None):
        logger.warning(
            f"Keyword argument feature is only valid for algorithm `hac`. "
            f"Given feature `{feature}` will be ignored."
        )

    n.determine_network_topology()

    n_clusters = distribute_clusters(
        n, n_clusters, focus_weights=focus_weights, solver_name=solver_name
    )


def clustering_for_n_clusters(
    n,
    n_clusters,
    custom_busmap=False,
    aggregate_carriers=None,
    line_length_factor=1.25,
    aggregation_strategies=dict(),
    solver_name="scip",
    algorithm="hac",
    feature=None,
    extended_link_costs=0,
    focus_weights=None,
):
    if not isinstance(custom_busmap, pd.Series):
        busmap = busmap_for_n_clusters(
            n, n_clusters, solver_name, focus_weights, algorithm, feature
        )
    else:
        busmap = custom_busmap

    line_strategies = aggregation_strategies.get("lines", dict())
    generator_strategies = aggregation_strategies.get("generators", dict())
    one_port_strategies = aggregation_strategies.get("one_ports", dict())

    clustering = get_clustering_from_busmap(
        n,
        busmap,
        aggregate_generators_weighted=True,
        aggregate_generators_carriers=aggregate_carriers,
        aggregate_one_ports=["Load", "StorageUnit"],
        line_length_factor=line_length_factor,
        line_strategies=line_strategies,
        generator_strategies=generator_strategies,
        one_port_strategies=one_port_strategies,
        scale_link_capital_costs=False,
    )

    if not n.links.empty:
        nc = clustering.network
        nc.links["underwater_fraction"] = (
            n.links.eval("underwater_fraction * length").div(nc.links.length).dropna()
        )
        nc.links["capital_cost"] = nc.links["capital_cost"].add(
            (nc.links.length - n.links.length)
            .clip(lower=0)
            .mul(extended_link_costs)
            .dropna(),
            fill_value=0,
        )

    return clustering


def cluster_regions(busmaps, input=None, output=None):
    busmap = reduce(lambda x, y: x.map(y), busmaps[1:], busmaps[0])

    for which in ("regions_onshore", "regions_offshore"):
        regions = gpd.read_file(getattr(input, which))
        regions = regions.reindex(columns=["name", "geometry"]).set_index("name")
        regions_c = regions.dissolve(busmap)
        regions_c.index.name = "name"
        regions_c = regions_c.reset_index()
        regions_c.to_file(getattr(output, which))