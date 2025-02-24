# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors, Lukas Franken
#
# SPDX-License-Identifier: MIT


import requests

from os.path import normpath, exists
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.utils import min_version
from datetime import datetime, timedelta

from scripts._helpers import (
    get_scenarios,
    to_total_seconds,
    get_outfiles,
    get_datelist,
    to_quarter,
    get_quarters,
)

HTTP = HTTPRemoteProvider()

min_version("7.7")

conf_file = os.path.join(workflow.current_basedir, "config/config.yaml")
conf_default_file = os.path.join(workflow.current_basedir, "config/config.yaml")
if not exists(conf_file) and exists(conf_default_file):
    copyfile(conf_default_file, conf_file)

configfile: "config/config.yaml"

run = config.get("run", {})
RDIR = run["name"] + "/" if run.get("name") else ""

LOGS = "logs/" + RDIR
BENCHMARKS = "benchmarks/" + RDIR

RESOURCES = "resources/"
RESULTS = "results/"

run = config["run"]
scenarios = get_scenarios(run)

wildcard_constraints:
    # year = r"\d{4}",
    date = r"\d{4}-\d{2}-\d{2}",
    period="[0-9]*",

include: "rules/gather.smk"
include: "rules/testing.smk"
include: "rules/misc.smk"


# Check if the workflow has access to the internet by trying to access the HEAD of specified url
def has_internet_access(url="www.zenodo.org") -> bool:
    import http.client as http_client

    # based on answer and comments from
    # https://stackoverflow.com/a/29854274/11318472
    conn = http_client.HTTPConnection(url, timeout=5)  # need access to zenodo anyway
    try:
        conn.request("HEAD", "/")
        return True
    except:
        return False
    finally:
        conn.close()


"""
if config["enable"].get("retrieve", "auto") == "auto":
    config["enable"]["retrieve"] = has_internet_access()

if config["enable"]["retrieve"] is False:
    print("Datafile downloads disabled in config[retrieve] or no internet access.")


if config["enable"]["retrieve"] and config["enable"].get("retrieve_databundle", True):
    datafiles = [
        "ch_cantons.csv",
        "je-e-21.03.02.xls",
        "eez/World_EEZ_v8_2014.shp",
        "hydro_capacities.csv",
        "naturalearth/ne_10m_admin_0_countries.shp",
        "NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp",
        "nama_10r_3popgdp.tsv.gz",
        "nama_10r_3gdp.tsv.gz",
        "corine/g250_clc06_V18_5.tif",
        "natura/Natura2000_end2015.shp",
        "GEBCO_2014_2D.nc",
    ]

    rule retrieve_databundle:
        output:
            protected(expand("data/bundle/{file}", file=datafiles)),
            # expand("data/bundle/{file}", file=datafiles),
        log:
            LOGS + "retrieve_databundle.log",
        resources:
            mem_mb=1000,
        # retries: 0
        conda:
            "envs/retrieve.yaml"
        script:
            "scripts/retrieve_databundle.py"


rule build_shapes:
    params:
        countries=config["countries"],
    input:
        naturalearth=ancient("data/bundle/naturalearth/ne_10m_admin_0_countries.shp"),
        eez=ancient("data/bundle/eez/World_EEZ_v8_2014.shp"),
        nuts3=ancient("data/bundle/NUTS_2013_60M_SH/data/NUTS_RG_60M_2013.shp"),
        nuts3pop=ancient("data/bundle/nama_10r_3popgdp.tsv.gz"),
        nuts3gdp=ancient("data/bundle/nama_10r_3gdp.tsv.gz"),
        ch_cantons=ancient("data/bundle/ch_cantons.csv"),
        ch_popgdp=ancient("data/bundle/je-e-21.03.02.xls"),
    output:
        country_shapes=RESOURCES + "country_shapes.geojson",
        offshore_shapes=RESOURCES + "offshore_shapes.geojson",
        europe_shape=RESOURCES + "europe_shape.geojson",
        nuts3_shapes=RESOURCES + "nuts3_shapes.geojson",
    log:
        LOGS + "build_shapes.log",
    threads: 1
    resources:
        mem_mb=1500,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/build_shapes.py"


rule base_network:
    params:
        countries=config["countries"],
        # snapshots={k: config["snapshots"][k] for k in ["start", "end", "inclusive"]},
        lines=config["lines"],
        links=config["links"],
        transformers=config["transformers"],
    input:
        eg_buses="data/entsoegridkit/buses.csv",
        eg_lines="data/entsoegridkit/lines.csv",
        eg_links="data/entsoegridkit/links.csv",
        eg_converters="data/entsoegridkit/converters.csv",
        eg_transformers="data/entsoegridkit/transformers.csv",
        parameter_corrections="data/parameter_corrections.yaml",
        links_p_nom="data/links_p_nom.csv",
        links_tyndp="data/links_tyndp.csv",
        # country_shapes=RESOURCES + "country_shapes.geojson",
        country_shapes="data/gb_shape.geojson",
        offshore_shapes=RESOURCES + "offshore_shapes.geojson",
        europe_shape=RESOURCES + "europe_shape.geojson",
    output:
        RESOURCES + "networks/base.nc",
    log:
        LOGS + "base_network.log",
    benchmark:
        BENCHMARKS + "base_network"
    threads: 1
    resources:
        mem_mb=1500,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/base_network.py"
"""


rule build_bus_regions:
    params:
        countries=config["countries"],
    input:
        # country_shapes=RESOURCES + "country_shapes.geojson",
        total_shape="data/gb_shape.geojson",
        offshore_shapes="data/offshore_shapes.geojson",
        # base_network=RESOURCES + "networks/base.nc",
        # base_network="data/ETYS_base.nc",
        base_network="data/base.nc"
        # if config["electricity"]["network_dataset"]
        # else "data/ETYS_base.nc",
    output:
        regions_onshore=RESOURCES + "regions_onshore.geojson",
        regions_offshore=RESOURCES + "regions_offshore.geojson",
    log:
        LOGS + "build_bus_regions.log",
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/build_bus_regions.py"


rule retrieve_wiki_data:
    output:
        wiki_data=RESOURCES + "wiki_data.csv",
    log:
        LOGS + "retrieve_wiki_data.log",
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/retrieve_wiki_data.py"


rule prepare_bmu_data:
    params:
        elexon=config["elexon"],
    input:
        wiki_data=RESOURCES + "wiki_data.csv",
        osuked_ids="data/ids.csv",
        osuked_plant_locations="data/plant-locations.csv",
        manual_bmus="data/manual/bmus.csv",
        carrier_mapper="data/manual/wiki_carrier_mapper.yaml",
    output:
        bmunits_loc=RESOURCES + "bmunits_loc.csv",
    log:
        LOGS + "prepare_bmu_data.log",
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/prepare_bmu_data.py"


rule build_load_weights:
    params:
        load=config["load"],
    input:
        regions_onshore=RESOURCES + "regions_onshore.geojson",
        gsp_regions="data/gsp_geometries.geojson",
        gsp_regions_lookup="data/gsp_gnode_directconnect_region_lookup.csv",
        fes_2021_lw_demandpeaks="data/FES-2021--Leading_the_Way--demandpk-all--gridsupplypoints.csv",
        fes_2021_fs_demandpeaks="data/FES-2021--Falling_Short--demandpk-all--gridsupplypoints.csv",
    output:
        load_weights=RESOURCES + "load_weights.csv",
    log:
        LOGS + "build_load_weights.log",
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/build_load_weights.py"


rule build_bmu_strike_prices:
    input:
        cfd_register="data/cfd/CfD Register 6-Jun-2024.xlsx",
        cfd_bmu_mapping="data/cfd/cfd_to_bm_unit_mapping.csv",
    output:
        strike_prices=RESOURCES + "bmu_strike_prices.csv",
    log:
        LOGS + "build_bmu_strike_prices.log",
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/build_bmu_strike_prices.py"


rule retrieve_live_bmu_data:
    output:
        elexon_bmus=RESOURCES + "live_data/{date}_{period}/elexon_bmus.csv",
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/retrieve_live_bmu_data.py"


rule retrieve_balancing_actions:
    output:
        real_balancing_actions=RESOURCES + "live_data/{date}_{period}/real_balancing_actions.csv",
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/retrieve_balancing_actions.py"


# improve this script
rule retrieve_live_prices:
    output:
        price_stats=RESOURCES + "live_data/{date}_{period}/price_stats.csv",
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/retrieve_live_prices.py"


rule build_dispatchable_costs:
    input:
        bmus_dispatch="data/bmu_operation_winter.csv",
        bmus_data=RESOURCES + "bmunits_loc.csv",
        market_prices="data/market_index_winter.csv",
    output:
        # bmu_cost_estimates=RESOURCES + "bmu_cost_estimates.csv",
    log:
        LOGS + "build_dispatchable_costs.log",
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/build_dispatchable_costs.py"


rule build_redispatch_costs:
    input:
        balancing_actions="data/winter_2024_balancing_actions.csv",
        wholesale_cost_estimates=RESOURCES + "bmu_cost_estimates.csv",
        network=RESOURCES + "networks/gen.nc",
    output:
        redispatch_cost=RESOURCES + "redispatch_cost.csv",
    log:
        LOGS + "build_dispatchable_costs.log",
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/build_redispatch_cost.py"


rule add_generators:
    params:
        elexon=config["elexon"],
        electricity=config["electricity"],
    input:
        base_network="data/base.nc",
        # base_network=RESOURCES + "networks/base.nc"
        # if config["electricity"]["network_dataset"]
        # else "data/ETYS_base.nc",
        gb_shape="data/gb_shape.geojson",
        bmunits_loc=RESOURCES + "bmunits_loc.csv",
        regions_onshore=RESOURCES + "regions_onshore.geojson",
        regions_offshore=RESOURCES + "regions_offshore.geojson",
        carrier_costs="data/manual/carrier_costs.yaml",
        bmu_cost_estimates=RESOURCES + "bmu_cost_estimates.csv",
    output:
        gen_network=RESOURCES + "networks/gen.nc",
        cost_estimated_generators=RESOURCES + "cost_estimated_generators.csv",
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/add_generators.py"


rule prepare_live_network:
    params:
        elexon=config["elexon"],
    input:
        network=RESOURCES + "networks/gen.nc",
        load_weights=RESOURCES + "load_weights.csv",
        elexon_bmus=RESOURCES + "live_data/{date}_{period}/elexon_bmus.csv",
        cost_estimated_generators=RESOURCES + "cost_estimated_generators.csv",
        price_stats=RESOURCES + "live_data/{date}_{period}/price_stats.csv",
    output:
        network=RESOURCES + "live_data/{date}_{period}/network.nc",
    resources:
        mem_mb=1500,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/prepare_live_network.py"


rule retrieve_live_constraint_flows:
    params:
        RESOURCES=RESOURCES,
    output:
        constraint_flows=RESOURCES + "live_data/{date}_{period}/constraint_flows.csv",
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/retrieve_live_constraint_flows.py"


rule simplify_network:
    params:
        simplify_network=config["clustering"]["simplify"],
        aggregation_strategies=config["clustering"]["aggregation_strategies"],
        renewable_carriers=config["electricity"]["renewable_carriers"],
        max_hours=config["electricity"]["max_hours"],
        length_factor=config["lines"]["length_factor"],
        p_max_pu=config["links"]["p_max_pu"],
        costs=config["costs"],
    input:
        network=RESOURCES + "live_data/{date}_{period}/network.nc",
        regions_onshore=RESOURCES + "regions_onshore.geojson",
        regions_offshore=RESOURCES + "regions_offshore.geojson",
        tech_costs="data/costs_2020.csv",
    output:
        network=RESOURCES + "live_data/{date}_{period}/network_s.nc",
        # regions_onshore=RESOURCES + "live_data/{date}_{period}/regions_onshore_s.geojson",
        regions_offshore=RESOURCES + "live_data/{date}_{period}/regions_offshore_s.geojson",
        busmap=RESOURCES + "live_data/{date}_{period}/busmap_s.csv",
        connection_costs=RESOURCES + "live_data/{date}_{period}/connection_costs_s.csv",
    resources:
        mem_mb=1500,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/simplify_network.py"


rule cluster_network:
    params:
        cluster_network=config["clustering"]["cluster_network"],
        aggregation_strategies=config["clustering"]["aggregation_strategies"],
        renewable_carriers=config["electricity"]["renewable_carriers"],
        conventional_carriers=config["electricity"]["conventional_carriers"],
        max_hours=config["electricity"]["max_hours"],
        length_factor=config["lines"]["length_factor"],
        costs=config["costs"],
    input:
        target_regions="data/{layout}_zones.geojson",
        # regions_onshore=RESOURCES + "live_data/{date}_{period}/regions_onshore_s.geojson",
        regions_onshore="data/regions_onshore_s.geojson",
        regions_offshore=RESOURCES + "live_data/{date}_{period}/regions_offshore_s.geojson",
        network=RESOURCES + "live_data/{date}_{period}/network_s.nc",
        tech_costs="data/costs_2020.csv",
    output:
        network=RESOURCES + "live_data/{date}_{period}/network_s_{layout}.nc",
    resources:
        mem_mb=1500,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/cluster_network.py"


rule solve_network:
    params:
        solving=config["solving"]["options"],
        boundaries=config["boundaries"],
    input:
        network=RESOURCES + "live_data/{date}_{period}/network_s_{layout}.nc",
        nodal_network=RESOURCES + "live_data/{date}_{period}/network_s_nodal.nc",
        network_constraints=RESOURCES + "live_data/{date}_{period}/constraint_flows.csv",
    output:
        network=RESOURCES + "live_data/{date}_{period}/network_s_{layout}_solved.nc",
    resources:
        mem_mb=1500,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/solve_network.py"


rule prepare_allowances:
    input:
        dno_regions="data/price_postprocessing/charge_restriction_regions.geojson",
        wholesale_allowances="data/price_postprocessing/Annex_2.xlsx",
        network_allowances="data/price_postprocessing/Annex_3.xlsx",
        policy_allowances="data/price_postprocessing/Annex_4.xlsx",
    output:
        zeroth_order=RESOURCES + "allowances/zeroth_order_allowances_{quarter}_{rate}.csv",
        first_order=RESOURCES + "allowances/first_order_allowances_{quarter}_{rate}.csv",
    resources:
        mem_mb=1500,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/prepare_allowances.py"


rule summarise_period:
    params:
        balancing=config["balancing"]["extra_cost"],
        policy_settings=config["policy_settings"],
    input:
        network_nodal=RESOURCES + "live_data/{date}_{period}/network_s_nodal_solved.nc",
        # regions_nodal=RESOURCES + "live_data/{date}_{period}/regions_onshore_s.geojson",
        regions_nodal="data/regions_onshore_s.geojson",
        # network_fti=RESOURCES + "live_data/{date}_{period}/network_s_fti_solved.nc",
        # regions_fti="data/fti_zones.geojson",
        network_eso=RESOURCES + "live_data/{date}_{period}/network_s_eso_solved.nc",
        regions_eso="data/eso_zones.geojson",
        network_national=RESOURCES + "live_data/{date}_{period}/network_s_national_solved.nc",
        regions_national="data/national_zones.geojson",
        # redispatch_cost=RESOURCES + "redispatch_cost.csv",
        redispatch_cost="data/d_balancing_cost_2022_2023_2024.csv",
        tariffs="data/octopus_12m_fixed_april_2024_v1.csv",
        elexon_demand_profiles="data/elexon_profiling_data_201314.xlsx",
        #allowance_single_standing=(
        #    lambda wildcards:
        #    RESOURCES +
        #    f"allowances/zeroth_order_allowances_{to_quarter(wildcards.date)}_single.csv"
        #),
        #allowance_multi_standing=(
        #    lambda wildcards:
        #    RESOURCES +
        #    f"allowances/zeroth_order_allowances_{to_quarter(wildcards.date)}_multi.csv"
        #),
        #allowance_single_linear=(
        #    lambda wildcards:
        #    RESOURCES +
        #    f"allowances/first_order_allowances_{to_quarter(wildcards.date)}_single.csv"
        #),
        #allowance_multi_linear=(
        #    lambda wildcards:
        #    RESOURCES +
        #    f"allowances/first_order_allowances_{to_quarter(wildcards.date)}_multi.csv"
        #),
        strike_prices=RESOURCES + "bmu_strike_prices.csv",
    output:
        summary=RESULTS + "periods/{date}_{period}.json",
    resources:
        mem_mb=1500,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/summarise_period.py"


rule plot_period:
    input:
        data=ancient(RESULTS + "periods/{date}_{period}.json"),
        regions_nodal="data/regions_onshore_s.geojson",
        regions_fti="data/fti_zones.geojson",
        regions_eso="data/eso_zones.geojson",
        regions_national="data/national_zones.geojson",
    output:
        plot=RESULTS + "plots/{date}_{period}.pdf",
    resources:
        mem_mb=1500,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/plot_period.py"


def apply_aggregation_mode(filelist, mode):
    
    assert mode in ['soft', 'hard', 'missing'], f"Unknown aggregation mode: {mode}, should be 'soft' or 'hard'"
    if mode == 'hard':
        return filelist
    
    import os
    path = Path(filelist[0]).parent

    if mode == 'soft':
        return [fn for fn in filelist if fn.split('/')[-1] in os.listdir(path)]

    elif mode == 'missing':
        return [fn for fn in filelist if not fn.split('/')[-1] in os.listdir(path)]


def get_halfhourly_input(path, template, date, mode):
    import pandas as pd

    infiles = pd.Index([
        path + template.format(date=date, period=period)
        for period in range(1, 49)
        ])

    return apply_aggregation_mode(infiles, mode)


rule aggregate_periods_to_halfhourly:
    params:
        mode=config["aggregation"],
    input:
        lambda wildcards: [
            snakemake.io.ancient(file) for file in
                get_halfhourly_input(
                RESULTS + "periods/",
                "{date}_{period}.json",
                wildcards.date,
                config["aggregation"]
            )],
    output:
        RESULTS + "half-hourly/{date}.json",
    resources:
        mem_mb=1500,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/aggregate_periods_to_halfhourly.py"


def get_daily_input(source_path, template, date, mode):
    import pandas as pd

    infiles = pd.Index([
        source_path + template.format(date=day.strftime('%Y-%m-%d'))
        for day in pd.date_range(
            start=date,
            end=pd.to_datetime(date) + pd.offsets.MonthEnd(1),
            freq='d')
        ])

    return apply_aggregation_mode(infiles, mode)


rule halfhourly_to_daily:
    params:
        mode=config["aggregation"],
    input:
        lambda wildcards: get_daily_input(
            RESULTS + "half-hourly/",
            "{date}.json",
            wildcards.month,
            config["aggregation"]
            ),
    output:
        RESULTS + "daily/{month}.json",
    log:
        LOGS + "halfhourly_to_daily_{month}.log",
    resources:
        mem_mb=1500,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/aggregate_halfhourly_to_daily.py"


def get_monthly_input(source_path, template, date, mode):
    import pandas as pd

    infiles = pd.Index([
        source_path + template.format(date=month.strftime('%Y-%m'))
        for month in pd.date_range(year, str(int(year) + 1), freq='MS')
        ])

    return apply_aggregation_mode(infiles, mode)


rule daily_to_monthly:
    params:
        mode=config["aggregation"],
    input:
        lambda wildcards: get_monthly_input(
            RESULTS + "daily/",
            "{date}.json",
            wildcards.year,
            config["aggregation"]
            ),
        # lambda wildcards: get_daily_input(
        #     RESULTS + "half-hourly/",
        #     "{date}.json",
        #     wildcards.month,
        #     config["aggregation"]
        #     ),
    output:
        RESULTS + "monthly/{year}.json",
    log:
        LOGS + "daily_to_monthly_{year}.log",
    resources:
        mem_mb=1500,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/aggregate_daily_to_monthly.py"


"""
rule aggregate_periods:
    params:
        date=config["scenario"]["aggregate"][0],
    input:
        expand(
            RESULTS + "periods/{date}_{period}.json",
            date=get_datelist(config["scenario"]["aggregate"][0]),
            period=list(range(1, 49))),
    output:
        RESULTS + get_outfiles(config["scenario"]["aggregate"][0]),
        expand(
            RESULTS + "half-hourly/{date}.json",
            date=get_datelist(config["scenario"]["aggregate"][0]),
            ),
    log:
        LOGS + "aggregate_periods.log",
    resources:
        mem_mb=1500,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/aggregate_periods.py"
"""


