# SPDX-FileCopyrightText: : 2017-2023 The PyPSA-Eur Authors, Lukas Franken
#
# SPDX-License-Identifier: MIT


import requests

from os.path import normpath, exists
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.utils import min_version
from datetime import datetime, timedelta

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
        country_shapes=RESOURCES + "country_shapes.geojson",
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


rule build_bus_regions:
    params:
        countries=config["countries"],
    input:
        country_shapes=RESOURCES + "country_shapes.geojson",
        offshore_shapes=RESOURCES + "offshore_shapes.geojson",
        base_network=RESOURCES + "networks/base.nc",
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


rule prepare_wiki_data:
    output:
        wiki_data=RESOURCES + "wiki_data.csv",
    log:
        LOGS + "prepare_wiki_data.log",
    threads: 1
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/prepare_wiki_data.py"