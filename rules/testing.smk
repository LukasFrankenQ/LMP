# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT


localrules:
    daily_to_monthly,


import os
from pathlib import Path


def get_available(wildcards):
    date = wildcards.date
    path = Path.cwd() / 'results' / 'periods'

    return [path / fn for fn in os.listdir(path) if date in fn]


# soft referring to the aggregation of available data, instead of forcing
# more files to be created.
rule soft_aggregate:
    input:
        get_available
    output:
        
        'results/summaries/{date}.csv'


def get_available_months(wildcards):

    year = wildcards.year
    # path = Path.cwd() / 'results' / 'daily'
    path = RESULTS / 'daily'

    return [path / fn for fn in os.listdir(path) if year in fn]


# daily to monthly is SOFT, i.e. only gathers for available months, instead of forcing
# more files to be created.
rule daily_to_monthly:
    input:
        get_available_months
    output:
        monthly_aggregate='results/monthly/{year}.json'
    log:
        '../logs/daily_to_monthly_{year}.log'
    resources:
        mem_mb=1500,
    conda:
        '../envs/environment.yaml'
    script:
        '../scripts/daily_to_monthly.py'


rule build_dno_overlap_matrices:
    input:
        dno_regions="data/price_postprocessing/charge_restriction_regions.geojson",
        nodal_network=RESOURCES + "live_data/2024-03-01_1/network_s_nodal.nc",
        nodal_regions="data/regions_onshore_s.geojson",
        fti_network=RESOURCES + "live_data/2024-03-01_1/network_s_fti.nc",
        fti_regions="data/fti_zones.geojson",
        eso_network=RESOURCES + "live_data/2024-03-01_1/network_s_eso.nc",
        eso_regions="data/eso_zones.geojson",
    output:
        nodal=RESOURCES + 'overlap_nodal.csv',
        fti=RESOURCES + 'overlap_fti.csv',
        eso=RESOURCES + 'overlap_eso.csv',
    log:
        '../logs/build_overlap_matrices.log',
    resources:
        mem_mb=1500,
    conda:
        '../envs/environment.yaml'
    script:
        '../scripts/build_dno_overlap_matrices.py'


rule build_dno_load_weights:
    params:
        load=config["load"],
    input:
        dno_regions="data/price_postprocessing/charge_restriction_regions.geojson",
        gsp_regions="data/gsp_geometries.geojson",
        gsp_regions_lookup="data/gsp_gnode_directconnect_region_lookup.csv",
        fes_2021_lw_demandpeaks="data/FES-2021--Leading_the_Way--demandpk-all--gridsupplypoints.csv",
        fes_2021_fs_demandpeaks="data/FES-2021--Falling_Short--demandpk-all--gridsupplypoints.csv",
    output:
        RESOURCES + "dno_load_weights.csv",
    log:
        '../logs/build_dno_load_weights.log',
    resources:
        mem_mb=1500,
    conda:
        '../envs/environment.yaml'
    script:
        '../scripts/build_dno_load_weights.py'


def get_year_quarters(year):
    if str(year) == '2024':
        months = [1, 4]
    else:
        months = [1, 4, 7, 10]

    return [f'{year}-{str(m).zfill(2)}' for m in months]


rule build_monthly_consumer_prices:
    params:
        year=config["scenario"]["year"]
    input:
        monthly_results='results/monthly/{year}.json',
        overlap_nodal=RESOURCES + 'overlap_nodal.csv',
        overlap_fti=RESOURCES + 'overlap_fti.csv',
        overlap_eso=RESOURCES + 'overlap_eso.csv',
        policy_allowance="data/price_postprocessing/Annex_4.xlsx",
        dno_load_weights=RESOURCES + "dno_load_weights.csv",
        zeroth_order=expand(
            RESOURCES + "allowances/zeroth_order_allowances_{quarter}_{rate}.csv",
            quarter=get_year_quarters(config["scenario"]["year"][0]), rate=["single", "multi"]
        ),
        first_order=expand(
            RESOURCES + "allowances/first_order_allowances_{quarter}_{rate}.csv",
            quarter=get_year_quarters(config["scenario"]["year"][0]), rate=["single", "multi"]
        ),
    output:
        monthly_results='results/summaries/{year}.json',
    log:
        '../logs/build_monthly_consumer_prices_{year}.log',
    resources:
        mem_mb=1500,
    conda:
        '../envs/environment.yaml'
    script:
        '../scripts/build_monthly_consumer_prices.py'
