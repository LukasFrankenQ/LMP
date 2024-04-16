# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT

from scripts._helpers import process_scenarios

localrules:
    all,
    gather_balancing,
    gather_bmu_data,
    gather_live_prices,


rule all:
    input:
        expand(
            RESOURCES + "live_data/{date}_{period}/maps.pdf",
            **process_scenarios(config["scenario"])
        )

rule gather_summaries:
    input:
        expand(
            RESULTS + "half-hourly/{date}_{period}.json",
            **process_scenarios(config["scenario"])
        )

rule gather_balancing:
    input:
        expand(
            RESOURCES + "live_data/{date}_{period}/real_balancing_actions.csv",
            **process_scenarios(config["scenario"])
        )

rule gather_bmu_data:
    input:
        expand(
            RESOURCES + "live_data/{date}_{period}/elexon_bmus.csv",
            **process_scenarios(config["scenario"])
        )

rule gather_live_prices:
    input:
        expand(
            RESOURCES + "live_data/{date}_{period}/price_stats.csv",
            **process_scenarios(config["scenario"])
        )

rule gather_constraint_flows:
    input:
        expand(
            RESOURCES + "live_data/{date}_{period}/constraint_flows.csv",
            **process_scenarios(config["scenario"])
        )

rule gather_model_plots:
    params:
        scenario=config["scenario"],
        RESOURCES=RESOURCES,
    input:
        networks=expand(
            RESOURCES + "live_data/{date}_{period}/network_s_{layout}_solved.nc",
            **process_scenarios(config["scenario"])
        ),
        # eso_networks=expand(
        #     RESOURCES + "live_data/{date}_{period}/network_s_eso_solved.nc",
        #     **process_scenarios(config["scenario"])
        # ),
        # fti_network=expand(
        #     RESOURCES + "live_data/{date}_{period}/network_s_fti_solved.nc",
        #     **process_scenarios(config["scenario"])
        # ),
        # national_network=expand(
        # RESOURCES + "live_data/{date}_{period}/network_s_national_solved.nc",
        #     **process_scenarios(config["scenario"])
        # )
        real_prices=expand(
            RESOURCES + "live_data/{date}_{period}/price_stats.csv",
            **process_scenarios(config["scenario"])
        )
    output:
        prices=RESOURCES + "plots/prices.pdf",
        # generation=RESOURCES + "plots/generation.pdf",
    log:
        RESOURCES + "logs/gather_model_plots.log"
    resources:
        mem_mb=1000,
    conda:
        "envs/environment.yaml"
    script:
        "scripts/gather_model_plots.py"
