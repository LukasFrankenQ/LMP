# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT

localrules:
    all,
    gather_balancing,
    gather_bmu_data,
    gather_live_prices,


rule all:
    input:
        expand(
            RESOURCES + "live_data/{date}_{period}/maps.pdf",
            **config["scenario"]
        )

rule gather_balancing:
    input:
        expand(
            RESOURCES + "live_data/{date}_{period}/real_balancing_actions.csv",
            **config["scenario"]
        )

rule gather_bmu_data:
    input:
        expand(
            RESOURCES + "live_data/{date}_{period}/elexon_bmus.csv",
            **config["scenario"]
        )

rule gather_live_prices:
    input:
        expand(
            RESOURCES + "live_data/{date}_{period}/price_stats.csv",
            **config["scenario"]
        )
