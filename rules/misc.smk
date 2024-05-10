# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT

"""
Unrelated to the main workflow, provides rules that are used to prevent more accessible
data formats for some data inputs
"""


rule build_demand_totals:
    input:
        elexon_demand_profiles="data/elexon_profiling_data_201314.xlsx",
    output:
        "data/demand_totals.json",
    conda:
        "../envs/environment.yaml"
    log:
        "../logs/build_demand_totals.log"
    script:
        "../scripts/build_demand_totals.py"
