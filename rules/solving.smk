# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later


def all_solved_networks(wildcards):
    return expand(
        "results/{scenario}/{year}/{esc}/{exporter}-{importer}/network.nc",
        scenario=wildcards.scenario,
        year=config["scenarios"][wildcards.scenario]["YEARS"],
        esc=config["scenarios"][wildcards.scenario]["ESCS"],
        exporter=config["scenarios"][wildcards.scenario]["EXPORTERS"],
        importer=config["scenarios"][wildcards.scenario]["IMPORTERS"],
    )


rule solve_scenario:
    input:
        inputs="results/{wildcards.scenario}/inputs.tar",
        networks=all_solved_networks,


def get_attempt(wildcards, attempt):
    return attempt


rule solve_network:
    input:
        network=(
            "resources/networks_ip_as/{scenario}/{year}/{esc}/{from}-{to}/network.nc"
        ),
        additional_components="resources/additional_components.pkl",
    output:
        network="results/{scenario}/{year}/{esc}/{from}-{to}/network.nc",
    group:
        "esc"
    threads: config["solver"]["default"]["threads"]
    params:
        scenario=lambda w: get_scenario(w["scenario"]),
    resources:
        mem_mb=lambda w, attempt: attempt * 8000,
        attempt=get_attempt,
    retries: 4
    log:
        python="logs/{scenario}/{year}/{esc}/{from}-{to}/solve_network.log",
    script:
        "../actions/solve_network.py"


rule backup_scenario:
    input:
        data="data/",
        costs="../technology-data/outputs/",
    output:
        tarchive="results/{scenario}/inputs.tar",
    threads: 1
    log:
        python="logs/{scenario}/backup_run.log",
    script:
        "../actions/backup_run.py"
