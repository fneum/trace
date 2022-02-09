# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: 2020-2021 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

from functools import lru_cache

@lru_cache
def calculate_annual_investment(name, r, fn):
    """Calculate the annual investment for an installation given a selected WACC.

    Annual investment for the selected installation 'name' is calculated
    using the given WACC as 'discountrate' and properties for the investment
    read from file (CAPEX, FOM in %, lifetime n in years).

    Parameters:
    -----------
    name : str
        of the installation investment, e.g. "onwind".
        Requires corresponding entry in file 'fn'.
    r : float
        discount rate or WACC used to calculate the annuity of the investment.
        In %.
    fn : str or pathlib.Path
        Filename or path to a pandas-readable csv containing the cost and investment
        specific information (investment, FOM, lifetime) for the installation.

    """

    import logging
    from pathlib import Path

    import pandas as pd

    fn = Path(fn)
    assert fn.exists()

    costs = pd.read_csv(fn, comment="#")

    costs = costs[costs["technology"] == name]

    if costs.empty is True:
        logging.info(f"No cost assumptions found for {name}.")
        return 0

    costs = costs.set_index("parameter")

    return calculate_annuity(
        costs.loc["investment", "value"],
        costs.loc["FOM", "value"],
        costs.loc["lifetime", "value"],
        r,
    )

@lru_cache
def calculate_annuity(invest, fom, lifetime, r):
    """
    Calculate annuity based on EAC.

    invest - investment
    fom    - annual FOM in percentage of investment
    lifetime - lifetime of investment in years
    r      - discount rate in percent
    """

    r = r / 100.0

    annuity_factor = r / (1.0 - 1.0 / (r + 1.0) ** (lifetime))

    return (annuity_factor + fom / 100.0) * invest


@lru_cache
def get_efficiency(process, fn, from_bus=None, to_bus=None):
    """Helper function to read the efficiency of 'process' 
    between 'from_bus' (optional) and 'to_bus' (optional'
    from 'fn'. To be used for file "data/efficiencies.csv ."""

    import pandas as pd

    efficiency = pd.read_csv(fn)
    efficiency = efficiency.query("process==@process")
    if from_bus:
        efficiency = efficiency.query("`from`==@from_bus")
    if to_bus:
        efficiency = efficiency.query("to==@to_bus")

    assert len(efficiency) == 1

    return efficiency["efficiency"].iloc[0]


def configure_logging(snakemake, skip_handlers=False):
    """
    Configure the basic behaviour for the logging module.

    Note: Must only be called once from the __main__ section of a script.

    The setup includes printing log messages to STDERR and to a log file defined
    by either (in priority order): snakemake.log.python, snakemake.log[0] or "logs/{rulename}.log".
    Additional keywords from logging.basicConfig are accepted via the snakemake configuration
    file under snakemake.config.logging.

    Parameters
    ----------
    snakemake : snakemake object
        Your snakemake object containing a snakemake.config and snakemake.log.
    skip_handlers : True | False (default)
        Do (not) skip the default handlers created for redirecting output to STDERR and file.
    """

    import logging
    from pathlib import Path

    kwargs = snakemake.config.get("logging", dict())
    kwargs.setdefault("level", "INFO")

    if skip_handlers is False:
        fallback_path = Path(__file__).parent.joinpath(
            "..", "logs", f"{snakemake.rule}.log"
        )
        logfile = snakemake.log.get(
            "python", snakemake.log[0] if snakemake.log else fallback_path
        )
        kwargs.update(
            {
                "handlers": [
                    # Prefer the 'python' log, otherwise take the first log for each
                    # Snakemake rule
                    logging.FileHandler(logfile),
                    logging.StreamHandler(),
                ],
                "datefmt": "%Y-%m-%d %H:%M:%S",
                "format": "%(asctime)s %(name)-14s %(levelname)-8s %(message)s",
            }
        )
    logging.basicConfig(**kwargs)

@lru_cache
def extract_technology(b):
    """Extract the technology from a bus name 'b' by removing trailing '(exp)' or '(imp)' and other content in same braket.

    Examples
    --------
    > extract_technology("hydrogen (g) storage (exp)")
    'hydrogen (g) storage (exp)'

    > extract_technology("hydrogen (g) storage")
    'hydrogen (g) storage'

    > extract_technology("battery inverter (charging, exp)")
    'battery inverter'

    > extract_technology("battery inverter (charging) (exp)")
    'battery inverter (charging)'

    """
    import re

    return re.sub("\([\w\s,\.]*?(?:exp|imp)\)$", "", b.strip()).strip()

@lru_cache
def get_bus_unit(b, n):
    """Get the 'unit' attribute of bus 'b' from PyPSA network 'n."""

    return n.buses.loc[b]["unit"]

@lru_cache
def extract_unit(b, n):
    """Extract the unit of a bus 'b' in the PyPSA network 'n' based on its carrier.

    Carrier name must be '<something> [unit]'.

    Examples
    --------
    > extract_unit('hydrogen (g) (exp)', network)
    't'

    > extract_unit('electricity (exp)', network)
    'MWh'
    """

    import re

    return re.match("^.*?\s*\[?(\w*)\]?$", n.buses.loc[b]["carrier"]).group(1)
