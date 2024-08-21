# SPDX-FileCopyrightText: 2023 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

import logging
import pickle

import numpy as np
import pypsa
from _helpers import configure_logging

logger = logging.getLogger(__name__)

from pypsa.optimization.compat import define_constraints, get_var, linexpr


def extra_functionality(n, snapshots):
    """... for solving the network; from pypsa-eur-sec."""
    add_battery_constraints(n)
    add_shipping_constraint(n)

    if n.name == "LOHC shipping":
        add_LOHC_shipping_constraints(n)
        add_LOHC_chemical_constraint(n)


def add_battery_constraints(n):
    """Constraint for battery inverter capacity.

    Assumed battery inverter is bidirectional.
    Enforce same capacity for charging/discharging link.

    From pypsa-eur-sec. (2024-08-19)
    """

    chargers = n.links.query(
        "index.str.startswith('battery inverter') "
        "and bus1.str.startswith('battery') "
        "and p_nom_extendable",
        engine="python",
    ).index
    dischargers = n.links.query(
        "index.str.startswith('battery inverter') "
        "and bus0.str.startswith('battery') "
        "and p_nom_extendable",
        engine="python",
    ).index

    eff = n.links.efficiency[dischargers].values
    lhs = (
        n.model["Link-p_nom"].loc[chargers]
        - n.model["Link-p_nom"].loc[dischargers] * eff
    )

    n.model.add_constraints(lhs == 0, name="Link-charger_ratio")


def add_shipping_constraint(n):
    """Constraint for coupling ship loading and unloading links capacities for improved solving times.

    Assumes loading and unloading links follow consistent naming convention used to automatically
    create them in create_network.py.ipynb.
    """

    loading_links = network.links.loc[
        (network.links.index.str.contains("\\d+\\sloading", regex=True))
        & (network.links["bus0"] == "berth (exp)")
    ].index
    unloading_links = network.links.loc[
        (network.links.index.str.contains("\\d+\\sunloading", regex=True))
        & (network.links["bus1"] == "berth (imp)")
    ].index

    lhs = (
        n.model["Link-p_nom"].loc[loading_links]
        - n.model["Link-p_nom"].loc[unloading_links]
    )

    n.model.add_constraints(lhs == 0, name="Link-shipping_constraint")


def add_LOHC_shipping_constraints(n):
    """Constraint for shipping of LOHC to ensure consistent cargo capacity.

    LOHC requires an additional cargo store per convoy to store used LOHC produced
    during the journey or transported from the importer back to the exporter during
    the inbound journey.
    This constraint ensures that this LOHC (used) store is of the same size as the cargo
    store, i.e. the ships size is determined by the cargo store and the LOHC (used) store
    may not exceed this capacity.
    """

    lohc_stores = n.stores.filter(like="cargo LOHC (used)", axis=0).index
    cargo_stores = lohc_stores.str.replace(
        "cargo LOHC (used)", "cargo (exp)", regex=False
    )

    stores_e_nom = get_var(n, "Store", "e_nom")

    lhs = linexpr(
        (1, stores_e_nom[cargo_stores]), (-1.0, stores_e_nom[lohc_stores].values)
    )

    define_constraints(n, lhs, "=", 0, "Store", "lohc_shipping_constraint")


def add_LOHC_chemical_constraint(n):
    # for DBT: 0.944t LOHC per 1 t of loaded LOHC
    lohc_dbt_share = n.links.loc["LOHC dehydrogenation (imp)", "efficiency2"]

    loaded_stores = n.stores.filter(like="LOHC unloaded", axis=0).index.union(
        n.stores.filter(
            regex="LOHC transport ship convoy \\d+ cargo \\(exp\\)", axis=0
        ).index
    )

    unloaded_stores = n.stores.filter(like="LOHC loaded", axis=0).index.union(
        n.stores.filter(
            regex="LOHC transport ship convoy \\d+ cargo LOHC \\(used\\)", axis=0
        ).index
    )

    generator_p_nom = get_var(n, "Generator", "p_nom")
    stores_e = get_var(n, "Store", "e")

    lhs = linexpr(
        (1, stores_e[unloaded_stores]), (lohc_dbt_share, stores_e[loaded_stores].values)
    ).sum(1)
    lhs += linexpr((-1, generator_p_nom["LOHC chemical (exp)"]))[0]

    define_constraints(
        n, lhs, "<=", 0, "Generator", "lohc_chemical__constraint"
    )


def clean_network(n):
    """
    Some additional housekeeping before the optimisation.
    """

    # Set cost for all components with 0 capital cost/capital_cost to defaults from config
    # to avoid solver shenannigans due to 0 +/- eps
    for components_name in ["links", "stores", "generators"]:
        for p in ["capital_cost", "marginal_cost"]:
            components = getattr(n, components_name)
            idx = components.query(f"{p} == 0").index
            # Set default value
            components.loc[idx, p] = scenario[p]["default"]
            # Add random noise (<10%) to default value to make solution space less flat,
            # should help the solver to find the optimum
            components.loc[idx, p] *= 1 + np.random.randint(100) / 1e3

    # clip very small values in time series
    for df in (
        n.generators_t.p_max_pu,
        n.generators_t.p_min_pu,
        n.links_t.p_max_pu,
        n.links_t.p_min_pu,
    ):
        df.where(df > 0.01, other=0.0, inplace=True)

    n.consistency_check()

    return n


def apply_modifiers(n):
    "Apply modifiers from config.yaml ."

    # Modifiers which are applied here:
    # type(s) of component and the component(s) each modifier affect
    mapping = [
        ("links", "electrolysis", "CAPEX_electrolysis"),
        ("links", "battery inverter", "CAPEX_battery"),
        ("stores", "battery", "CAPEX_battery"),
        ("generators", "pvplant|wind", "CAPEX_RES"),
        ("links", "pipeline", "CAPEX_pipeline"),
        ("links", "methanolisation", "CAPEX_MeOHSynthesis"),
        ("links", "direct air capture", "CAPEX_DAC"),
        ("stores", "hydrogen storage tank", "CAPEX_H2storage"),
    ]
    for components_name, search_string, modifier_name in mapping:
        c = getattr(n, components_name)
        c.loc[c.index.str.contains(search_string), "capital_cost"] *= scenario[
            "modifiers"
        ][modifier_name]

    # Modifiers which are applied here:
    # type(s) of component and the component(s) each modifier affect
    mapping = [
        ("generators", "biogenic co2", "OPEX_bioco2"),
        ("generators", "electricity sold", "WTP_excess_power"),
        ("generators", "heat vent", "WTP_waste_heat"),
    ]
    for components_name, search_string, modifier_name in mapping:
        c = getattr(n, components_name)
        c.loc[c.index.str.contains(search_string), "marginal_cost"] *= scenario[
            "modifiers"
        ][modifier_name]

    # Apply modifier to loads (static and time-dependent equally)
    n.loads["p_set"] *= scenario["modifiers"]["import_demand"]
    n.loads_t["p_set"] *= scenario["modifiers"]["import_demand"]

    return n


def average_every_nhours(n, offset):
    m = n.copy(with_time=False)#

    snapshot_weightings = n.snapshot_weightings.resample(offset).sum()
    m.set_snapshots(snapshot_weightings.index)
    m.snapshot_weightings = snapshot_weightings

    for c in n.iterate_components():
        pnl = getattr(m, c.list_name + "_t")
        for k, df in c.pnl.items():
            if not df.empty:
                pnl[k] = df.resample(offset).mean()

    return m

# In[ ]:


if __name__ == "__main__":
    configure_logging(snakemake)
    scenario = snakemake.params["scenario"]

    for op in snakemake.config["solver"].keys():
        solver_options = snakemake.config["solver"][op].copy()
        solver_name = solver_options.pop("name")

        # Load additional components
        with open(snakemake.input["additional_components"], "rb") as f:
            override_component_attrs = pickle.load(f)
        network = pypsa.Network(override_component_attrs=override_component_attrs)

        network.import_from_netcdf(snakemake.input["network"])
        network = apply_modifiers(network)
        network = clean_network(network)

        offset = snakemake.config["time_resolution"]
        network = average_every_nhours(network, offset)

        if solver_name == "gurobi":
            logging.getLogger("gurobipy").setLevel(logging.CRITICAL)

        logger.info(
            f'Solving network using solver options "{op}":\n' f"{solver_options}"
        )
        logger.info("Starting optimization.")
        status, termination_condition = network.optimize(
            snapshots=network.snapshots,
            extra_functionality=extra_functionality,
            solver_name=solver_name,
            solver_options=solver_options,
            log_fn=snakemake.log["python"],
        )
        logger.info("End of optimisation.")

        if status == "ok":
            network.export_to_netcdf(snakemake.output["network"])
            break
        elif status == "warning" and termination_condition == "suboptimal":
            logger.info(
                "Suboptimal optimisation result. "
                'Saving as "_suboptimal.nc in case you want to use it.'
            )
            network.export_to_netcdf(
                snakemake.output["network"].replace(".nc", "_suboptimal.nc")
            )
        else:
            logger.info(
                f"Optimsation ended with \n"
                f"\tstatus: {status}\n"
                f"\ttermination condition: {termination_condition}\n"
                f"Unoptimised network not saved."
                f"Retrying with different solver config."
            )
