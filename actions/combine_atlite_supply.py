#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import logging
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from _helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    configure_logging(snakemake)

    # Keys = naming of technologies/snakemake.input.keys() until now
    # Values = Name of the technologies from costs.csv / for downstream starting here
    technology_mapping = {
        "wind_onshore": "onwind",
        "wind_offshore": "offwind",
        "pvplant": "solar-utility",
        "csp_tower": "csp-tower",
    }

    dsl = []  # Used to merge all profile files

    for fp in snakemake.input["profiles"]:
        # determine technology from each file path (could probably be done nicer)
        # example fp: "resources/profiles/Saudi_Arabia_wind_onshore.nc"
        technology = fp.replace(".nc", "").replace(
            f"resources/profiles/{snakemake.wildcards['region']}_", ""
        )
        assert technology in technology_mapping.keys()
        technology = technology_mapping[technology]

        ds = xr.open_dataset(fp)

        # Skipt on empty potentials/profiles to next technology
        if not ds:
            continue

        ds = ds.assign_coords({"technology": technology}).expand_dims("technology")
        ds = ds.transpose("technology", "time", "class")
        dsl.append(ds)

    ds = xr.merge(dsl)

    ## Nice plots, keep for later
    # import hvplot.xarray
    # ds["capacities"].hvplot(kind="bar", by="technology", subplots=True).cols(2)
    # ds["profiles"].mean("time").hvplot(kind="bar", ylabel="Mean CF (p.u.)", by="technology", subplots=True).cols(2)
    # (ds["profiles"].mean("time")*ds["capacities"]/1e6*8780).rename("Available generation (TWh per year)").cumsum(dim="class").hvplot(by="technology")

    ds.to_netcdf(snakemake.output["supply"])

