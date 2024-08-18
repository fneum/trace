#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import logging

import atlite
import geopandas as gpd
from dask.distributed import Client
from _helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    configure_logging(snakemake)

    nprocesses = int(snakemake.threads)
    if nprocesses > 1:
        client = Client(n_workers=nprocesses, threads_per_worker=1)
    else:
        client = None
    dask_kwargs = {"scheduler": client}

    region = gpd.read_file(snakemake.input["gpkg"])

    # Fix: Explode does not work if region has no offshore geometry ("None" for offshore)
    region = region.dropna()

    # .buffer(0) to fix some invalid geometries
    # exploding before buffer significantly increases performance
    region = region.explode(ignore_index=True).buffer(0)

    region = region.union_all()

    # TODO: make sure cutout contains full bounds,
    # might need to substract / add in lower/upper bounds some <=0.25°lat/lon

    cutout = atlite.Cutout(
        snakemake.output["cutout"],
        module="era5",
        bounds=region.bounds,
        time=str(snakemake.params["era5_year"]),
    )

    cutout.prepare(
        ["height", "wind", "influx", "temperature"],
        monthly_requests=False,
        concurrent_requests=True,
        dask_kwargs=dask_kwargs,
        compression=None
    )

