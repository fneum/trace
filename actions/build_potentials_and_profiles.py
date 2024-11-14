#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import progressbar as pgb

pgb.streams.wrap_stderr()

import functools
import logging
import os
import time
from tempfile import NamedTemporaryFile
from datetime import datetime

import atlite
import fiona
import geopandas as gpd
import hvplot
import hvplot.pandas
import hvplot.xarray
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib import colors
from rasterio.plot import show
from dask.distributed import Client

logger = logging.getLogger(__name__)

from _helpers import configure_logging


def get_wdpa_layer_name(wdpa_fn, layer_substring):
    """Get layername from file "wdpa_fn" whose name contains "layer_substring"."""
    l = fiona.listlayers(wdpa_fn)
    return [_ for _ in l if layer_substring in _][0]


# In[ ]:


def do_main():
    configure_logging(snakemake)

    nprocesses = int(snakemake.threads)

    if nprocesses > 1:
        client = Client(n_workers=nprocesses, threads_per_worker=1)
        dask_kwargs = {"scheduler": client}
    else:
        client = None
        dask_kwargs = {}

    technology_details = snakemake.params["technology_details"]

    cutout = atlite.Cutout(snakemake.input["cutout"])

    # Full region recognised
    region = gpd.read_file(snakemake.input["region"]).set_index("index")
    # Region in which RES technology can be build
    build_region = region.loc[technology_details["build_in"]]
    build_region = build_region.dropna()

    # Some regions don't have specific build_regions, i.e. offshore regions
    # don't exist for landlocked regions. In this case save an empty DataSet
    if build_region.empty:
        logging.info(
            "No region to build in for this technology. No potentials and profiles generated. Exiting."
        )
        xr.Dataset().to_netcdf(snakemake.output["profiles"])

        # Need to save empty figures otherwise snakemake complains about missing output files
        plt.figure()
        plt.title(
            "Region does not contain the specified type of build_region. No data to plot."
        )
        plt.axis("off")
        plt.savefig(snakemake.output["area_mask"])

        hvplot.save(
            pd.Series().hvplot(
                title="Region does not contain the specified type of build_region. No data to plot."
            ),
            snakemake.output["capacity_factor_map"],
        )
        hvplot.save(
            pd.Series().hvplot(
                title="Region does not contain the specified type of build_region. No data to plot."
            ),
            snakemake.output["potential_map"],
        )

        return  # we don't need to do anything else

    # Initialise excluder with Mollweide CRS
    # Could be made more accurate by choosing a more accurate (local) CRS
    crs_ea = "ESRI:54009"

    # Same CRS for region and exclusion calculation required
    region = region.to_crs(crs_ea)
    build_region = build_region.to_crs(crs_ea)
    excluder = atlite.ExclusionContainer(crs=crs_ea, res=100)
    if technology_details["wdpa"]:
        logger.info("Adding WDPA information...")
        # WDPA layers "points" and "polygons" are treated separately

        # Polygons are excluded with their presented shape
        if "polygons" in technology_details["wdpa"]["layers"]:
            wdpa = gpd.read_file(
                snakemake.input["wdpa"],
                # Only load geometries intersecting the cutout
                bbox=tuple(cutout.bounds),
                # polygon layer, points layer optionally later
                layer=get_wdpa_layer_name(snakemake.input["wdpa"], "polygons"),
            )

            # Limit considered protected areas to configured statuses to include
            wdpa[wdpa["STATUS"].isin(technology_details["wdpa"]["include_status"])]

            # Ensure correct CRS
            wdpa = wdpa.to_crs(crs_ea)

            logging.info(f"Excluding {len(wdpa)} PA polygon areas.")
            # temporary file needed for parallelization
            with NamedTemporaryFile(suffix=".gpkg", delete=False) as f:
                plg_tmp_fn = f.name
            if not wdpa.empty:
                wdpa[["geometry"]].to_file(plg_tmp_fn, driver="GPKG")
                while not os.path.exists(plg_tmp_fn):
                    time.sleep(1)
                excluder.add_geometry(plg_tmp_fn)

        # Points are assumed to be the center of an protected area
        # A circular area around this assumed center with the reported
        # area "REP_AREA" is excluded for points to estimate the protected area
        if "points" in technology_details["wdpa"]["layers"]:
            wdpa = gpd.read_file(
                snakemake.input["wdpa"],
                # Only load geometries intersecting the cutout
                # for the points layer treated here this
                # excludes points outside the bounds for which the
                # modelled area would overlap into the cutout.bounds
                # methodological weakness: accepted here
                bbox=tuple(cutout.bounds),
                layer=get_wdpa_layer_name(snakemake.input["wdpa"], "points"),
            )

            # Limit considered protected areas to configured statuses to include
            wdpa[wdpa["STATUS"].isin(technology_details["wdpa"]["include_status"])]

            # Ensure correct CRS
            wdpa = wdpa.to_crs(crs_ea)

            # Only model areas for points with relevant reported areas (> 1 km)
            wdpa = wdpa[wdpa["REP_AREA"] > 1]

            # Calculate radius in km to buffer around central point
            wdpa["buffer_radius"] = np.sqrt(wdpa["REP_AREA"] / np.pi)

            wdpa = wdpa.set_geometry(wdpa["geometry"].buffer(wdpa["buffer_radius"]))

            logging.info(f"Excluding {len(wdpa)} PA point areas.")
            # temporary file needed for parallelization
            with NamedTemporaryFile(suffix=".gpkg", delete=False) as f:
                pts_tmp_fn = f.name
            if not wdpa.empty:
                wdpa[["geometry"]].to_file(pts_tmp_fn, driver="GPKG")
                while not os.path.exists(pts_tmp_fn):
                    time.sleep(1)
                excluder.add_geometry(pts_tmp_fn)

    if technology_details["wdpa_marine"]:
        logger.info("Adding WDPA marine information...")
        # WDPA layers "points" and "polygons" are treated separately

        # Polygons are excluded with their presented shape
        if "polygons" in technology_details["wdpa_marine"]["layers"]:
            wdpa_marine = gpd.read_file(
                snakemake.input["wdpa_marine"],
                # Only load geometries intersecting the cutout
                bbox=tuple(cutout.bounds),
                # polygon layer, points layer optionally later
                layer=get_wdpa_layer_name(snakemake.input["wdpa_marine"], "polygons"),
            )

            # Limit considered protected areas to configured statuses to include
            wdpa_marine[
                wdpa_marine["STATUS"].isin(
                    technology_details["wdpa_marine"]["include_status"]
                )
            ]

            # Ensure correct CRS
            wdpa_marine = wdpa_marine.to_crs(crs_ea)

            logging.info(f"Excluding {len(wdpa_marine)} MPA polygon areas.")
            with NamedTemporaryFile(suffix=".gpkg", delete=False) as f:
                pts_marine_tmp_fn = f.name
            if not wdpa_marine.empty:
                wdpa_marine["geometry"].to_file(pts_marine_tmp_fn, driver="GPKG")
                while not os.path.exists(pts_marine_tmp_fn):
                    time.sleep(1)
                excluder.add_geometry(pts_marine_tmp_fn)

        # Points are assumed to be the center of an protected area
        # A circular area around this assumed center with the reported
        # area "REP_AREA" is excluded for points to estimate the protected area
        if "points" in technology_details["wdpa_marine"]["layers"]:
            wdpa_marine = gpd.read_file(
                snakemake.input["wdpa_marine"],
                # Only load geometries intersecting the cutout
                # for the points layer treated here this
                # excludes points outside the bounds for which the
                # modelled area would overlap into the cutout.bounds
                # methodological weakness: accepted here
                bbox=tuple(cutout.bounds),
                layer=get_wdpa_layer_name(snakemake.input["wdpa_marine"], "points"),
            )

            # Limit considered protected areas to configured statuses to include
            wdpa_marine[
                wdpa_marine["STATUS"].isin(
                    technology_details["wdpa_marine"]["include_status"]
                )
            ]

            # Ensure correct CRS
            wdpa_marine = wdpa_marine.to_crs(crs_ea)

            # Only model areas for points with relevant reported areas (> 1 km)
            wdpa_marine = wdpa_marine[wdpa_marine["REP_AREA"] > 1]

            # Calculate radius in km to buffer around central point
            wdpa_marine["buffer_radius"] = np.sqrt(wdpa_marine["REP_AREA"] / np.pi)

            wdpa_marine = wdpa_marine.set_geometry(
                wdpa_marine["geometry"].buffer(wdpa_marine["buffer_radius"])
            )

            logging.info(f"Excluding {len(wdpa_marine)} MPA point areas.")
            with NamedTemporaryFile(suffix=".gpkg", delete=False) as f:
                plg_marine_tmp_fn = f.name
            if not wdpa_marine.empty:
                wdpa_marine["geometry"].to_file(plg_marine_tmp_fn, driver="GPKG")
                while not os.path.exists(plg_marine_tmp_fn):
                    time.sleep(1)
                excluder.add_geometry(plg_marine_tmp_fn)

    if technology_details["shipping_routes"]:
        logger.info("Adding shipping route density information...")
        # Exclude routes with ship densities >= the threshold
        # atlite requires the 'partial'-ised function
        func = functools.partial(
            np.less,
            # Config value is average per hour. Multiply by number of hours for which AIS was observed (Jan 2015 - Feb 2021)
            technology_details["shipping_routes"]["max_density_threshold"]
            * (datetime(2021, 3, 1) - datetime(2015, 1, 1, 0, 0)).total_seconds()
            / (60 * 60),
        )
        excluder.add_raster(
            snakemake.input["shipping_routes"],
            codes=func,
            crs="EPSG:4326",
        )

    if technology_details["gebco"]:
        logger.info("Adding GEBCO height information...")
        # GEBCO ships with this CRS, but it is not always properly
        # recorded/read from the datafiles
        crs_gebco = "EPSG:4326"

        if technology_details["gebco"]["max_depth"]:
            # lambda not supported for atlite + multiprocessing
            # use named function np.greater with partially frozen argument instead
            # and exclude areas where: max_depth > grid cell depth
            func = functools.partial(
                np.greater_equal, technology_details["gebco"]["max_depth"]
            )
            excluder.add_raster(
                snakemake.input["gebco"], codes=func, crs=crs_gebco, nodata=-9999
            )

        if technology_details["gebco"]["max_altitude"]:
            # lambda not supported for atlite + multiprocessing
            # exclude areas where: max_altitude < grid cell depth
            func = functools.partial(
                np.less_equal, technology_details["gebco"]["max_altitude"]
            )
            excluder.add_raster(
                snakemake.input["gebco"], codes=func, crs=crs_gebco, nodata=9999
            )

        # if technology_details["gebco"]["max_slope"]:
        #     func = functools.partial(
        #         np.less_equal, technology_details["gebco"]["max_slope"]
        #     )
        #     # Hand temporary raster file to excl. calculator
        #     excluder.add_raster(
        #         snakemake.input["gebco_slope"], crs="ESRI:54009", codes=func
        #     )

    if technology_details["copernicus"]:
        logger.info("Adding Copernicus Land Cover information...")
        # Eligible land use/cover codes for building RES
        excluder.add_raster(
            snakemake.input["copernicus"],
            crs="EPSG:4326",
            codes=technology_details["copernicus"]["include_codes"],
            invert=True,
        )

        # Buffer distance around these land use/cover codes blocked
        # for RES build-up
        if technology_details["copernicus"]["distancing_codes"]:
            excluder.add_raster(
                snakemake.input["copernicus"],
                crs="EPSG:4326",
                codes=technology_details["copernicus"]["distancing_codes"],
                buffer=technology_details["copernicus"]["distance_to_codes"],
            )

    if technology_details["max_shore_distance"]:
        logger.info("Adding 'max_offshore_distance' constraint...")
        # Include only offshore/eez areas within the configured shore distance
        with NamedTemporaryFile(suffix=".gpkg", delete=False) as f:
            onshore_fn = f.name
        region.loc[["onshore"]].to_file(onshore_fn, driver="GPKG")
        while not os.path.exists(onshore_fn):
            time.sleep(1)
        excluder.add_geometry(
            onshore_fn,
            buffer=technology_details["max_shore_distance"],
            invert=True,
        )

    ## Determine eligible/excluded areas on cutout raster (grid cell level)
    logger.info("Calculating eligible RES build areas...")
    availability = cutout.availabilitymatrix(
        build_region, excluder, nprocesses=nprocesses, disable_progressbar=True
    )
    availability = availability.rename("Grid cell share available for RES")
    availability.attrs["units"] = "p.u."

    # Combine all indices from build_region into one availability
    availability = availability.sum(dim="index")

    logger.info("Determining eligible areas on exclude raster resolution...")
    # Determine eligible/excluded areas on exclude raster (100m)
    masked, transform = atlite.gis.shape_availability(
        build_region["geometry"], excluder
    )

    # atlite v0.2.8: masked is now boolean; convert to np.float64 for compatibility with workflow
    masked = masked.astype(np.float64)

    eligible_share = (
        masked.sum() * excluder.res**2 / build_region.geometry.union_all().area
    )

    ## Plot mask for technology of all eligible areas in the relevant parts of the region
    masked[
        masked == 0
    ] = np.nan  # Otherwise plotting masks is not transparent for plot layers below

    # create colormap with single color (= Green) to plot the mask
    cmap = colors.ListedColormap(["green"])
    fig, ax = plt.subplots(figsize=(16, 12))
    # Region to not build in
    region[~region.index.isin(technology_details["build_in"])].plot(
        ax=ax, edgecolor="red", color="None", hatch="///", zorder=1
    )
    # Region to build in
    build_region.plot(ax=ax, edgecolor="k", fc="lightgray", alpha=0.5, lw=1, zorder=5)
    # Mask for eligible region
    show(
        masked,
        transform=transform,
        interpolation="nearest",
        cmap=cmap,
        ax=ax,
        zorder=20,
    )
    ax.set_title(
        f'{snakemake.wildcards["region"]}: Eligible area (green, {eligible_share * 100:2.2f}%)'
    )
    ax.set_xlabel("Longitude [km]")
    ax.set_ylabel("Latitude [km]")
    fig.savefig(snakemake.output["area_mask"])

    # DataArray representing the area of each cutout grid cell
    logger.info("Calculating cutout areas and installable potentials...")
    area = cutout.grid.to_crs("ESRI:54009").area / 1e6
    area = xr.DataArray(
        area.values.reshape(cutout.shape), [cutout.coords["y"], cutout.coords["x"]]
    )
    area.attrs["units"] = "km^2"

    potential = technology_details["technical_potential"] * availability * area
    potential = potential.rename("Technical potential per grid cell")
    potential.attrs["units"] = "MW"

    ## Calculate each grid cells capacity factor to create quality classes
    logger.info("Calculating Capacity Factor...")
    func = getattr(cutout, technology_details["atlite"].pop("method"))
    capacity_factor = func(
        capacity_factor=True,
        dask_kwargs=dask_kwargs,
        show_progress=False,
        **technology_details["atlite"],
    )

    ## Plot/save technical potential (map) and capacity factor (map)
    plot = potential.hvplot(
        title=f"{potential.name} [{potential.attrs['units']}]", cmap="viridis"
    ) * build_region.to_crs(cutout.crs).dissolve().hvplot(color="None")
    hvplot.save(plot, snakemake.output["potential_map"])

    plot = capacity_factor.drop(["lat", "lon"]).hvplot(
        x="x",
        y="y",
        title="Annual capacity factor per grid cell [p.u.]",
        cmap="inferno",
    ) * build_region.to_crs(cutout.crs).dissolve().hvplot(color="None")
    hvplot.save(plot, snakemake.output["capacity_factor_map"])

    ## Create masks for capacity layouts by binning into "number_quality_classes"
    logger.info("Creating masks for all quality classes...")
    cf_min, cf_max, nbins = (
        capacity_factor.min(),
        capacity_factor.max(),
        technology_details["number_quality_classes"],
    )
    bins = np.linspace(cf_min, cf_max, nbins + 1)

    masks = []
    for i in range(nbins):
        m = np.logical_and(bins[i] <= capacity_factor, capacity_factor < bins[i + 1])

        masks.append(m)

    masks = xr.concat(masks, pd.Index(range(nbins), name="class"))

    ## Calculate capacities and profiles for each quality class
    logger.info("Generating time-series...")
    profile, capacities = func(
        matrix=masks.stack(spatial=["y", "x"]),
        layout=potential,
        index=masks.coords["class"],
        per_unit=True,
        return_capacity=True,
        dask_kwargs=dask_kwargs,
        show_progress=False,
        **technology_details["atlite"],
    )

    ## Merge and save
    logger.info("Saving results...")
    ds = xr.merge([profile.rename("profiles"), capacities.rename("capacities")])
    ds.to_netcdf(snakemake.output["profiles"])
    logger.info("Done.")

    if client is not None:
        client.shutdown()


if __name__ == "__main__":
    # Required so we can "return" prematurely to end the script
    # sys.exit() is unsuitable for usage in workflows
    do_main()
