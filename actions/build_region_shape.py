#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import logging
import re

import geopandas as gpd
import pandas as pd
import pycountry
import shapely.validation as shpval

logger = logging.getLogger(__name__)

from _helpers import configure_logging


if __name__ == "__main__":
    configure_logging(snakemake)

    # Split for each region member the different levels of GADM and remove whitespaces
    region_members = [re.split("\s*->\s*", m) for m in snakemake.params.region_members]

    # Sort based on length; Lists with less elements are for a higher GADM level and are consumed first
    region_members = sorted(region_members, key=len)

    # Read at maximum until level4
    layers_to_read = [f"level{i}" for i in range(5)]

    # Construct dict indicating in which GADM level to search for each region_member
    # and the respective field names to compare against
    search_targets = {layer: [] for layer in layers_to_read}
    for member in region_members:
        member_entry = {
            f"NAME_{position}": element for position, element in enumerate(member)
        }
        search_layer = f"level{len(member)-1}"
        search_targets[search_layer].append(member_entry)

    # Holds all found geometries for region members from different layers
    member_geometries = []

    # Load GADM and search for all region members in the correct layers
    for search_layer, search_target in search_targets.items():
        # Only search layer if we expect a region_member to be in this layer
        if search_target:
            gadm = gpd.read_file(snakemake.input["gadm"], layer=search_layer)

            # Collapse dict into dict of lists for easier use with pandas
            search_target = pd.DataFrame(search_target).to_dict("list")

            # Select all entries with matching columns/values
            gadm = gadm.loc[
                gadm[search_target.keys()].isin(search_target).all(axis="columns")
            ]

            member_geometries.append(gadm)

    if len(member_geometries) == 0:
        logger.error(
            f"No matching entries on GADM found for region '{snakemake.wildcards['region']}'. "
            f"Requested region members: {region_members}."
            f"Check correct and matching spelling with https://gadm.org ."
        )

    member_geometries = pd.concat(member_geometries)

    member_geometries["country"] = member_geometries["NAME_0"]

    if snakemake.wildcards.region == "US-Alaska":
        bbox = (-179.5, 50, -125, 73)
        member_geometries = member_geometries.clip(bbox)
    # elif snakemake.wildcards.region == "BR-Southeast":
    #     bbox = (-54 -34, -28, 6)
    #     member_geometries = member_geometries.clip(bbox)

    # Read EEZs for potential offshore locations and
    # add EEZs for all involved countries
    # (neglecting proximity to specified members for now)
    eez = gpd.read_file(snakemake.input["eez"])

    # Drop entries without ISO_TER1 entry
    # These are mostly small island states + Hawaii
    eez = eez.dropna(axis=0, how="any", subset=["ISO_TER1"])

    # Determine associated country names for EEZs
    eez["country"] = eez["ISO_TER1"].map(
        lambda c: pycountry.countries.get(alpha_3=c).name
    )
    country_renames = {
        "Viet Nam": "Vietnam",
        "Türkiye": "Turkey",
        "Tanzania, United Republic of": "Tanzania",
        "Venezuela, Bolivarian Republic of": "Venezuela",
        "Syrian Arab Republic": "Syria",
        'Iran, Islamic Republic of': 'Iran',
    }
    eez["country"] = eez["country"].replace(country_renames)

    # Drop unnecessary columns
    eez = eez[["country", "geometry"]]

    # Relevant countries
    eez = eez[eez["country"].isin(member_geometries["NAME_0"].unique())]

    # .buffer(...) operation used later is significantly faster and more efficient
    # on exploded shape compared to MultiPolygon.
    # -> explode -> buffer -> union
    member_geometries = member_geometries.explode(ignore_index=True)

    ## Switch to CRS with m[etres] as unit

    # Use Mollweide CRS for estimating offshore distance and adjacency in m[etre]
    # Mollweide is not very accurate at high latitudes, but sufficient for
    # this initial step and large offshore_proximity values.
    # See: https://epsg.io/54009
    # and Usery and Seong (2001), doi:10.1559/152304001782153053
    crs_m = "ESRI:54009"
    crs_org = gadm.crs

    member_geometries = member_geometries.to_crs(crs_m)
    eez = eez.to_crs(crs_m)

    # Merge onshore regions into one
    onshore = member_geometries.union_all()
    onshore = onshore.simplify(0)
    onshore = shpval.make_valid(onshore)

    # For determining offshore regions we only need to buffer the boundary
    # of the onshore region + simplification doesn't hurt/might improve niceiness of shapes
    mgb = member_geometries.copy(deep=True)
    mgb["geometry"] = mgb.boundary
    mgb["geometry"] = mgb.simplify(100)

    # First only consider offshore regions which are adjacent to any onshore region
    # Use buffer to prevent small gaps to overeagerly exclude an offshore region
    eez = eez[eez.geometry.intersects(mgb.buffer(100).union_all())]

    # logger.info(
    #     f"{len(eez)} offshore region(s) found for region '{snakemake.wildcards.region}'."
    # )
    # if len(eez) > 0:
    #     ## Select only offshore locations within <offshore_proximity> m[eters] of
    #     ## an onshore location which is part of the region
    #     ## (=offshore locations accessible from the region under consideration)
    #     # offshore = eez.intersection(mgb.buffer(snakemake.params.offshore_proximity).union_all())

    #     # Combine offshore region into MultiPolygon
    #     offshore = eez.union_all()
    #     offshore = offshore.simplify(0)
    #     offshore = shpval.make_valid(offshore)

    #     # clip offshore region to area within 200NM of onshore region if region
    #     # consists of multiple subregions
    #     if len(snakemake.params.region_members) > 1:
    #         eez_range = 370400  # 200 nautical miles in meters
    #         offshore = onshore.simplify(5e3).buffer(eez_range).intersection(offshore)

    # else:
    #     offshore = None

    offshore = None

    ## Combine resulting region masks and convert back to original CRS for saving
    region_masks = gpd.GeoDataFrame(
        gpd.GeoSeries(
            [onshore, offshore],
            index=["onshore", "offshore"],
            name="geometry",
            crs=crs_m,
        )
    )
    region_masks = region_masks.to_crs(crs_org)

    # Save all region geometries into single file
    region_masks.to_file(snakemake.output["gpkg"], driver="GPKG")

