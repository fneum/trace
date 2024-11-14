# SPDX-FileCopyrightText: 2022 Johannes Hampp
#
# SPDX-License-Identifier: GPL-3.0-or-later

# - Rules for determining renewables potentials and time-series - #

import requests
from pathlib import Path
from shutil import unpack_archive
from datetime import datetime, timedelta

rule download_eez:
    message:
        """Trying to download EEZ files. If this fail, manually download from https://www.marineregions.org/download_file.php?name=World_EEZ_v12_20231025_LR.zip and extract them into the 'resources/' folder."""
    params:
        zip="resources/World_EEZ_v12_20231025_LR.zip",
    output:
        gpkg="resources/World_EEZ_v12_20231025_LR/eez_v12_lowres.gpkg",
    run:
        import os
        import requests
        from uuid import uuid4

        name = str(uuid4())[:8]
        org = str(uuid4())[:8]

        response = requests.post(
            "https://www.marineregions.org/download_file.php",
            params={"name": "World_EEZ_v12_20231025_LR.zip"},
            data={
                "name": name,
                "organisation": org,
                "email": f"{name}@{org}.org",
                "country": "Germany",
                "user_category": "academia",
                "purpose_category": "Research",
                "agree": "1",
            },
        )

        with open(params["zip"], "wb") as f:
            f.write(response.content)
        output_folder = Path(params["zip"]).parent
        unpack_archive(params["zip"], output_folder)
        os.remove(params["zip"])


# Downloading Copernicus Global Land Cover for land cover and land use:
# Website: https://land.copernicus.eu/global/products/lc
rule download_land_cover:
    input:
        storage(
            "https://zenodo.org/record/3939050/files/PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        ),
    output:
        "resources/Copernicus_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
    shell:
        "mv {input} {output}"


# Downloading bathymetry data (GEBCO)
# Website: https://www.gebco.net/data_and_products/gridded_bathymetry_data/#global
rule download_gebco:
    input:
        storage(
            "https://www.bodc.ac.uk/data/open_download/gebco/gebco_2021/zip/",
        ),
    output:
        zip="resources/gebco_2021.zip",
        gebco="resources/gebco/GEBCO_2021.nc",
    run:
        shell("mv {input} {output.zip}")
        output_folder = Path(output["gebco"]).parent
        shell("unzip {output.zip} -d {output_folder}")


# Transform the GEBCO dataset into a second dataset
# with georeferenced slope values (in percent) for all locations
# Projection of the output dataset is Global Mollweide (ESRI:54009)
rule calculate_gebco_slope:
    input:
        gebco="resources/gebco/GEBCO_2021.nc",
    output:
        mollweide=temp("resources/gebco/GEBCO_2021_mollweide.nc"),
        slope_inflated=temp("resources/gebco/GEBCO_2021_mollweide_inflated.nc"),
        slope="resources/gebco/GEBCO_2021_slope.nc",
    run:
        shell(
            "gdalwarp -multi -wo NUM_THREADS={threads} -of netCDF -co FORMAT=NC4 -s_srs 'EPSG:4326' -t_srs 'ESRI:54009' {input.gebco} {output.mollweide}"
        )
        shell(
            "gdaldem slope -p -of netCDF -co FORMAT=NC4 {output.mollweide} {output.slope_inflated}"
        )
        shell(
            "gdal_translate -of netCDF -co FORMAT=NC4 -co COMPRESS=DEFLATE -co ZLEVEL=1 {output.slope_inflated} {output.slope}"
        )

# Preparation for WDPA downloads

def check_file_exists(url):
    response = requests.head(url)
    return response.status_code == 200

# Basic pattern where WDPA files can be found
url_pattern = (
    "https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_{bYYYY}_Public_shp.zip"
)

# 3-letter month + 4 digit year for current/previous/next month to test
current_monthyear = datetime.now().strftime("%b%Y")
prev_monthyear = (datetime.now() - timedelta(30)).strftime("%b%Y")
next_monthyear = (datetime.now() + timedelta(30)).strftime("%b%Y")

# Test prioritised: current month -> previous -> next
for bYYYY in [current_monthyear, prev_monthyear, next_monthyear]:
    if check_file_exists(url := url_pattern.format(bYYYY=bYYYY)):
        break
    else:
        # If None of the three URLs are working
        url = False

assert (
    url
), f"No WDPA files found at {url_pattern} for bY='{current_monthyear}, {prev_monthyear}, or {next_monthyear}'"

# Downloading Marine protected area database from WDPA
# extract the main zip and then merge the contained 3 zipped shapefiles
# Website: https://www.protectedplanet.net/en/thematic-areas/marine-protected-areas
rule download_wdpa_marine:
    input:
        storage(
            f"https://d1gam3xoknrgr2.cloudfront.net/current/WDPA_WDOECM_{bYYYY}_Public_marine_shp.zip",
            keep_local=True,
        ),
    output:
        zip="resources/WDPA_WDOECM_marine.zip",
        folder=directory("resources/WDPA_WDOECM_marine"),
        gpkg="resources/WDPA_WDOECM_marine.gpkg",
    run:
        shell("mv {input} {output.zip}")
        shell("unzip {output.zip} -d {output.folder}")
        for i in range(3):
            # vsizip is special driver for directly working with zipped shapefiles in ogr2ogr
            layer_path = (
                f"/vsizip/{output.folder}/WDPA_WDOECM_{bYYYY}_Public_marine_shp_{i}.zip"
            )
            print(f"Adding layer {i+1} of 3 to combined output file.")
            shell("ogr2ogr -f gpkg -update -append {output.gpkg} {layer_path}")


# Downloading protected area database from WDPA
# extract the main zip and then merge the contained 3 zipped shapefiles
# Website: https://www.protectedplanet.net/en/thematic-areas/wdpa
rule download_wdpa:
    input:
        storage(url, keep_local=True),
    output:
        zip="resources/WDPA_shp.zip",
        folder=directory("resources/WDPA"),
        gpkg="resources/WDPA.gpkg",
    run:
        shell("mv {input} {output.zip}")
        shell("unzip {output.zip} -d {output.folder}")
        for i in range(3):
            # vsizip is special driver for directly working with zipped shapefiles in ogr2ogr
            layer_path = f"/vsizip/{output.folder}/WDPA_{bYYYY}_Public_shp_{i}.zip"
            print(f"Adding layer {i+1} of 3 to combined output file.")
            shell("ogr2ogr -f gpkg -update -append {output.gpkg} {layer_path}")


# Downloading global shipping traffic density map
# with ship movement between 2015-2020
# Specific dataset: "Global Ship Density"
# Website: https://datacatalog.worldbank.org/search/dataset/0037580/Global-Shipping-Traffic-Density
rule download_shipping_density:
    output:
        zip="resources/shipdensity_global.zip",
        shipdensity="resources/shipdensity/shipdensity_global.tif",
    run:
        shell(
            "curl 'https://datacatalogapi.worldbank.org/ddhxext/ResourceDownload?resource_unique_id=DR0045406&version_id=2022-06-28T08:26:08.9439370Z' --output {output.zip} --silent"
        )
        output_folder = Path(output["shipdensity"]).parent
        shell("unzip {output.zip} -d {output_folder}")


rule build_region_shape:
    message:
        "Creating region definition (off and onshore) for: {wildcards.region}."
    input:
        gadm="resources/gadm/gadm36_levels.gpkg",
        eez="resources/World_EEZ_v12_20231025_LR/eez_v12_lowres.gpkg",
    output:
        gpkg="resources/regions/{region}.gpkg",
    params:
        region_members=lambda w: config["regions"][w["region"]],
    log:
        python="logs/build_region_shape/{region}.log",
    resources:
        mem_mb=lambda wc, attempt: 32000 * 2 ** (attempt - 1),
    retries:
        2
    script:
        "../actions/build_region_shape.py"


rule build_cutout:
    message:
        "Downloading cutout for: {wildcards.region}."
    input:
        gpkg="resources/regions/{region}.gpkg",
    output:
        cutout="resources/cutouts/{region}.nc",
    params:
        era5_year=config["renewables"]["era5_year"],
    log:
        python="logs/build_cutout/{region}.log",
    retries: 2
    resources:
        mem_mb=32000,
        runtime="6h",
        cdsapi_tokens=1,
        # snakemake --resources cdsapi_tokens=1 to limit to one parallel execution
        # https://stackoverflow.com/questions/51977436/restrict-number-of-jobs-by-a-rule-in-snakemake
    threads: 8
    script:
        "../actions/build_cutout.py"


rule build_potentials_and_profiles:
    message:
        "Determining {wildcards.technology} potentials and generation profiles for {wildcards.region}."
    input:
        copernicus="resources/Copernicus_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        gebco="resources/gebco/GEBCO_2021.nc",
        # gebco_slope="resources/gebco/GEBCO_2021_slope.nc",
        wdpa="resources/WDPA.gpkg",
        wdpa_marine="resources/WDPA_WDOECM_marine.gpkg",
        shipping_routes="resources/shipdensity/shipdensity_global.tif",
        cutout="resources/cutouts/{region}.nc",
        region="resources/regions/{region}.gpkg",
    output:
        profiles="resources/profiles/{region}_{technology}.nc",
        area_mask="figures/{region}/{technology}/area_mask.png",
        capacity_factor_map="figures/{region}/{technology}/capacity_factors.html",
        potential_map="figures/{region}/{technology}/potential_map.html",
    params:
        technology_details=lambda w: config["renewables"][w.technology],
    wildcard_constraints:
        technology="(pvplant|wind_onshore|wind_offshore|csp_tower)",
    threads: 8
    resources:
        mem_mb=lambda wc, attempt: 32000 * 2 ** (attempt - 1),
    retries: 2
    log:
        python="logs/build_potentials_and_profiles/{region}_{technology}.log",
    benchmark:
        "benchmarks/build_potentials_and_profiles/{region}_{technology}.csv"
    script:
        "../actions/build_potentials_and_profiles.py"


# Convert atlite RES supply files into a single file
rule combine_atlite_supply:
    input:
        profiles=expand(
            "resources/profiles/{region}_{technology}.nc",
            technology=["wind_offshore", "wind_onshore", "pvplant"], # also supported: "csp-tower"
            allow_missing=True,
        ),
    output:
        supply="resources/supply_{region}.nc",
    log:
        python="logs/combine_atlite_supply/{region}.log",
    script:
        "../actions/combine_atlite_supply.py"


# Downloading GADM database for country/region shapes which are used to define
# export region extents.
# Website: https://gadm.org/
rule download_gadm:
    input:
        storage(
            "https://geodata.ucdavis.edu/gadm/gadm3.6/gadm36_levels_gpkg.zip",
            keep_local=True,
        ),
    output:
        "resources/gadm/gadm36_levels.gpkg",
    run:
        output_folder = Path(output[0]).parent
        shell("unzip {input} -d {output_folder}")


# Downloading Global LNG Terminals database
# Website: https://globalenergymonitor.org
rule download_lng_terminals:
    output:
        "resources/gem/GEM-GGIT-LNG-Terminals-2024-01.xlsx",
    run:
        import requests
        response = requests.get(
            "https://globalenergymonitor.org/wp-content/uploads/2024/03/GEM-GGIT-LNG-Terminals-2024-01.xlsx",
            headers={'User-Agent': 'Mozilla/5.0'}
        )
        with open(output[0], "wb") as f:
            f.write(response.content)


rule build_distances:
    input:
        lng="resources/gem/GEM-GGIT-LNG-Terminals-2024-01.xlsx",
    output:
        "resources/distances.csv",
    log:
        python="logs/build_distances.log",
    script:
        "../actions/build_distances.py"


rule download_synde:
    input:
        storage(
            "https://zenodo.org/records/6569890/files/resources.zip",
            keep_local=True,
        ),
    output:
        directory("resources/synde"),
        expand(
            "resources/synde/resources/ssp2-2.6/2050/era5_2013/{region}.csv",
            region=["Africa", "Asia", "Europe", "SouthAmerica"]
        ),
    run:
        output_folder = Path(output[0])
        shell("unzip {input} -d {output_folder}")


rule build_demand:
    input:
        synde=expand(
            "resources/synde/resources/ssp2-2.6/2050/era5_2013/{region}.csv",
            region=["Africa", "Asia", "Europe", "SouthAmerica"]
        ),
    output:
        "resources/demand.csv",
    log:
        python="logs/build_demand.log",
    script:
        "../actions/build_demand.py"