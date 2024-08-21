import logging

from _helpers import configure_logging


import pandas as pd
from searoute import searoute
from haversine import haversine
from itertools import product

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    configure_logging(snakemake)

    lng = (
        pd.read_excel(snakemake.input.lng, na_values=["--"])
        .sort_values(by="CapacityInMtpa", ascending=False)
        .fillna({"CapacityInMtpa": 1})
    )

    destinations = [
        "T0557",  # South Hook LNG terminal, UK
        "T0492",  # Gate LNG terminal, Netherlands
        "T0498",  # Swinoujscie LNG terminal, Poland
        "T0462",  # Revithoussa LNG terminal, Greece
        "T0466",  # Adriatic LNG terminal, Italy
        "T0522",  # Barcelona LNG terminal, Spain
        "T0500",  # Sines LNG terminal, Portugal
    ]

    origins = {
        "AR-South": "T0832",
        "AR-North": "T0540",
        "BR-Northeast": "T0935",
        "BR-Southeast": "T0546",
        "AU-West": "T0341",
        "AU-East": "T0333",
        "CA-East": "T0349",
        "CL": "T0560",
        "CN-Northeast": "T0409",
        "CN-Southeast": "T0384",
        "DZ": "T0243",
        "EG": "T0621",
        "LY": "T0255",
        "MA": "T0256",
        "EH": "T0530",
        "NA": "T0622",
        "OM": "T0593",
        "SA": "T0591",
        "TN": "T0244",
        "TR": "T0543",
        "UA": "T0561",
        "US-Northeast": "T0217",
        "US-Southeast": "T0221",
        "US-Northwest": "T0652",
        "US-Southwest": "T0655",
        "US-South": "T0240",
        "US-Alaska": "T0206",
        "MX": "T0570",
        "GE": "T0910",
        "QA": "T0597",
        "AE": "T0602",
        "KW": "T0617",
        "PK": "T0299",
        "TH": "T0321",
        "LK": "T0316",
        "ZA": "T0264",
        "KE": "T0254",
        "TZ": "T0263",
        "NG": "T0251",
        "AO": "T0245",
        "ET": "T0667",
        "MZ": "T0257",
        "ER": "T0691",
        "MR": "T1091",
        "SN": "T1092",
        "UY": "T0989",
        "CO": "T0679",
        "PE": "T0576",
        "IN-Northwest": "T0427",
        "IN-South": "T0541",
        "IN-East": "T0431",
        "IS": "T0733",
        "IL": "T0589",
        "SD": "T0691",
        "SY": "T0963",
        "YE": "T0605",
        "JO": "T0591",
        "IR": "T0586",
        "IQ": "T0592",
        "SO": "T0667",
        "MG": "T0789",
        "DJ": "T0667",
        "BO": "T0990",
        "VE": "T0626",
        "VN": "T0325",
        "MM": "T1065",
    }

    sea_distances = []
    for (region, origin), destination in product(origins.items(), destinations):
        o = lng.query("TerminalID == @origin").iloc[0]
        d = lng.query("TerminalID == @destination").iloc[0]

        route = searoute(
            origin=(o.Longitude, o.Latitude),
            destination=(d.Longitude, d.Latitude),
            # restrictions=["northwest", "suez", "panama"],
        )
        distance = int(route["properties"]["length"])

        o_location = o.Location if o.Location else o["State/Province"]
        d_location = d.Location if d.Location else d["State/Province"]

        sea_distances.append(
            pd.Series(
                {
                    "region_a": region,
                    "region_b": d.TerminalID,
                    "type": "sea route",
                    "value": distance,
                    "unit": "km",
                    "source": f"Distance between {o.TerminalName}, {o_location}, {o.Country} and {d.TerminalName}, {d_location}, {d.Country}",
                }
            )
        )
    sea_distances = pd.concat(sea_distances, axis=1).T

    destinations = {
        "EUSW": (36.39, -5.64, "Spain"),
        "EUS": (37.18, 14.53, "Italy"),
        "EUSE": (41.09, 25.46, "Greece"),
        "EUE": (48.54, 22.14, "Slovakia"),
        "EUNW": (58.28, -4.43, "Scotland"),
    }

    origins = {
        "CN-West": (37.44, 95.45),
        "CN-Southeast": (36.39, 108.62),
        "CN-Northeast": (38.05, 112.11),
        "EG": (26.61, 30.24),
        "DZ": (32.83, 3.75),
        "KZ": (48.74, 60.84),
        "LY": (27.20, 17.91),
        "MA": (28.69, -8, 82),
        "NA": (-21.91, 18.16),
        "SA": (24.14, 44.60),
        "TN": (34.09, 9.66),
        "TR": (38.93, 35.57),
        "UA": (48.66, 31.27),
        "OM": (21.14, 56.58),
        "GE": (41.92, 43.74),
        "AZ": (40.57, 47.74),
        "QA": (25.20, 51.20),
        "AE": (24.25, 54.69),
        "KW": (29.51, 47.49),
        "TM": (39.79, 59.22),
        "MN": (46.53, 103.54),
        "UZ": (42.33, 63.58),
        "MR": (20.57, -10.71),
        "IL": (31.77, 34.98),
        "SD": (16.51, 30.27),
        "ER": (15.32, 38.62),
        "EH": (24.89, -13.69),
        "SY": (35.00, 38.54),
        "YE": (16.02, 47.67),
        "IR": (32.06, 54.29),
        "IQ": (33.05, 42.77),
        "AF": (33.86, 64.88),
        "JO": (30.99, 36.23),
    }

    crowfly_distances = []
    for (s, s_latlon), (d, (d_lat, d_lon, d_name)) in product(
        origins.items(), destinations.items()
    ):
        distance = int(haversine((s_latlon[0], s_latlon[1]), (d_lat, d_lon)))
        crowfly_distances.append(
            pd.Series(
                {
                    "region_a": s,
                    "region_b": d,
                    "type": "as-the-crow-flies",
                    "value": distance,
                    "unit": "km",
                    "source": f"Haversine between {s} ({s_latlon[0]}, {s_latlon[1]}; latlon) and existing gas pipeline entrypoint in {d_name}  ({d_lat}, {d_lon}; latlon)",
                }
            )
        )
        crowfly_distances.append(
            pd.Series(
                {
                    "region_a": s,
                    "region_b": d,
                    "type": "as-the-hake-swims",
                    "value": 0,
                    "unit": "km",
                    "source": "All distances are assumed to be onland. Set underwater fraction to 0 accordingly.",
                }
            )
        )
    crowfly_distances = pd.concat(crowfly_distances, axis=1).T

    pd.concat([sea_distances, crowfly_distances]).to_csv(
        snakemake.output[0], index=False
    )
