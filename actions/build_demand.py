import logging

from _helpers import configure_logging

import pandas as pd

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    configure_logging(snakemake)
        
    df = (
        pd.concat(
            [
                pd.read_csv(fn, sep=";")
                for fn in snakemake.input.synde
            ],
            ignore_index=True
        )
        .groupby("region_code")["Electricity demand"]
        .sum()
        .div(1e3)
    )

    regions = list(snakemake.config["regions"].keys())

    substitute_values = pd.Series({
        "US": 4650375,
        "MX": 1767142,
        "NA": 6279,
        "CA": 624010,
        "AU": 258647,
        "MG": 40573,
        "MR": 15197,
        "DJ": 312,
        "AF": 16748,
        "SO": 4897,
        "EH": 1925,
        "AE": 225800,
        "QA": 64661,
        "KW": 102427,
    })
    df = substitute_values.combine_first(df)
    subregion_shares = {
        "AR-South": 0.06,
        "AR-North": 0.94,
        "AU-West": 0.11,
        "AU-East": 0.89,
        "BR-Northeast": 0.28,
        "BR-Southeast": 0.62,
        "CA-East": 0.71,
        "CN-Northeast": 0.48,
        "CN-Southeast": 0.38,
        "CN-West": 0.14,
        "US-Northeast": 0.43,
        "US-Southeast": 0.19,
        "US-Northwest": 0.06,
        "US-Southwest": 0.18,
        "US-South": 0.13,
        "US-Alaska": 0.01,
        "IN-Northwest": 0.41,
        "IN-South": 0.31,
        "IN-East": 0.28,
    }

    for subregion, share in subregion_shares.items():
        df[subregion] = df[subregion.split("-")[0]] * share

    factor = snakemake.config["demand_factor"]
    df = df.reindex(regions).dropna().astype(int).sort_index() * factor # GWh

    df.index.name = "region"
    df.name = "demand [GWh]"

    df.to_csv(snakemake.output[0])