#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import csv
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from _helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":
    configure_logging(snakemake)

    if len(snakemake.input["results"]) != len(set(snakemake.input["results"])):
        ds = pd.Series(list(snakemake.input["results"])).value_counts()
        raise ValueError(
            "List of input files contains duplicates. "
            "Script would later fail, so stopping here gracefully. "
            f"Affected files: {list(ds.loc[ds > 1].index)}"
        )

    dfs = []

    for fn in snakemake.input["results"]:
        fn = Path(fn)

        df = pd.read_csv(fn, keep_default_na=False, sep=";")

        dfs.append(df)

    df = pd.concat(dfs)
    df = df.set_index(
        [
            "scenario",
            "year",
            "wacc",
            "esc",
            "exporter",
            "importer",
            "category",
            "subcategory",
        ],
        verify_integrity=True,
    )

    df.to_csv(
        snakemake.output["results_csv"], sep=";", quotechar='"', quoting=csv.QUOTE_ALL
    )

    df.to_parquet(snakemake.output["results_parquet"], compression="snappy")
