"""
Parses tabular output file with COG annotations.

Because the output of Diamond can be fairly large, and the COG database files are also not small,
the current version of this script only uses the required columns from both files. This greatly
reduces memory usage, as lots of unnecessary columns are required.

In the future, a configurable option may be added to enable the output of all columns.
"""

import sys
import logging
import argparse
import subprocess
from pathlib import Path
from functools import lru_cache
from utils import cog_csv_names

import pandas as pd


def main(args):

    (dmnd_out, tax_out,) = (
        args.dmnd_out,
        args.tax_out,
    )

    # Load data
    logging.info(f"Loading annotation data: '{dmnd_out}'.")
    df = pd.read_csv(
        dmnd_out,
        sep="\t",
        usecols=["qseqid", "sseqid", "bitscore", "staxids", "sscinames"],
    )

    # Keep the best score
    df = df.sort_values("bitscore", ascending=False)
    df = df.drop_duplicates("qseqid")

    logging.info(f"Loaded {len(df)} records.")

    df = df[["taxid", "taxname"]].value_counts()
    df.name = "absolute"
    df.to_csv(tax_out, sep="\t")
    logging.info(f"Wrote {len(df)} rows to '{tax_out}'.")
    del df


def parse_args():
    # TODO: improve these arguments and add parser description
    parser = argparse.ArgumentParser()
    parser.add_argument("--dmnd_out")
    parser.add_argument("--tax_out")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    # The driver function is standardized across scripts in this workflow
    # Please check the workflow/scripts/utils.py module for reference
    from utils import driver

    if "snakemake" not in locals():
        snakemake = None
        parse_args_fn = parse_args
    else:
        parse_args_fn = None
    driver(main, snakemake, __file__, parse_args_fn)
