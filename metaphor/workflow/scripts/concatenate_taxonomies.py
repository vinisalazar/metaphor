#!/usr/bin/env python
"""
Concatenates COG outputs (long format for each individual sample)
into a single table (wide with all samples).
"""
import argparse
import logging
from pathlib import Path
import pandas as pd


def main(args):
    rank = args.rank
    counts = args.taxonomy_counts
    outfile = args.concatenated_taxonomy_counts

    # If it's a single file, make it into a list
    if isinstance(counts, str):
        counts = counts.split()
    logging.info(f"Loading '{rank}' files:\n" + "\n".join(counts) + "\n")

    # Create dictionary of sample_names: dataframe
    df = {
        str(Path(file).stem).replace(f"_{rank}", ""): pd.read_csv(
            file, sep="\t", index_col=0
        )
        for file in counts
    }

    # This handles repeated TaxIDs with different tax names
    # It only keeps the first tax name
    try:
        df = pd.concat(df.values()).fillna(0)
    except pd.errors.InvalidIndexError:
        df = (
            v.groupby(v.index).agg({v.columns[0]: "first", v.columns[1]: sum})
            for v in df.values()
        )
        df = pd.concat(df).fillna(0)

    df.to_csv(outfile, sep="\t")
    logging.info(f"Wrote {len(df)} records to '{outfile}'.")


def parse_args():
    # Unfortunately this ugly block of code is required due to standardization of argument parsing across the workflow
    # 'Simple is better than complex.'
    # 'Special cases aren't special enough to break the rules.'
    parser = argparse.ArgumentParser()

    parser.add_argument("--taxonomy-counts")
    parser.add_argument("--concatenated-taxonomy-counts")
    parser.add_argument("--rank")
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
