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
    kind = args.kind
    files = args.functional_counts
    if isinstance(files, str):
        files = files.split()
    logging.info(f"Loading '{kind}' files:\n")
    logging.info("\n".join(files) + "\n")
    if kind == "codes":
        ix_col = [0, 1]
    else:
        ix_col = [
            0,
        ]
    df = {
        str(Path(file).stem).replace(f"_{kind}", ""): pd.read_csv(
            file, sep="\t", index_col=ix_col
        )
        for file in files
    }

    df = {
        k: v.rename(columns={"absolute": f"{k}_absolute", "relative": f"{k}_relative"})
        for k, v in df.items()
    }

    try:
        df = pd.concat(df.values(), axis=1).fillna(0)
    except pd.errors.InvalidIndexError:
        df = (
            v.groupby(v.index).agg({v.columns[0]: "first", v.columns[1]: sum})
            for v in df.values()
        )
        df = pd.concat(df, axis=1)

    # Transpose to have samples as rows
    df = df.T
    outfile = args.concatenated_functional_counts
    df.to_csv(outfile, sep="\t")
    logging.info(f"Wrote concatenated output to '{outfile}'.")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--functional-counts")
    parser.add_argument("--concatenated-functional-counts")
    parser.add_argument("--kind")
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
