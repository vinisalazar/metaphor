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

    # This handles repeated TaxIDs with different tax names
    # It only keeps the first tax name
    try:
        df = pd.concat(df.values(), axis=1)
    except pd.errors.InvalidIndexError:
        df = (
            v.groupby(v.index).agg({v.columns[0]: "first", v.columns[1]: sum})
            for v in df.values()
        )
        df = pd.concat(df, axis=1)

    for count in ("absolute", "relative"):
        outdf = df[[i for i in df.columns if count in i]]
        outdf = outdf.rename(lambda s: s.replace(f"_{count}", ""), axis="columns")
        outdf = outdf.round(6).reset_index()
        try:
            attr = f"{kind}_{count}"
            outfile = getattr(args, attr)
            outdf.to_csv(outfile, index=False, sep="\t")
            if outfile:
                logging.info(f"Wrote {len(outdf)} records to '{outfile}'.")
        except AttributeError:
            logging.error(
                f"Attribute {attr} wasn't found, please check the passed args:\n{args}"
            )
            pass


def parse_args():
    # Unfortunately this ugly block of code is required due to standardization of argument parsing across the workflow
    # 'Simple is better than complex.'
    # 'Special cases aren't special enough to break the rules.'
    parser = argparse.ArgumentParser()
    parser.add_argument("--functional-counts")
    parser.add_argument("--functional-relative-counts")
    parser.add_argument("--functional-absolute-counts")
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
