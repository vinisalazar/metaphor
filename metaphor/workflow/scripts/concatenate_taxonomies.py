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
    files = args.files

    # If it's a single file, make it into a list
    if isinstance(files, str):
        files = files.split()
    logging.info(f"Loading '{rank}' files:\n")
    logging.info("\n".join(files) + "\n")

    # Create dictionary of sample_names: dataframe
    df = {
        str(Path(file).stem).replace(f"_{rank}", ""): pd.read_csv(
            file, sep="\t", index_col=0
        )
        for file in files
    }

    # Rename dataframe columns
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
        # Use integers for taxids
        if outdf.index.dtype == "float64":
            # Ignore in case index is all NaN
            if not (all(outdf.index.isnull())):
                outdf.index = outdf.index.astype(int)
        outdf = outdf.round(6).T.fillna(0.0)
        outdf.index.name = "samples"
        try:
            attr = f"{count}_counts"
            outfile = getattr(args, attr)
            outdf.to_csv(outfile, sep="\t")
            if outfile:
                logging.info(f"Wrote {outdf.shape[1]} records to '{outfile}'.")
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

    parser.add_argument("--files")
    parser.add_argument("--absolute-counts")
    parser.add_argument("--relative-counts")
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
