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

    kinds = "categories", "codes"
    idx = (0, [0, 1])
    for kind, ix_col in zip(kinds, idx):
        logging.info(f"Loading '{kind}' files.")
        files = getattr(args, kind)
        if isinstance(files, str):
            files = files.split()
        df = {
            str(Path(file).stem).replace(f"_{kind}", ""): pd.read_csv(
                file, sep="\t", index_col=ix_col
            )
            for file in files
        }

        df = {
            k: v.rename(
                columns={"absolute": f"{k}_absolute", "relative": f"{k}_relative"}
            )
            for k, v in df.items()
        }

        df = pd.concat(df.values(), axis=1)

        for count in ("absolute", "relative"):
            outdf = df[[i for i in df.columns if count in i]]
            outdf = outdf.rename(lambda s: s.replace(f"_{count}", ""), axis="columns")
            outdf = outdf.round(6).reset_index()
            outfile = getattr(args, f"concat_{kind}_{count}")
            outdf.to_csv(outfile, index=False, sep="\t")
            if outfile:
                logging.info(f"Wrote {len(outdf)} records to '{outfile}'.")


def parse_snakemake_args(snakemake):
    args = argparse.Namespace()
    args_dict = vars(args)

    for directive in "input", "output", "params":
        try:
            for k, v in getattr(snakemake, directive).items():
                args_dict[k] = v
        except AttributeError:
            pass

    return args


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--categories")
    parser.add_argument("--codes")
    parser.add_argument("--concat_categories_absolute")
    parser.add_argument("--concat_categories_relative")
    parser.add_argument("--concat_codes_absolute")
    parser.add_argument("--concat_codes_relative")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%m/%d/%Y %H:%M:%S",
    )
    logging.info(f"Starting script '{__file__.split('/')[-1]}'.")
    logging.debug(f"Full script path: '{__file__}'.")
    if "snakemake" in locals():
        logging.basicConfig(filename=str(snakemake.log))
        args = parse_snakemake_args(snakemake)
    else:
        args = parse_args()
    try:
        main(args)
        logging.info("Done.")
    except Exception as e:
        logging.error(e)