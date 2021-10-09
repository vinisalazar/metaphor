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

    for rule_input in ("categories", "codes"):
        args_dict[rule_input] = snakemake.input[rule_input]

    for rule_output in (
        "concat_categories_absolute",
        "concat_categories_relative",
        "concat_codes_absolute",
        "concat_codes_relative",
    ):
        args_dict[rule_output] = snakemake.output[rule_output]

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


if "snakemake" in locals():
    logging.basicConfig(
        filename=str(snakemake.log),
        encoding="utf-8",
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%m/%d/%Y %H:%M:%S",
    )
    logging.info(f"Starting script {__file__.split('/')[-1]}.")
    logging.debug(f"Full script path: {__file__}")
    args = parse_snakemake_args(snakemake)
    main(args)
    logging.info("Done.")
elif __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%m/%d/%Y %H:%M:%S",
    )
    logging.info(f"Starting script '{__file__.split('/')[-1]}'.")
    logging.debug(f"Full script path: '{__file__}'.")
    args = parse_args()
    main(args)
    logging.info("Done.")
