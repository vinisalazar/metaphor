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

    cog = [
        "categories",
        "codes",
        "taxs",
        "pathways",
    ]
    ranks = "species genus family order class phylum kingdom domain".split()

    kinds = cog + ranks
    for kind in kinds:
        try:
            files = getattr(args, kind)
        except AttributeError:
            logging.error(
                f"Attribute {kind} wasn't found, please check the passed args:\n{args}"
            )
            continue
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
            k: v.rename(
                columns={"absolute": f"{k}_absolute", "relative": f"{k}_relative"}
            )
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


def parse_snakemake_args(snakemake):
    args = argparse.Namespace()
    args_dict = vars(args)

    for directive in "input", "output", "params":
        try:
            # The following block prevents conflict with Python's reserved keyword 'class'
            # An old problem for Pythonistas that work with taxonomy :)
            # See the rule 'lineage_parser' output directive to see class is spelled with a 'k'
            for k, v in getattr(snakemake, directive).items():
                if "klass" in k:
                    k = k.replace("klass", "class")
                args_dict[k] = v
        except AttributeError:
            pass

    return args


def parse_args():
    # Unfortunately this ugly block of code is required due to standardization of argument parsing across the workflow
    # 'Simple is better than complex.'
    # 'Special cases aren't special enough to break the rules.'
    parser = argparse.ArgumentParser()
    parser.add_argument("--categories")
    parser.add_argument("--codes")
    parser.add_argument("--taxs")
    parser.add_argument("--pathways")
    parser.add_argument("--species")
    parser.add_argument("--genus")
    parser.add_argument("--family")
    parser.add_argument("--order")
    parser.add_argument("--class")
    parser.add_argument("--phylum")
    parser.add_argument("--kingdom")
    parser.add_argument("--domain")
    parser.add_argument("--categories_absolute")
    parser.add_argument("--categories_relative")
    parser.add_argument("--codes_absolute")
    parser.add_argument("--codes_relative")
    parser.add_argument("--taxs_absolute")
    parser.add_argument("--taxs_relative")
    parser.add_argument("--pathways_absolute")
    parser.add_argument("--pathways_relative")
    parser.add_argument("--species_absolute")
    parser.add_argument("--species_relative")
    parser.add_argument("--genus_absolute")
    parser.add_argument("--genus_relative")
    parser.add_argument("--family_absolute")
    parser.add_argument("--family_relative")
    parser.add_argument("--order_absolute")
    parser.add_argument("--order_relative")
    parser.add_argument("--class_absolute")
    parser.add_argument("--class_relative")
    parser.add_argument("--phylum_absolute")
    parser.add_argument("--phylum_relative")
    parser.add_argument("--kingdom_absolute")
    parser.add_argument("--kingdom_relative")
    parser.add_argument("--domain_absolute")
    parser.add_argument("--domain_relative")
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
        parse_args_fn = parse_snakemake_args
    driver(main, snakemake, __file__, parse_args_fn)
