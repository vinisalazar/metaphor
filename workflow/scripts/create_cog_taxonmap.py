#!/usr/bin/env python
import argparse
import logging
from pathlib import Path

import pandas as pd

from utils import cog_csv_names


def main(args):
    cog, org, outfile = args.cog_csv, args.org_csv, args.taxonmap
    cog_df = pd.read_csv(
        cog,
        names=cog_csv_names,
        usecols=["NCBI Assembly ID", "Protein ID"],
    )

    org_names = [
        "NCBI Assembly ID",
        "Organism (genome) name",
        "taxid",
        "Taxonomic category used in COGs",
    ]
    org_df = pd.read_csv(org, names=org_names, usecols=["NCBI Assembly ID", "taxid"])
    df = cog_df.merge(org_df, on="NCBI Assembly ID")
    if outfile is None:
        outfile = cog.replace("cog.csv", "taxonmap.tsv")

    # Add required columns
    df["accession"] = df["Protein ID"].str.split(".", expand=True).iloc[:, 0]
    df["gi"] = ""

    # Rename and reorder
    df = df[["accession", "Protein ID", "taxid", "gi"]].rename(
        columns={"Protein ID": "accession.version"}
    )
    df = df.drop_duplicates(subset=["accession"])
    df.to_csv(outfile, index=False, sep="\t")
    logging.info(f"Generated Accession to TaxID map: '{outfile}'.")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cog-csv")
    parser.add_argument("--org-csv")
    parser.add_argument("--taxonmap")
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    # The driver function is standardized across scripts in this workflow
    # Please check the workflow/scripts/utils.py module for reference
    from utils import driver

    if "snakemake" not in locals():
        snakemake = None
    driver(main, snakemake, __file__, parse_args)
