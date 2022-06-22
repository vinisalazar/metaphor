"""
Parses tabular output file with COG annotations.

Because the output of Diamond can be fairly large, and the COG database files are also not small,
the current version of this script only uses the required columns from both files. This greatly
reduces memory usage, as lots of unnecessary columns are required.

In the future, a configurable option may be added to enable the output of all columns.
"""

import logging
import argparse

import pandas as pd

from utils import write_dfs


def main(args):

    (dmnd_out, cov_depths, tax_out_absolute, tax_out_relative) = (
        args.dmnd_out,
        args.coverage_depths,
        args.tax_out_absolute,
        args.tax_out_relative,
    )

    # Load data
    logging.info(f"Loading annotation data: '{dmnd_out}'.")
    df = pd.read_csv(
        dmnd_out,
        sep="\t",
        usecols=["qseqid", "sseqid", "bitscore", "staxids", "sscinames"],
        index_col="qseqid",
    )
    logging.info(f"Loaded {len(df)} records.")

    # Keep the best score, format taxids
    df = df.sort_values("bitscore", ascending=False)
    df = df[df.index.duplicated(keep="first")]
    df["staxids"] = df["staxids"].astype(int, errors="ignore")

    logging.info("Merging coverage depths.")
    covs = pd.read_csv(cov_depths, sep="\t", index_col=0)
    df["contig"] = df.index.map(lambda n: "_".join(n.split("_")[:-1]))
    df = df.merge(covs, left_on="contig", right_index=True)

    # Calculating counts from cov depths
    samples = [i for i in df.columns if i.endswith(".sorted.bam")]

    df = round(
        df[["staxids", "sscinames"] + samples].groupby(["staxids", "sscinames"]).sum(),
        4,
    )
    write_dfs(df, tax_out_absolute, tax_out_relative)


def parse_args():
    # TODO: improve these arguments and add parser description
    parser = argparse.ArgumentParser()
    parser.add_argument("--dmnd_out")
    parser.add_argument("--coverage_depths")
    parser.add_argument("--tax_out_relative")
    parser.add_argument("--tax_out_absolute")
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
