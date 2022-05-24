#!/usr/bin/env python
import argparse
import logging
from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def create_heatmap(args):
    logging.info(f"Processing COG categories file: '{args.categories_file}'.")
    dataframe = pd.read_csv(args.categories_file, sep="\t", index_col=0).T
    logging.info(f"{len(dataframe)} categories detected.")

    if args.filter_categories:
        logging.info("Filtering categories.")
        dataframe = dataframe.drop("Function unknown")
        dataframe = dataframe.drop("General function prediction only")
        dataframe = dataframe[dataframe > args.categories_cutoff]
        dataframe = dataframe.dropna()
        filtered = (abs(dataframe.sum() - 1)).mean()
        logging.info(f"{len(dataframe)} categories left after filtering.")
        dataframe = dataframe / dataframe.sum()
        logging.info("Normalising data after filtering.")
        vmax, vmin = None, None
    else:
        nlargest = dataframe.sum(axis=1).nlargest().index.to_list()
        for ix in nlargest:
            if ix not in ("Function unknown", "General function prediction only"):
                vmax = dataframe.loc[ix].max()
                break

        vmin = args.categories_cutoff

    # TODO: improve this path construction
    outfile = str(Path(args.categories_file).with_suffix(".png")).replace(
        "tables", "plots"
    )
    fig, ax = plt.subplots(figsize=(3 + len(dataframe.columns), 6))
    sns.heatmap(dataframe, cmap="viridis", vmax=vmax, vmin=vmin, ax=ax)
    plt.savefig(outfile, bbox_inches="tight")
    logging.info(f"Generated plot: '{outfile}'.")


def main(args):
    create_heatmap(args)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--categories_file")
    parser.add_argument("--categories_plot")
    parser.add_argument("--filter_categories", action="store_true")
    parser.add_argument("--categories_cutoff", type=float)
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
