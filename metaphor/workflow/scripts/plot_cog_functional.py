#!/usr/bin/env python
import argparse
import logging


import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def create_heatmap(args):
    logging.info(f"Processing COG categories file: '{args.categories_file}'.")
    dataframe = pd.read_csv(args.categories_file, sep="\t", index_col=0)
    logging.info(f"{len(dataframe)} categories detected.")

    if args.filter_categories:
        logging.info("Filtering categories.")
        dataframe = dataframe.drop("Function unknown", errors="ignore")
        dataframe = dataframe.drop("General function prediction only", errors="ignore")
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

    # Sort rows by abundance
    reidx = dataframe.mean(axis=1).sort_values(ascending=False).index
    dataframe = dataframe.reindex(reidx)
    fig, ax = plt.subplots(figsize=(1 + int(len(dataframe.columns) * 0.5), 6))
    _ = sns.heatmap(dataframe, cmap="viridis", vmax=vmax, vmin=vmin, ax=ax, )
    # _ = ax.set_ylabel("COG categories", labelpad=25, rotation=0)
    outfile = args.categories_plot
    transparent = not getattr(args, "white_background", False)
    plt.savefig(outfile, dpi=args.dpi, bbox_inches="tight", transparent=transparent)
    logging.info(f"Generated plot: '{outfile}'.")


def main(args):
    create_heatmap(args)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--categories-file")
    parser.add_argument("--categories-plot")
    parser.add_argument("--filter-categories", action="store_true")
    parser.add_argument("--categories-cutoff", type=float)
    parser.add_argument("--white-background", action="store_true", default=False)
    parser.add_argument("--dpi", type=int, default=600)
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
