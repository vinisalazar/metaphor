#!/usr/bin/env python
import argparse
import logging
from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


#######################################
# Category heatmap
#######################################


def create_heatmap(categories_file, filter_categories=True, cutoff=0.01):
    logging.info(f"Processing COG categories file: '{categories_file}'.")
    dataframe = pd.read_csv(categories_file, sep="\t", index_col=0)
    logging.info(f"{len(dataframe)} categories detected.")

    if filter_categories:
        logging.info("Filtering categories.")
        dataframe = dataframe.drop("Function unknown")
        dataframe = dataframe.drop("General function prediction only")
        dataframe = dataframe[dataframe > 0.01]
        dataframe = dataframe.dropna()
        filtered = (abs(dataframe.sum() - 1)).mean()
        logging.info(f"{len(dataframe)} categories left after filtering.")
        dataframe = dataframe / dataframe.sum()
        logging.info("Normalising data after filtering.")
        vmax, vmin = None, None
    else:
        nlargest = categories.sum(axis=1).nlargest().index.to_list()
        for ix in nlargest:
            if ix not in ("Function unknown", "General function prediction only"):
                vmax = dataframe.loc[ix].max()
                break

        vmin = cutoff

    outfile = Path(categories_file).with_suffix(".png")
    fig, ax = plt.subplots(figsize=(3 + len(dataframe.columns), 6))
    sns.heatmap(dataframe, cmap="viridis", vmax=vmax, vmin=vmin, ax=ax)
    plt.savefig(outfile, bbox_inches="tight")
    logging.info(f"Generated plot: '{outfile}'.")


#######################################
# Taxa bar plots
#######################################


def calculate_legend_width(index):
    max_len = max(len(str(i)) for i in index)
    if max_len > 20:
        return (3, 1)
    elif max_len < 10:
        return (8, 1)
    else:
        return (4, 1)


def format_index(index):
    new_index = []
    for s in index:
        try:
            s = str(s)
            if len(s.split()) >= 2 and "Undetermined/other" not in s:
                s = s.split()[0][0] + ". " + " ".join(s.split()[1:])
                s = "$" + s + "$"
        except (TypeError, ValueError):
            pass
        new_index.append(s)

    return pd.Index(new_index, name=index.name)


def create_tax_barplot(dataframe, save=True):
    figsize = (10, len(dataframe.columns))
    dataframe.index = format_index(dataframe.index)
    rank = dataframe.index.name
    logging.info(f"Generating plot for rank '{rank}'.")
    fig, axs = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=figsize,
        gridspec_kw={"width_ratios": calculate_legend_width(dataframe.index)},
    )
    dataframe.T.plot(kind="barh", stacked=True, ax=axs[0], cmap="tab20c")
    handles, labels = axs[0].get_legend_handles_labels()
    axs[0].get_legend().remove()
    axs[0].set_xlabel("Relative abundance")
    axs[1].legend(handles, labels, title=rank.capitalize())
    axs[1].set_xticks([])
    axs[1].set_yticks([])
    axs[1].axis("off")
    if save:
        outfile = f"output/annotation/cog/COG_{rank}_relative.png"
        plt.savefig(outfile, bbox_inches="tight")
        logging.info(f"Generated plot: '{outfile}'.")


def process_rank_files(args):
    ranks = "species genus family order class phylum domain".split()
    rank_files = (
        getattr(args, rank) if rank != "class" else getattr(args, "klass")
        for rank in ranks
    )
    rank_df_list, rank_df = [], []  # rank_df is going to become a dataframe

    cutoff = 0.05
    for rank, file in zip(ranks, rank_files):
        sample = file.split("/")[-1].replace(f"_{rank}.tsv", "")
        series = pd.read_csv(file, sep="\t", index_col=0)["relative"]
        series.name = sample
        rank_df.append(series)
        rank_df = pd.concat(rank_df, axis=1)
        rank_df = rank_df[rank_df.sum(axis=1) > args.tax_cutoff]
        rank_df = rank_df.loc[[i for i in rank_df.index if not isinstance(i, float)]]
        rank_df.index.name = rank

        # Sort in descending abundance
        rank_df = rank_df.loc[rank_df.sum(axis=1).sort_values(ascending=False).index]

        # Group low abundance and undetermined taxa
        filtered = abs(rank_df.sum() - 1)
        if filtered.sum() < args.tax_cutoff:
            pass
        else:
            rank_df.loc["Undetermined/other"] = filtered

        rank_df = rank_df[sorted(rank_df.columns, reverse=True)]

        rank_df_list.append(rank_df)

    return rank_df_list


def main(args):
    categories_file, filter_categories, categories_cutoff = (
        args.categories_file,
        args.filter_categories,
        args.categories_cutoff,
    )
    create_heatmap(categories_file, filter_categories, categories_cutoff)
    rank_df_list = process_rank_files(args)
    for dataframe in rank_df_list:
        create_tax_barplot(dataframe, save=True)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--categories_file")
    parser.add_argument("--species")
    parser.add_argument("--genus")
    parser.add_argument("--family")
    parser.add_argument("--order")
    parser.add_argument("--klass")
    parser.add_argument("--phylum")
    parser.add_argument("--kingdom")
    parser.add_argument("--domain")
    parser.add_argument("--categories_plt")
    parser.add_argument("--filter_categories")
    parser.add_argument("--categories_cutoff")
    parser.add_argument("--tax_cutoff")
    args = parser.parse_args()


if __name__ == "__main__":
    # The driver function is standardized across scripts in this workflow
    # Please check the workflow/scripts/utils.py module for reference
    from utils import driver

    if "snakemake" not in locals():
        snakemake = None
    driver(main, snakemake, __file__)
