#!/usr/bin/env python
import argparse
from audioop import avg
import logging
from pathlib import Path
from textwrap import wrap

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

#######################################
# Taxa bar plots
#######################################

# String used to describe the low abundance taxa that are filtered/bundled together.
filtered_label = "Undetermined/Low abundance"
plt.style.use("seaborn")


def calculate_legend_width(index):
    max_len = max(len(str(i)) for i in index)
    if max_len > 20:
        return (3, 1)
    elif max_len < 10:
        return (8, 1)
    else:
        return (4, 1)


def format_index(index, abbreviate_genus=True):
    new_index = []
    for s in index:
        try:
            s = str(s)
            if len(s.split()) >= 2 and (s != filtered_label):
                if abbreviate_genus:
                    s = s.split()[0][0] + ". " + " ".join(s.split()[1:])
                s = "$" + s + "$"
        except (TypeError, ValueError):
            pass
        new_index.append(s)

    return pd.Index(new_index, name=index.name)


def create_tax_barplot(
    dataframe, tax_cutoff, colormap, save=True, outfile=None, abbreviate_genus=True
):
    figsize = (10, len(dataframe.columns))
    if dataframe.index.name == "species":
        dataframe.index = format_index(dataframe.index, abbreviate_genus)
    rank = dataframe.index.name
    logging.info(f"Generating plot for rank '{rank}'.")
    fig, axs = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=figsize,
        gridspec_kw={"width_ratios": calculate_legend_width(dataframe.index)},
    )
    if (undetermined_label := "Unnamed: 1") in dataframe.index:
        dataframe.loc[filtered_label] = (
            dataframe.loc[filtered_label] + dataframe.loc[undetermined_label]
        )
        dataframe = dataframe.drop(undetermined_label)
    dataframe.T.plot(kind="barh", stacked=True, ax=axs[0], cmap=colormap)
    handles, labels = axs[0].get_legend_handles_labels()
    axs[0].get_legend().remove()
    axs[0].set_xlabel("Relative abundance")
    axs[0].set_ylabel("")

    # Bundle unnamed taxa with other filtered taxa
    if filtered_label in dataframe.index:
        # Get legend for filtered average
        avg_filtered = dataframe.loc[filtered_label].mean()
        avg_filtered = round(avg_filtered * 100, 2)
        if avg_filtered > 1e-6:
            if tax_cutoff:
                handle, label = [
                    Rectangle(
                        (0, 0), 1, 1, fc="w", fill=False, edgecolor="none", linewidth=0
                    )
                ], [f"Showing the {tax_cutoff} most abundant taxa."]
                handles += handle
                labels += label
            handle, label = [
                Rectangle(
                    (0, 0), 1, 1, fc="w", fill=False, edgecolor="none", linewidth=0
                )
            ], [f"Avg. of {avg_filtered}% reads were filtered as low abundance."]
            handles += handle
            labels += label
            labels = ["\n".join(wrap(l, 20)) for l in labels]
        else:
            dataframe = dataframe.drop(filtered_label)

    axs[1].legend(handles, labels, title=rank.capitalize() if rank else "")
    axs[1].set_xticks([])
    axs[1].set_yticks([])
    axs[1].axis("off")
    if save:
        if not outfile:
            outfile = f"output/annotation/cog/plots/COG_{rank}_relative.png"
        plt.savefig(outfile, dpi=600, bbox_inches="tight", transparent=True)
        logging.info(f"Generated plot: '{outfile}'.")


def process_rank_file(args):
    cutoff = int(args.tax_cutoff)
    rank_df = pd.read_csv(args.taxonomy_relative_counts, sep="\t", index_col=0)
    rank_df = rank_df.loc[[i for i in rank_df.index if not isinstance(i, float)]]

    # Sort in descending abundance
    rank_df = rank_df.loc[rank_df.sum(axis=1).sort_values(ascending=False).index]

    # Apply cutoff
    if cutoff:
        rank_df = rank_df.iloc[:cutoff]
        # Group low abundance and undetermined taxa
        filtered = abs(rank_df.sum() - 1)
        rank_df.loc[filtered_label] = filtered

        rank_df = rank_df[sorted(rank_df.columns, reverse=True)]

    return rank_df


def main(args):
    rank_df = process_rank_file(args)
    create_tax_barplot(
        rank_df,
        args.tax_cutoff,
        args.colormap,
        save=True,
        outfile=args.taxonomy_barplot,
    )


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--taxonomy-relative-counts")
    parser.add_argument("--tax-cutoff")
    parser.add_argument("--taxonomy-barplot")
    parser.add_argument("--rank")
    parser.add_argument("--colormap", default="tab20c")
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
