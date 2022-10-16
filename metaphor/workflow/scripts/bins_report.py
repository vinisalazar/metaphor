#!/usr/bin/env python
import argparse
import logging
from inspect import currentframe

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.colors import ListedColormap

sns.set_palette("colorblind")


def create_df(file, score_threshold):
    df = pd.read_csv(file, sep="\t")

    rename_dict = {
        "SCG_completeness": "Completeness (%)",
        "SCG_redundancy": "Redundancy (%)",
        "size": "Bin size",
        "bin_set": "Binning software",
        "SCG_set": "Domain",
        "n50": "N50",
        "bin_score": "Bin score",
    }

    df["Quality threshold"] = df["bin_score"].apply(
        lambda x: "Pass" if x >= score_threshold else "Fail"
    )
    df["SCG_set"] = df["SCG_set"].str.capitalize()
    df = df.rename(columns=rename_dict)
    qc_pass = "Quality threshold"
    domain = "Domain" if len(df["Domain"].value_counts()) > 1 else None

    return df, rename_dict, qc_pass, domain


def bin_quality(
    df, rename_dict, score_threshold, qc_pass, domain, binning_group, save=True
):
    fig, ax = plt.subplots()
    sns.scatterplot(
        x="Completeness (%)",
        y="Redundancy (%)",
        data=df,
        hue="Binning software",
        style=qc_pass,
        ax=ax,
    )

    ax.set_ylim(-5, 105)
    ax.set_xlim(-5, 105)
    _ = ax.set_title("Bin quality scatterplot")

    if save:
        func_name = currentframe().f_code.co_name
        outfile = f"output/binning/plots/{binning_group}/{func_name}.png"
        plt.savefig(outfile, bbox_inches="tight")
        logging.info(f"Generated plot: '{outfile}'.")

    return ax


def bin_scores(
    df, rename_dict, score_threshold, qc_pass, domain, binning_group, save=True
):
    fig, ax = plt.subplots()
    hp = sns.histplot(
        x=rename_dict["bin_score"], data=df, hue="Binning software", multiple="stack"
    )
    ax.set_xlim(-2.5, 1)

    if score_threshold:
        plt.axvline(
            score_threshold, 0.0, 1, alpha=1, zorder=2, color="k", linestyle="--"
        )
        trans = ax.get_xaxis_transform()
        plt.text(
            score_threshold * 1.01,
            0.60,
            qc_pass + f"({score_threshold})",
            transform=trans,
        )

    if save:
        func_name = currentframe().f_code.co_name
        outfile = f"output/binning/plots/{binning_group}/{func_name}.png"
        plt.savefig(outfile, bbox_inches="tight")
        logging.info(f"Generated plot: '{outfile}'.")

    return ax


def bin_quantity(
    df, rename_dict, score_threshold, qc_pass, domain, binning_group, save=True
):
    fig, ax = plt.subplots()
    plot_data = pd.crosstab(df["Binning software"], df[qc_pass])

    plot_data.plot(
        kind="bar",
        stacked=True,
        ax=ax,
        cmap=ListedColormap(
            ["#4196C5", "#E6AB44"]
        ),  # Change this if switching color palette
        edgecolor="k",
    )
    # From https://stackoverflow.com/questions/41296313/stacked-bar-chart-with-centered-labels
    for c in ax.containers:
        # Optional: if the segment is small or 0, customize the labels
        labels = [v.get_height() if v.get_height() > 0 else "" for v in c]
        # Remove the labels parameter if it's not needed for customized labels
        ax.bar_label(c, labels=labels, label_type="center")
    _ = ax.set_title("Number of bins")
    _ = ax.set_ylabel("")
    _ = ax.set_xlabel(rename_dict["bin_set"])
    plt.xticks(rotation=0)

    if save:
        func_name = currentframe().f_code.co_name
        outfile = f"output/binning/plots/{binning_group}/{func_name}.png"
        plt.savefig(outfile, bbox_inches="tight")
        logging.info(f"Generated plot: '{outfile}'.")

    return ax


def bin_sizes(
    df, rename_dict, score_threshold, qc_pass, domain, binning_group, save=True
):
    fig, ax = plt.subplots()
    bp = sns.boxplot(
        x=rename_dict["bin_set"],
        y=rename_dict["size"],
        data=df,
        ax=ax,
        showfliers=False,
    )
    sp = sns.stripplot(x=rename_dict["bin_set"], y=rename_dict["size"], data=df, ax=ax)
    _ = ax.set_title("Size of bins (# nucleotides)")
    _ = ax.set_ylabel("")

    for patch in bp.patches:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, 0.5))

    if save:
        func_name = currentframe().f_code.co_name
        outfile = f"output/binning/plots/{binning_group}/{func_name}.png"
        plt.savefig(outfile, bbox_inches="tight")
        logging.info(f"Generated plot: '{outfile}'.")

    return ax


def bin_N50(
    df, rename_dict, score_threshold, qc_pass, domain, binning_group, save=True
):
    fig, ax = plt.subplots()
    bp = sns.boxplot(
        x=rename_dict["bin_set"], y=rename_dict["n50"], data=df, ax=ax, showfliers=False
    )
    sp = sns.stripplot(x=rename_dict["bin_set"], y=rename_dict["n50"], data=df, ax=ax)
    _ = ax.set_title("N50")
    _ = ax.set_ylabel("")
    for patch in bp.patches:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, 0.3))

    if save:
        func_name = currentframe().f_code.co_name
        outfile = f"output/binning/plots/{binning_group}/{func_name}.png"
        plt.savefig(outfile, bbox_inches="tight")
        logging.info(f"Generated plot: '{outfile}'.")

    return ax


def main(args):
    bins_eval = args.bins_eval
    score_threshold = args.score_threshold
    binning_group = args.binning_group
    save = not getattr(args, "skip_save", False)
    df, rename_dict, qc_pass, domain = create_df(bins_eval, score_threshold)
    bin_quality(df, rename_dict, score_threshold, qc_pass, domain, binning_group, save)
    bin_scores(df, rename_dict, score_threshold, qc_pass, domain, binning_group, save)
    bin_quantity(df, rename_dict, score_threshold, qc_pass, domain, binning_group, save)
    bin_sizes(df, rename_dict, score_threshold, qc_pass, domain, binning_group, save)
    bin_N50(df, rename_dict, score_threshold, qc_pass, domain, binning_group, save)
    print()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--score-threshold", type=float)
    parser.add_argument("--bins-eval", type=str)
    parser.add_argument("--binning-group", type=str)
    parser.add_argument("--skip-save", action="store_true", default=False)
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
