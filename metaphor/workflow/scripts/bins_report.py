#!/usr/bin/env python
import argparse
import logging
from inspect import currentframe

import numpy as np
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

    # Add total values to legends
    for column in ("Quality threshold", "bin_set"):
        for ix, value in df[column].value_counts().items():
            df[column] = df[column].str.replace(ix, f"{ix} ({value})")

    df["SCG_set"] = df["SCG_set"].str.capitalize()
    df = df.rename(columns=rename_dict)
    qc_pass = "Quality threshold"
    domain = "Domain" if len(df["Domain"].value_counts()) > 1 else None

    return df, rename_dict, qc_pass, domain


def bin_quality(
    df,
    rename_dict,
    score_threshold,
    qc_pass,
    domain,
    binning_group,
    save=True,
    transparent=True,
    dpi=600,
    output_format="png",
):
    fig, ax = plt.subplots(figsize=(6, 6))

    plot_df = df

    def jitter(values, j):
        return values + np.random.normal(j, 0.2, values.shape)

    plot_df["Completeness (%)"] = jitter(plot_df["Completeness (%)"], 0)
    plot_df["Redundancy (%)"] = jitter(plot_df["Redundancy (%)"], 0)

    sns.scatterplot(
        x="Completeness (%)",
        y="Redundancy (%)",
        data=plot_df,  # [plot_df[qc_pass] == "Pass"],
        hue="Binning software",
        hue_order=df["Binning software"].value_counts().sort_index().index.to_list(),
        style=qc_pass,
        ax=ax,
        alpha=0.75,
        s=20,
    )

    _ = ax.set_ylim(-5, 105)
    _ = ax.set_xlim(-5, 105)
    _ = ax.set_ylabel("Redundancy (%)", labelpad=0, rotation=90)
    _ = ax.set_title(f"Bin quality scatterplot: {binning_group}")

    if save:
        func_name = currentframe().f_code.co_name
        outfile = f"output/binning/plots/{binning_group}/{func_name}.{output_format}"
        plt.savefig(outfile, dpi=dpi, bbox_inches="tight", transparent=transparent)
        logging.info(f"Generated plot: '{outfile}'.")

    return ax


def bin_scores(
    df,
    rename_dict,
    score_threshold,
    qc_pass,
    domain,
    binning_group,
    save=True,
    transparent=True,
    dpi=600,
    output_format="png",
):
    fig, ax = plt.subplots(figsize=(8, 8))
    hp = sns.histplot(
        x=rename_dict["bin_score"],
        data=df,
        hue="Binning software",
        hue_order=df["Binning software"].value_counts().sort_index().index.to_list(),
        multiple="stack",
        bins=np.linspace(-0.5, 1, 31),
    )
    _ = ax.set_xlim(-0.5, 1)
    _ = ax.set_ylabel("Number of bins", labelpad=0, rotation=90)

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
        outfile = f"output/binning/plots/{binning_group}/{func_name}.{output_format}"
        plt.savefig(outfile, dpi=dpi, bbox_inches="tight", transparent=transparent)
        logging.info(f"Generated plot: '{outfile}'.")

    return ax


def bin_quantity(
    df,
    rename_dict,
    score_threshold,
    qc_pass,
    domain,
    binning_group,
    save=True,
    transparent=True,
    dpi=600,
    output_format="png",
):
    fig, ax = plt.subplots(figsize=(6, 6))
    plot_data = pd.crosstab(df["Binning software"], df[qc_pass])

    plot_data[sorted(plot_data.columns)].plot(
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
    _ = ax.set_title(f"Number of bins: {binning_group}")
    _ = ax.set_ylabel("")
    _ = ax.set_xlabel(rename_dict["bin_set"])
    plt.xticks(rotation=0)

    if save:
        func_name = currentframe().f_code.co_name
        outfile = f"output/binning/plots/{binning_group}/{func_name}.{output_format}"
        plt.savefig(outfile, dpi=dpi, bbox_inches="tight", transparent=transparent)
        logging.info(f"Generated plot: '{outfile}'.")

    return ax


def bin_sizes(
    df,
    rename_dict,
    score_threshold,
    qc_pass,
    domain,
    binning_group,
    save=True,
    transparent=True,
    dpi=600,
    output_format="png",
):
    fig, ax = plt.subplots(figsize=(6, 6))
    bp = sns.boxplot(
        x=rename_dict["bin_set"],
        y=rename_dict["size"],
        order=df["Binning software"].value_counts().sort_index().index.to_list(),
        data=df,
        ax=ax,
        showfliers=False,
    )
    sp = sns.stripplot(
        x=rename_dict["bin_set"],
        y=rename_dict["size"],
        order=df["Binning software"].value_counts().sort_index().index.to_list(),
        data=df,
        ax=ax,
        jitter=True,
        edgecolor="k",
        linewidth=0.5,
        alpha=0.5,
    )

    # Axes formatting
    _ = ax.set_title(f"Size of bins: {binning_group}")
    _ = ax.set_ylabel("No. basepairs", rotation=90, labelpad=0)
    _ = ax.set_yscale("log")

    for patch in bp.patches:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, 0.5))

    if save:
        func_name = currentframe().f_code.co_name
        outfile = f"output/binning/plots/{binning_group}/{func_name}.{output_format}"
        plt.savefig(outfile, dpi=dpi, bbox_inches="tight", transparent=transparent)
        logging.info(f"Generated plot: '{outfile}'.")

    return ax


def bin_N50(
    df,
    rename_dict,
    score_threshold,
    qc_pass,
    domain,
    binning_group,
    save=True,
    transparent=True,
    dpi=600,
    output_format="png",
):
    fig, ax = plt.subplots(figsize=(6, 6))
    bp = sns.boxplot(
        x=rename_dict["bin_set"],
        y=rename_dict["n50"],
        data=df,
        ax=ax,
        showfliers=False,
        order=df["Binning software"].value_counts().sort_index().index.to_list(),
    )
    sp = sns.stripplot(
        x=rename_dict["bin_set"],
        y=rename_dict["n50"],
        order=df["Binning software"].value_counts().sort_index().index.to_list(),
        data=df,
        ax=ax,
        jitter=True,
        edgecolor="k",
        linewidth=0.5,
        alpha=0.5,
    )
    _ = ax.set_title(f"N50: {binning_group}")
    _ = ax.set_ylabel("")
    _ = ax.set_yscale("log")
    for patch in bp.patches:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, 0.3))

    if save:
        func_name = currentframe().f_code.co_name
        outfile = f"output/binning/plots/{binning_group}/{func_name}.{output_format}"
        plt.savefig(outfile, dpi=dpi, bbox_inches="tight", transparent=transparent)
        logging.info(f"Generated plot: '{outfile}'.")

    return ax


def main(args):
    bins_eval = args.bins_eval
    score_threshold = args.score_threshold
    binning_group = args.binning_group
    transparent = not getattr(args, "white_background", False)
    save = not getattr(args, "skip_save", False)
    dpi = getattr(args, "dpi", 600)
    output_format = getattr(args, "output_format", "png")
    df, rename_dict, qc_pass, domain = create_df(bins_eval, score_threshold)
    bin_quality(
        df,
        rename_dict,
        score_threshold,
        qc_pass,
        domain,
        binning_group,
        save,
        transparent=transparent,
        dpi=dpi,
        output_format=output_format,
    )
    bin_scores(
        df,
        rename_dict,
        score_threshold,
        qc_pass,
        domain,
        binning_group,
        save,
        transparent=transparent,
        dpi=dpi,
        output_format=output_format,
    )
    bin_quantity(
        df,
        rename_dict,
        score_threshold,
        qc_pass,
        domain,
        binning_group,
        save,
        transparent=transparent,
        dpi=dpi,
        output_format=output_format,
    )
    bin_sizes(
        df,
        rename_dict,
        score_threshold,
        qc_pass,
        domain,
        binning_group,
        save,
        transparent=transparent,
        dpi=dpi,
        output_format=output_format,
    )
    bin_N50(
        df,
        rename_dict,
        score_threshold,
        qc_pass,
        domain,
        binning_group,
        save,
        transparent=transparent,
        dpi=dpi,
        output_format=output_format,
    )
    print()


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--score-threshold", type=float)
    parser.add_argument("--bins-eval", type=str)
    parser.add_argument("--binning-group", type=str)
    parser.add_argument("--skip-save", action="store_true", default=False)
    parser.add_argument("--white-background", action="store_true", default=False)
    parser.add_argument("--output-format", type=str, default="png")
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
