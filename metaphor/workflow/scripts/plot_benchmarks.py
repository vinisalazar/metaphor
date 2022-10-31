#!/usr/bin/env python
import argparse
import logging
from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


def load_data(filepath):
    df = pd.read_csv(filepath)
    df = df.sort_values("s", ascending=False)
    df["m"] = df["s"] / 60
    df["h"] = df["m"] / 60
    return df


def runtime_barplot_sum(df, **kwargs):
    """
    Generates barplot of sum of rules' runtimes.
    """
    fig, ax = plt.subplots(figsize=(4, 8))
    outfile, time_unit, cutoff = (
        kwargs["outfile"],
        kwargs["time_unit"],
        kwargs["time_cutoff"],
    )
    plot_data = (
        (df.groupby(by=["module", "rule"])[time_unit].sum())
        .reset_index()
        .sort_values(time_unit, ascending=False)
        .query(f"{time_unit} > {cutoff}")
    )
    sns.barplot(
        x=time_unit,
        y="rule",
        hue="module",
        dodge=False,
        data=plot_data,
        palette="Dark2",
    )
    plt.legend(loc="lower right", title="Module")

    xlabels = {"h": "(hours)", "m": "(minutes)", "s": "(seconds)"}

    plt.title("Total rule runtime")
    ax.set_xlabel("Runtime " + xlabels[time_unit], fontsize=12)
    ax.set_ylabel("Rule name", rotation=0, fontsize=12)
    if outfile:
        plt.savefig(outfile, dpi=600, bbox_inches="tight")


def runtime_barplot_errorbar(df, n_samples=None, **kwargs):
    """
    Generates barplot of rules' runtimes for each sample.

    The time-per-sample is corrected in rules that run once for all samples.
    """
    fig, ax = plt.subplots(figsize=(4, 8))
    outfile, time_unit, cutoff = (
        kwargs["outfile"],
        kwargs["time_unit"],
        kwargs["time_cutoff"],
    )
    df[time_unit + "_corrected"] = df.apply(
        lambda row: row["m"] / df["sample"].str[:6].value_counts().shape[0]
        if not isinstance(row["sample"], str)
        else row["m"],
        axis=1,
    )
    plot_data = df.sort_values(time_unit + "_corrected", ascending=False).query(
        f"{time_unit} > {cutoff}"
    )

    # Divide rules that run once for all samples by the number of samples
    if n_samples:
        plot_data[time_unit + "_corrected"] = plot_data.apply(
            lambda row: row[time_unit + "_corrected"] / float(n_samples)
            if row["sample"] != row["sample"]
            else row[time_unit + "_corrected"],
            axis=1,
        )

    sns.barplot(
        x=time_unit + "_corrected",
        y="rule",
        hue="module",
        dodge=False,
        data=plot_data,
        palette="Dark2",
    )
    plt.legend(loc="lower right", title="Module")

    xlabels = {"h": "(hours)", "m": "(minutes)", "s": "(seconds)"}

    plt.title("Mean runtime per sample")
    ax.set_xlabel("Runtime " + xlabels[time_unit], fontsize=12)
    ax.set_ylabel("Rule name", rotation=0, fontsize=12, position=(0, 0.7))
    if outfile:
        plt.savefig(outfile, dpi=600, bbox_inches="tight")


def memory_barplot_sum(df, **kwargs):

    """
    Generates barplot of sum of rules' memory usage.
    """
    outfile, memory_unit, cutoff, gb = (
        kwargs["outfile"],
        kwargs["memory_unit"],
        kwargs["memory_cutoff"],
        kwargs["gb"],
    )
    if df[memory_unit].all() == "-":
        logging.info("Skipping plot, no memory stats recorded.")
        return
    fig, ax = plt.subplots(figsize=(4, 8))
    df_ = df.copy()
    df_[memory_unit]
    if gb:
        df_[memory_unit] = df_[memory_unit] / 2**10
    plot_data = (
        (df_.groupby(by=["module", "rule"])[memory_unit].sum())
        .reset_index()
        .sort_values(memory_unit, ascending=False)
        .query(f"{memory_unit} > {cutoff}")
    )
    sns.barplot(
        x=memory_unit,
        y="rule",
        hue="module",
        dodge=False,
        data=plot_data,
        palette="Dark2",
    )
    plt.legend(loc="lower right", title="Module")

    xlabels = {
        "max_vms": "Virtual Memory Size (MB)",
        "max_uss": "Unique Set Size (MB)",
        "max_pss": "Proportional Set Size (MB)",
    }

    plt.title("Sum of rules' memory usage")
    xlabel = xlabels[memory_unit]
    if gb:
        xlabel = xlabel.replace("MB", "GB")
    ax.set_xlabel(xlabel, labelpad=15, fontsize=12)
    ax.set_ylabel("Rule name", rotation=0, fontsize=12)
    if outfile:
        plt.savefig(outfile, dpi=600, bbox_inches="tight")


def memory_barplot_errorbar(df, n_samples=None, **kwargs):
    """
    Generates barplot of sum of rules' memory usage.
    """
    outfile, memory_unit, cutoff, gb = (
        kwargs["outfile"],
        kwargs["memory_unit"],
        kwargs["memory_cutoff"],
        kwargs["gb"],
    )
    if df[memory_unit].all() == "-":
        logging.info("Skipping plot, no memory stats recorded.")
        return
    fig, ax = plt.subplots(figsize=(4, 8))
    df_ = df.copy()
    if gb:
        df_[memory_unit] = df_[memory_unit] / 2**10
    plot_data = df_.query(f"{memory_unit} > {cutoff}").sort_values(
        memory_unit, ascending=False
    )

    # Divide rules that run once for all samples by the number of samples
    if n_samples:
        plot_data[memory_unit] = plot_data.apply(
            lambda row: row[memory_unit] / n_samples
            if row["sample"] != row["sample"]
            else row[memory_unit],
            axis=1,
        )

    sns.barplot(
        x=memory_unit,
        y="rule",
        hue="module",
        dodge=False,
        data=plot_data,
        palette="Dark2",
    )
    plt.legend(loc="lower right", title="Module")

    xlabels = {
        "max_vms": "Virtual Memory Size (MB)",
        "max_uss": "Unique Set Size (MB)",
        "max_pss": "Proportional Set Size (MB)",
    }

    plt.title("Mean memory usage per sample")
    xlabel = xlabels[memory_unit]
    if gb:
        xlabel = xlabel.replace("MB", "GB")
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel("Rule name", rotation=0, fontsize=12)
    if outfile:
        plt.savefig(outfile, dpi=600, bbox_inches="tight")


def main(args):
    df = load_data(args.benchmarks_df)
    plots = (
        runtime_barplot_sum,
        runtime_barplot_errorbar,
        memory_barplot_sum,
        memory_barplot_errorbar,
    )
    for plot in plots:
        try:
            outfile = getattr(args, plot.__name__)
            plot(
                df,
                n_samples=args.n_samples,
                outfile=outfile,
                time_unit=args.time_unit,
                memory_unit=args.memory_unit,
                time_cutoff=args.time_cutoff,
                memory_cutoff=args.memory_cutoff,
                gb=args.gb,
            )
        except Exception as error:
            logging.info(f"Could not generate {outfile}. The following error occurred:")
            logging.info(error)
            logging.info("\n")
            if "memory" in plot.__name__:
                logging.info(
                    "It is possible that your OS does not support capture of memory usage."
                )
                logging.info("Therefore, only runtime plots will be generated.")
                logging.info(
                    "Memory plot files will be touched so the workflow will finish correctly.\n"
                )
                Path(outfile).touch()
            pass


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmarks_df")
    parser.add_argument("--n_samples")
    parser.add_argument("--time_unit", default="m")
    parser.add_argument("--memory_unit", default="max_vms")
    parser.add_argument("--time_cutoff", default=0)
    parser.add_argument("--memory_cutoff", default=0)
    parser.add_argument("--gb", action="store_true")
    parser.add_argument("--runtime_barplot_sum")
    parser.add_argument("--runtime_barplot_errorbar")
    parser.add_argument("--memory_barplot_sum")
    parser.add_argument("--memory_barplot_errorbar")
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
