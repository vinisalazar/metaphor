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
        plt.savefig(outfile, bbox_inches="tight")


def runtime_barplot_errorbar(df, **kwargs):
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
        plt.savefig(outfile, bbox_inches="tight")


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
        print("Skipping plot, no memory stats recorded.")
        return
    fig, ax = plt.subplots(figsize=(4, 8))
    df_ = df.copy()
    df_[memory_unit]
    if gb:
        df_[memory_unit] = df_[memory_unit] / 2 ** 10
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
        plt.savefig(outfile, bbox_inches="tight")


def memory_barplot_errorbar(df, **kwargs):
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
        print("Skipping plot, no memory stats recorded.")
        return
    fig, ax = plt.subplots(figsize=(4, 8))
    df_ = df.copy()
    if gb:
        df_[memory_unit] = df_[memory_unit] / 2 ** 10
    plot_data = df_.query(f"{memory_unit} > {cutoff}").sort_values(
        memory_unit, ascending=False
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
        plt.savefig(outfile, bbox_inches="tight")


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
                outfile=outfile,
                time_unit=args.time_unit,
                memory_unit=args.memory_unit,
                time_cutoff=args.time_cutoff,
                memory_cutoff=args.memory_cutoff,
                gb=args.gb,
            )
        except Exception as error:
            print(f"Could not generate {outfile}. The following error occurred:")
            print(error)
            if "memory" in plot.__name__:
                print(
                    "It is possible that your OS does not support capture of memory usage."
                )
                print("Therefore, only runtime plots will be generated.\n")
            pass


def parse_snakemake_args(snakemake):
    args = argparse.Namespace()
    args_dict = vars(args)

    for directive in "input", "output", "params":
        try:
            for k, v in getattr(snakemake, directive).items():
                args_dict[k] = v
        except AttributeError:
            pass

    return args


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmarks_df")
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
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%m/%d/%Y %H:%M:%S",
    )
    logging.info(f"Starting script '{__file__.split('/')[-1]}'.")
    logging.debug(f"Full script path: '{__file__}'.")
    if "snakemake" in locals():
        logging.basicConfig(filename=str(snakemake.log))
        args = parse_snakemake_args(snakemake)
    else:
        args = parse_args()
    try:
        main(args)
        logging.info("Done.")
    except Exception as e:
        logging.error(e)
