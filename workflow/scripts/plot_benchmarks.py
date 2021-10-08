import sys
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


def runtime_barplot_sum(df, outfile=None, time_unit="m", cutoff=5):
    """
    Generates barplot of sum of rules' runtimes.
    """
    fig, ax = plt.subplots(figsize=(4, 8))
    plot_data = (df.groupby(by=["module", "rule"])[time_unit].sum()).reset_index().sort_values(time_unit, ascending=False).query(f"{time_unit} > {cutoff}")
    sns.barplot(x=time_unit, y="rule", hue="module", dodge=False, data=plot_data, palette="Dark2")
    plt.legend(loc="lower right", title="Module")

    xlabels = {
        "h": "(hours)",
        "m": "(minutes)",
        "s": "(seconds)"
    }

    plt.title("Total rule runtime")
    ax.set_xlabel("Runtime " + xlabels[time_unit], fontsize=12)
    ax.set_ylabel("Rule name", rotation=0, fontsize=12)
    if outfile:
        plt.savefig(outfile, bbox_inches="tight")


def runtime_barplot_errorbar(df, outfile=None, time_unit="m", cutoff=1):
    """
    Generates barplot of rules' runtimes for each sample.

    The time-per-sample is corrected in rules that run once for all samples.
    """
    fig, ax = plt.subplots(figsize=(4, 8))
    df[time_unit + "_corrected"] = df.apply(lambda row: row['m'] / df['sample'].str[:6].value_counts().shape[0] if not isinstance(row['sample'], str) else row['m'], axis=1)
    plot_data = df.sort_values(time_unit + "_corrected", ascending=False).query(f"{time_unit} > {cutoff}")
    sns.barplot(x=time_unit + "_corrected", y="rule", hue="module", dodge=False, data=plot_data, palette="Dark2")
    plt.legend(loc="lower right", title="Module")

    xlabels = {
        "h": "(hours)",
        "m": "(minutes)",
        "s": "(seconds)"
    }

    plt.title("Mean runtime per sample")
    ax.set_xlabel("Runtime " + xlabels[time_unit], fontsize=12)
    ax.set_ylabel("Rule name", rotation=0, fontsize=12, position=(0, .7))
    if outfile:
        plt.savefig(outfile, bbox_inches="tight")


def memory_barplot_sum(df, outfile=None, memory_unit="max_vms", cutoff=1, gb=True):

    """
    Generates barplot of sum of rules' memory usage.
    """
    if df[memory_unit].all() == '-':
        print("Skipping plot, no memory stats recorded.")
        return
    fig, ax = plt.subplots(figsize=(4, 8))
    df_ = df.copy()
    if gb:
        df_[memory_unit] = df_[memory_unit] / 2**10
    plot_data = (df_.groupby(by=["module", "rule"])[memory_unit].sum()).reset_index().sort_values(memory_unit, ascending=False).query(f"{memory_unit} > {cutoff}")
    sns.barplot(x=memory_unit, y="rule", hue="module", dodge=False, data=plot_data, palette="Dark2")
    plt.legend(loc="lower right", title="Module")

    xlabels = {
        "max_vms": "Virtual Memory Size (MB)",
        "max_uss": "Unique Set Size (MB)",
        "max_pss": "Proportional Set Size (MB)"
    }

    plt.title("Sum of rules' memory usage")
    xlabel = xlabels[memory_unit]
    if gb:
        xlabel = xlabel.replace("MB", "GB")
    ax.set_xlabel(xlabel, labelpad=15, fontsize=12)
    ax.set_ylabel("Rule name", rotation=0, fontsize=12)
    if outfile:
        plt.savefig(outfile, bbox_inches="tight")


def memory_barplot_errorbar(df, outfile=None, memory_unit="max_vms", cutoff=0, gb=True):
    """
    Generates barplot of sum of rules' memory usage.
    """
    if df[memory_unit].all() == '-':
        print("Skipping plot, no memory stats recorded.")
        return
    fig, ax = plt.subplots(figsize=(4, 8))
    df_ = df.copy()
    if gb:
        df_[memory_unit] = df_[memory_unit] / 2**10
    plot_data = df_.query(f"{memory_unit} > {cutoff}").sort_values(memory_unit, ascending=False)
    sns.barplot(x=memory_unit, y="rule", hue="module", dodge=False, data=plot_data, palette="Dark2")
    plt.legend(loc="lower right", title="Module")

    xlabels = {
        "max_vms": "Virtual Memory Size (MB)",
        "max_uss": "Unique Set Size (MB)",
        "max_pss": "Proportional Set Size (MB)"
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
    df = load_data(args[1])
    plots = (
        runtime_barplot_sum,
        runtime_barplot_errorbar,
        memory_barplot_sum,
        memory_barplot_errorbar,
    )
    for plot in plots:
        try:
            outfile = Path(args[1]).parent.joinpath(plot.__name__ + ".png")
            plot(df, outfile)
        except Exception as error:
            print(f"Could not generate {outfile}. The following error occurred:\n")
            print(error)
            pass


if __name__ == "__main__":
    args = sys.argv
    main(args)
