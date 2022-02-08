#!/usr/bin/env python
import argparse
import logging
from pathlib import Path
from subprocess import Popen, PIPE

import pandas as pd
from matplotlib import pyplot as plt


def run_seqstats(fasta):
    """
    Runs seqstats and returns dictionary of metrics.

    Formats metrics to numeric type as well.

    fasta: path to FASTA file

    returns: dict
    """
    seqstats = Popen(f"seqstats {fasta}", shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = seqstats.communicate()
    stdout, stderr = (
        stdout.decode("utf-8"),
        stderr.decode("utf-8"),
    )
    assert (
        "command not found" not in stderr
    ), "seqstats was not found. Please make sure it's installed and on PATH."
    metrics = stdout.strip().split("\n")
    metrics = [i.split("\t") for i in metrics]
    metrics = {i[0][:-1]: i[-1].replace("bp", "").strip() for i in metrics}

    # Convert to numeric values
    metrics = dict(
        (k, float(v)) if "." in v else (k, int(v)) for k, v in metrics.items()
    )
    return metrics


def metrics_to_df(fastas):
    samples = [str(Path(file).name).replace(".contigs.fa", "") for file in fastas] * 3
    metrics = [run_seqstats(file) for file in fastas] * 3
    df = pd.DataFrame(metrics)
    df.index = samples
    rename_dict = {
        "Total n": "# contigs",
        "Total seq": "# bp",
        "Avg. seq": "Avg. length",
        "Median seq": "Median length",
        "N 50": "N50",
        "Min seq": "Min. length",
        "Max seq": "Max. length",
    }
    df = df.rename(columns=rename_dict)
    outfile = "output/assembly/assembly_report.tsv"
    df.to_csv(outfile, sep="\t")
    logging.info(f"Generated assembly report: '{outfile}'.")
    return df


def plot_column(df, column, outfile):
    fig, ax = plt.subplots(figsize=(4, 4))
    df[column].plot(kind="barh", ax=ax)
    ax.set_xlabel(column)
    plt.ticklabel_format(axis="x", style="sci", scilimits=(0, 4))
    fig.savefig(outfile, bbox_inches="tight")


def plot_columns(df, outdir="./"):
    for column in df.columns:
        outfile = Path(outdir).joinpath(
            column.replace(" ", "_").replace(".", "").replace("#", "n").lower() + ".pdf"
        )
        plot_column(df, column, outfile)
        logging.info(f"Generated plot: '{outfile}'.")


def main(args):
    df = metrics_to_df(args.fastas)
    plot_columns(df)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastas", nargs="+")
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
