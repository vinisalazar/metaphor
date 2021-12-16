#!/usr/bin/env python

"""
Concatenates benchmarks generated by the workflow into a single table.

Usage:

    python cat_benchmarks.py <benchmarks-dir>
"""


import argparse
import logging
import traceback
from glob import glob
from pathlib import Path
import pandas as pd


def create_benchmark_df(file, benchmarks_dir):
    """
    Returns a DataFrame from benchmark files. Contains no header.
    """
    df = pd.read_csv(file, sep="\t", skiprows=1, header=None)
    df["file"] = file.replace(benchmarks_dir + "/", "")
    rules = df["file"].str.split("/", expand=True)
    rules.iloc[:, -1] = rules.iloc[:, -1].apply(lambda s: str(Path(s).stem))
    if rules.shape[1] < 3:
        rules[2] = None
    df = pd.concat((df, rules), axis=1)
    return df


def cat_benchmark_dfs(list_of_dfs):
    df = pd.concat(list_of_dfs)
    df.columns = "s h:m:s max_rss max_vms max_uss max_pss io_in io_out mean_load cpu_time file module rule sample".split()
    return df


def get_benchmark_files(benchmarks_dir):
    files = glob(benchmarks_dir + "/*/*.txt") + glob(benchmarks_dir + "/*/*/*.txt")
    return files


def main(args):
    benchmarks_dir, outfile = args.benchmarks_dir, args.outfile
    files = get_benchmark_files(benchmarks_dir)
    list_of_dfs = [create_benchmark_df(file, benchmarks_dir) for file in files]
    df = cat_benchmark_dfs(list_of_dfs)

    print(f"Writing {len(files)} benchmarks to {outfile}.\n")
    df.to_csv(outfile, index=False)


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
    parser.add_argument("--benchmarks_dir")
    parser.add_argument("--outfile")
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
        logging.error(traceback.format_exc())
