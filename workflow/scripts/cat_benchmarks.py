#!/usr/bin/env python

"""
Concatenates benchmarks generated by the workflow into a single table.

Usage:

    python cat_benchmarks.py <benchmarks-dir>
"""


import sys
from glob import glob
import pandas as pd
    

def create_benchmark_df(file, benchmarks_dir):
    """
    Returns a DataFrame from benchmark files. Contains no header.
    """
    df = pd.read_csv(file, sep="\t", skiprows=1, header=None)
    df["file"] = file.replace(benchmarks_dir + "/", "")
    return df


def cat_benchmark_dfs(list_of_dfs):
    df = pd.concat(list_of_dfs)
    df.columns = "s h:m:s max_rss max_vms max_uss max_pss io_in io_out mean_load cpu_time file".split()
    return df


def get_benchmark_files(benchmarks_dir):
    files = glob(benchmarks_dir + "/*/*.txt") + glob(benchmarks_dir + "/*/*/*.txt")
    return files


def main(benchmarks_dir):
    files = get_benchmark_files(benchmarks_dir)
    list_of_dfs = [create_benchmark_df(file, benchmarks_dir) for file in files]
    df = cat_benchmark_dfs(list_of_dfs)
    try:
        outfile = sys.argv[2]
    except IndexError:
        outfile = benchmarks_dir + "/all_benchmarks.csv"

    print(f"Writing {len(files)} benchmarks to {outfile}.")
    df.to_csv(outfile, index=False)


if __name__ == "__main__":
    main(sys.argv[1])