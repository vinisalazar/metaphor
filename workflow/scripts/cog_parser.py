"""
Parses tabular output file with COG annotations.
"""

import sys
import logging
import argparse
from pathlib import Path
from functools import lru_cache

import pandas as pd


def main(args):

    dmnd_out, cog_csv, fun_tab, def_tab, categories_out, codes_out = (
        args.dmnd_out,
        args.cog_csv,
        args.fun_tab,
        args.def_tab,
        args.categories_out,
        args.codes_out,
    )

    for file in (dmnd_out, cog_csv, fun_tab, def_tab):
        assert Path(file).exists(), f"File '{file}' was not found."

    # Column names
    cog_csv_names = [
        "Gene ID (GenBank or ad hoc)",
        "NCBI Assembly ID",
        "Protein ID",
        "Protein length",
        "COG footprint coordinates",
        "Length of the COG footprint on the proteins",
        "COG ID",
        "reserved",
        "COG membership class",
        "PSI-BLAST bit score",
        "PSI-BLAST e-value",
        "COG profile length",
        "Protein footprint coordinates on the COG profile",
    ]

    def_tab_names = [
        "COG ID",
        "COG functional category",
        "COG name",
        "Gene",
        "Functional pathway",
        "PubMed ID",
        "PDB ID",
    ]
    fun_tab_names = [
        "ID",
        "RGB",
        "Description",
    ]

    # Load data
    logging.info(f"Loading annotation data: '{dmnd_out}'.")
    df = pd.read_csv(dmnd_out, sep="\t").drop_duplicates("qseqid")
    df["Protein ID"] = (
        df["sseqid"].str[::-1].str.split("_", 1).apply(lambda l: ".".join(l)).str[::-1]
    )
    logging.info(f"Loaded {len(df)} records.")

    logging.info(f"Loading main COG dataframe: '{cog_csv}'.")
    cog_csv = pd.read_csv(
        cog_csv, names=cog_csv_names, index_col="Protein ID"
    ).drop_duplicates()
    logging.info(f"Loading additional dataframes: '{def_tab}', '{fun_tab}'.")

    try:
        def_tab = pd.read_csv(def_tab, sep="\t", names=def_tab_names, index_col=0)
        fun_tab = pd.read_csv(fun_tab, sep="\t", names=fun_tab_names, index_col=0)
    except UnicodeDecodeError:
        logging.info("Couldn't load as utf-8. Using 'latin-1' encoding.")
        def_tab = pd.read_csv(
            def_tab, sep="\t", names=def_tab_names, index_col=0, encoding="latin-1"
        )
        fun_tab = pd.read_csv(
            fun_tab, sep="\t", names=fun_tab_names, index_col=0, encoding="latin-1"
        )

    # Merge data
    logging.info("Merging dataframes.")
    merged_df = df.merge(cog_csv, left_on="Protein ID", right_index=True).reset_index(
        drop=True
    )
    del cog_csv
    del df
    merged_df = merged_df.merge(
        def_tab, left_on="COG ID", right_index=True
    ).drop_duplicates("qseqid")
    logging.info(f"{len(merged_df)} records after merging.")

    # Formatting
    logging.info("Applying final formatting.")

    # Cache function lookup to go faster
    @lru_cache(128)
    def get_function(func_id):
        return fun_tab.loc[func_id, "Description"]

    merged_df["COG categories"] = merged_df["COG functional category"].apply(
        lambda s: [[get_function(c) for c in list(string)] for string in s]
    )
    merged_df["COG categories"] = merged_df["COG categories"].apply(
        lambda list_: [item for sublist in list_ for item in sublist]
    )

    cat_counts = merged_df["COG categories"].explode().value_counts()
    cat_counts = pd.concat(
        (cat_counts, cat_counts / cat_counts.sum()), axis=1
    ).reset_index()
    cat_counts.columns = "COG category", "absolute", "relative"
    cat_counts.to_csv(categories_out, index=False, sep="\t")
    logging.info(f"Wrote category counts to '{categories_out}.'")

    cog_counts = merged_df["COG ID"].value_counts().reset_index()

    @lru_cache(128)
    def get_cog_name(cog_name):
        return def_tab.loc[cog_name, "COG name"]

    cog_counts["COG name"] = cog_counts["index"].apply(get_cog_name)
    cog_counts.columns = "COG code", "absolute", "COG name"
    cog_counts["relative"] = cog_counts["absolute"] / cog_counts["absolute"].sum()
    cog_counts = cog_counts[["COG code", "COG name", "absolute", "relative"]]
    cog_counts.to_csv(codes_out, index=False, sep="\t")
    logging.info(f"Wrote category counts to '{codes_out}.'")


def parse_snakemake_args(snakemake):
    args = argparse.Namespace()
    args_dict = vars(args)

    for rule_input in ("dmnd_out",):
        args_dict[rule_input] = snakemake.input[rule_input]

    for rule_param in ("cog_csv", "fun_tab", "def_tab"):
        args_dict[rule_param] = snakemake.params[rule_param]

    for rule_output in ("cog_csv", "fun_tab", "def_tab"):
        args_dict[rule_output] = snakemake.output[rule_output]

    return args


def parse_args():
    # TODO: improve these arguments and add parser description
    parser = argparse.ArgumentParser()
    parser.add_argument("--dmnd_out")
    parser.add_argument("--cog_csv")
    parser.add_argument("--fun_tab")
    parser.add_argument("--def_tab")
    parser.add_argument("--categories_out")
    parser.add_argument("--codes_out")
    args = parser.parse_args()
    return args


if "snakemake" in locals():
    logging.basicConfig(
        filename=str(snakemake.log),
        encoding="utf-8",
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%m/%d/%Y %H:%M:%S",
    )
    logging.info(f"Starting script {__file__.split('/')[-1]}.")
    logging.debug(f"Full script path: {__file__}")
    args = parse_snakemake_args(snakemake)
    main(args)
    logging.info("Done.")
elif __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%m/%d/%Y %H:%M:%S",
    )
    logging.info(f"Starting script '{__file__.split('/')[-1]}'.")
    logging.debug(f"Full script path: '{__file__}'.")
    args = parse_args()
    main(args)
    logging.info("Done.")
