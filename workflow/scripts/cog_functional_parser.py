"""
Parses tabular output file with COG annotations.

Because the output of Diamond can be fairly large, and the COG database files are also not small,
the current version of this script only uses the required columns from both files. This greatly
reduces memory usage, as lots of unnecessary columns are required.

In the future, a configurable option may be added to enable the output of all columns.
"""

import logging
import argparse
from pathlib import Path
from functools import lru_cache
from utils import cog_csv_names

import pandas as pd


def main(args):
    (dmnd_out, cog_csv, fun_tab, def_tab, pathways_out, codes_out, categories_out) = (
        args.dmnd_out,
        args.cog_csv,
        args.fun_tab,
        args.def_tab,
        args.pathways_out,
        args.codes_out,
        args.categories_out,
    )

    for file in (dmnd_out, cog_csv, fun_tab, def_tab):
        assert Path(file).exists(), f"File '{file}' was not found."

    # The following functions use the merged_df variable
    # to write different informations, respectively the COG pathways, categories,
    # and codes.
    #
    # Notice that each function starts and ends with a logging call,
    # with the latter logging call being followed by deletion of the
    # DataFrame that was just processed.
    merged_df, def_tab = create_merged_df(dmnd_out, cog_csv, def_tab)
    merged_df = write_cog_pathways(merged_df, pathways_out)
    merged_df = write_cog_categories(merged_df, fun_tab, categories_out)
    merged_df = write_cog_codes(merged_df, def_tab, codes_out)


def create_merged_df(dmnd_out, cog_csv, def_tab):
    """
    Create dataframe merging COG information with Diamond output.
    """

    def_tab_names = [
        "COG ID",
        "COG functional category",
        "COG name",
        "Gene",
        "Functional pathway",
        "PubMed ID",
        "PDB ID",
    ]

    # Load data
    logging.info(f"Loading annotation data: '{dmnd_out}'.")
    df = pd.read_csv(
        dmnd_out,
        sep="\t",
        usecols=[
            "qseqid",
            "sseqid",
            "bitscore",
        ],
    )

    # Keep the best score
    df = df.sort_values("bitscore", ascending=False)
    df = df.drop_duplicates("qseqid")

    df["Protein ID"] = (
        df["sseqid"].str[::-1].str.split("_", 1).apply(lambda l: ".".join(l)).str[::-1]
    )
    logging.info(f"Loaded {len(df)} records.")

    logging.info(f"Loading main COG dataframe: '{cog_csv}'.")
    cog_csv = pd.read_csv(
        cog_csv,
        names=cog_csv_names,
        index_col="Protein ID",
        usecols=["Protein ID", "COG ID"],
    ).drop_duplicates()

    # Merge COG database information (cog_csv, def_tab) with Diamond output (df)
    logging.info("Merging dataframes.")
    merged_df = df.merge(cog_csv, left_on="Protein ID", right_index=True).reset_index(
        drop=True
    )
    def_tab = load_dataframe(
        def_tab,
        sep="\t",
        names=def_tab_names,
        index_col=0,
        usecols=[
            "COG ID",
            "COG functional category",
            "COG name",
            "Functional pathway",
        ],
    )
    merged_df = merged_df.merge(
        def_tab, left_on="COG ID", right_index=True
    ).drop_duplicates("qseqid")
    logging.info(f"{len(merged_df)} records after merging.")
    return merged_df, def_tab


# COG pathways
def write_cog_pathways(merged_df, pathways_out):
    logging.info("Writing COG pathways.")
    print(merged_df.head())
    print(merged_df.columns)
    pathways = merged_df["Functional pathway"].fillna("Unknown").value_counts()
    pathways.index.name = "Functional pathway"
    pathways = pd.concat((pathways, pathways / pathways.sum()), axis=1)
    pathways.columns = "absolute", "relative"
    pathways.to_csv(pathways_out, sep="\t")
    logging.info(f"Wrote {len(pathways)} rows to '{pathways_out}'.")
    return merged_df


# COG categories
def write_cog_categories(merged_df, fun_tab, categories_out):
    fun_tab_names = [
        "ID",
        "RGB",
        "Description",
    ]
    logging.info("Writing COG categories.")
    fun_tab = load_dataframe(fun_tab, sep="\t", names=fun_tab_names, index_col=0)

    @lru_cache(1024)
    def get_COG_function(func_id):
        return fun_tab.loc[func_id, "Description"]

    merged_df["COG categories"] = merged_df["COG functional category"].apply(
        lambda s: [[get_COG_function(c) for c in list(string)] for string in s]
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
    logging.info(f"Wrote {len(cat_counts)} rows to '{categories_out}.'")
    return merged_df


def write_cog_codes(merged_df, def_tab, codes_out):
    # COG codes
    logging.info("Writing COG codes.")
    cog_counts = merged_df["COG ID"].value_counts().reset_index()

    @lru_cache(1024)
    def get_COG_name(cog_name):
        return def_tab.loc[cog_name, "COG name"]

    cog_counts["COG name"] = cog_counts["index"].apply(get_COG_name)
    cog_counts.columns = "COG code", "absolute", "COG name"
    cog_counts["relative"] = cog_counts["absolute"] / cog_counts["absolute"].sum()
    cog_counts = cog_counts[["COG code", "COG name", "absolute", "relative"]]
    cog_counts.to_csv(codes_out, index=False, sep="\t")
    logging.info(f"Wrote {len(cog_counts)} rows to '{codes_out}.'")
    return merged_df


def load_dataframe(file, **kwargs):
    logging.info(f"Loading dataframe: '{file}'.")
    try:
        df_ = pd.read_csv(file, **kwargs)
    except UnicodeDecodeError:
        logging.info("Couldn't load as utf-8. Using 'latin-1' encoding.")
        df_ = pd.read_csv(file, **kwargs, encoding="latin-1")

    return df_


def parse_args():
    # TODO: improve these arguments and add parser description
    parser = argparse.ArgumentParser()
    parser.add_argument("--dmnd_out")
    parser.add_argument("--cog_csv")
    parser.add_argument("--fun_tab")
    parser.add_argument("--def_tab")
    parser.add_argument("--categories_out")
    parser.add_argument("--codes_out")
    parser.add_argument("--pathways_out")
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
