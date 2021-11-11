"""
Parses tabular output file with COG annotations.
"""

import sys
import logging
import argparse
import subprocess
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

    # Save this variable for later
    cog_dir = Path(cog_csv).parent.joinpath("fasta")

    logging.info(f"Loading main COG dataframe: '{cog_csv}'.")
    cog_csv = pd.read_csv(
        cog_csv, names=cog_csv_names, index_col="Protein ID"
    ).drop_duplicates()

    # Merge data
    logging.info("Merging dataframes.")
    merged_df = df.merge(cog_csv, left_on="Protein ID", right_index=True).reset_index(
        drop=True
    )
    def_tab = load_dataframe(def_tab, sep="\t", names=def_tab_names, index_col=0)
    merged_df = merged_df.merge(
        def_tab, left_on="COG ID", right_index=True
    ).drop_duplicates("qseqid")
    logging.info(f"{len(merged_df)} records after merging.")

    # Write COG categories
    logging.info("Applying final formatting.")
    fun_tab = load_dataframe(fun_tab, sep="\t", names=fun_tab_names, index_col=0)

    # Cache COG function lookup to go faster
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
    logging.info(f"Wrote category counts to '{categories_out}.'")

    # Write COG codes
    cog_counts = merged_df["COG ID"].value_counts().reset_index()

    @lru_cache(1024)
    def get_cog_name(cog_name):
        return def_tab.loc[cog_name, "COG name"]

    cog_counts["COG name"] = cog_counts["index"].apply(get_cog_name)
    cog_counts.columns = "COG code", "absolute", "COG name"
    cog_counts["relative"] = cog_counts["absolute"] / cog_counts["absolute"].sum()
    cog_counts = cog_counts[["COG code", "COG name", "absolute", "relative"]]
    cog_counts.to_csv(codes_out, index=False, sep="\t")
    logging.info(f"Wrote code counts to '{codes_out}.'")

    # Write taxonomies
    @lru_cache(None)
    def get_taxid_and_name(protein_id, cog_code):
        try:
            suffix = ".tsv.gz"
            cog_tsv_file = cog_dir.joinpath(Path(cog_code).with_suffix(suffix))
            if not cog_tsv_file.exists():
                suffix = ".tsv"
                cog_tsv_file = cog_dir.joinpath(Path(cog_code).with_suffix(suffix))
            protid, plen, taxid, name, footprint = subprocess.getoutput(
                f"zgrep '{protein_id}' {cog_tsv_file}"
            ).split("\t")
            fmt_name = name.replace("_", " ") + f" {footprint}"
            return taxid, fmt_name
        except:
            logging.info(f"Couldn't find tax ID for {cog_code} : {protein_id}.")
            logging.info(
                "Please check you have all correct files in the fasta/ directory of the COG database."
            )
            # raise
            return None, None

    logging.info("Matching Protein IDs to Tax IDs.")
    merged_df[["taxid", "taxname"]] = merged_df.apply(
        lambda row: get_taxid_and_name(row["Protein ID"], row["COG ID"]), axis=1
    )
    # breakpoint()


def load_dataframe(file, **kwargs):
    logging.info(f"Loading dataframe: '{file}'.")
    try:
        df_ = pd.read_csv(file, **kwargs)
    except UnicodeDecodeError:
        logging.info("Couldn't load as utf-8. Using 'latin-1' encoding.")
        df_ = pd.read_csv(file, **kwargs, encoding="latin-1")

    return df_


def parse_snakemake_args(snakemake):
    args = argparse.Namespace()
    args_dict = vars(args)

    for rule_input in ("dmnd_out",):
        args_dict[rule_input] = snakemake.input[rule_input]

    for rule_param in ("cog_csv", "fun_tab", "def_tab"):
        args_dict[rule_param] = snakemake.params[rule_param]

    for rule_output in ("categories_out", "codes_out"):
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
