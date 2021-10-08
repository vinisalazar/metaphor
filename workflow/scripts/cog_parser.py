"""
Parses tabular output file with COG annotations.
"""

import sys
import pandas as pd

from pathlib import Path
from functools import lru_cache


def main(dmnd_out, cog_csv, fun_tab, def_tab):
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
    print(f"Loading annotation data: '{dmnd_out}'.")
    df = pd.read_csv(dmnd_out, sep="\t").drop_duplicates("qseqid")
    df["Protein ID"] = (
        df["sseqid"].str[::-1].str.split("_", 1).apply(lambda l: ".".join(l)).str[::-1]
    )
    print(f"Loaded {len(df)} records.")

    print(f"Loading main COG dataframe: '{cog_csv}'.")
    cog_csv = pd.read_csv(
        cog_csv, names=cog_csv_names, index_col="Protein ID"
    ).drop_duplicates()
    print(f"Loading additional dataframes: '{def_tab}', '{fun_tab}'.")
    def_tab = pd.read_csv(def_tab, sep="\t", names=def_tab_names, index_col=0)
    fun_tab = pd.read_csv(fun_tab, sep="\t", names=fun_tab_names, index_col=0)

    # Merge data
    print("Merging dataframes.")
    merged_df = df.merge(cog_csv, left_on="Protein ID", right_index=True).reset_index(
        drop=True
    )
    del cog_csv
    del df
    merged_df = merged_df.merge(
        def_tab, left_on="COG ID", right_index=True
    ).drop_duplicates("qseqid")
    print(f"{len(merged_df)} records after merging.")

    # Formatting
    print("Applying final formatting.")

    @lru_cache(128)
    def get_function(func_id):
        return fun_tab.loc[func_id, "Description"]

    merged_df["COG categories"] = merged_df["COG functional category"].apply(
        lambda s: [[get_function(c) for c in list(string)] for string in s]
    )
    merged_df["COG categories"] = merged_df["COG categories"].apply(
        lambda list_: [item for sublist in list_ for item in sublist]
    )

    cat_outfile = Path(dmnd_out.replace("_dmnd.out", "_COG_categories")).with_suffix(
        ".csv"
    )
    cat_counts = merged_df["COG categories"].explode().value_counts()
    cat_counts = pd.concat(
        (cat_counts, cat_counts / cat_counts.sum()), axis=1
    ).reset_index()
    cat_counts.columns = "COG category", "absolute", "relative"
    breakpoint()
    cat_counts.to_csv(cat_outfile, index=False)
    print(f"Wrote category counts to '{cat_outfile}.'")

    cog_counts = merged_df["COG ID"].value_counts().reset_index()

    @lru_cache(128)
    def get_cog_name(cog_name):
        return def_tab.loc[cog_name, "COG name"]

    cog_counts["COG name"] = cog_counts["index"].apply(get_cog_name)
    cog_counts.columns = "COG code", "absolute", "COG name"
    cog_counts["relative"] = cog_counts["absolute"] / cog_counts["absolute"].sum()
    cog_counts_outfile = Path(dmnd_out.replace("_dmnd.out", "_COG_codes")).with_suffix(
        ".csv"
    )
    cog_counts.to_csv(cog_counts_outfile, index=False)
    print(f"Wrote category counts to '{cog_counts_outfile}.'")


def parse_snakemake_args(snakemake):
    args = argparse.Namespace()
    args_dict = vars(args)

    for rule_input in ("dmnd_out", "cog_csv", "fun_tab", "def_tab"):
        args_dict[rule_input] = snakemake.input[rule_input]

    return args


if "snakemake" in locals():
    args = run(snakemake)
    main(*args)

if __name__ == "__main__":
    main(*sys.argv[1:])
