#!/usr/bin/env python

import argparse
import logging
from collections import defaultdict

import pandas as pd
from Bio import SearchIO


def hmmer_to_df(hmmTbl, only_top_hit=False):
    """
    Takes a table from HMMER 3 and converts it to a Pandas Dataframe

    Adapted from https://stackoverflow.com/a/62021471
    """
    attribs = ["id", "evalue"]  # can add in more columns
    hits = defaultdict(list)

    # Prevent empty dataframe
    for attr in attribs + [
        "KO",
    ]:
        hits[attr]

    prev_hit_id = None
    ## open hmmTbl and extract hits
    with open(hmmTbl) as handle:
        for queryresult in SearchIO.parse(handle, "hmmer3-tab"):
            for hit in queryresult.hits:
                # Only record the top hit
                if only_top_hit and hit.id == prev_hit_id:
                    continue

                for attrib in attribs:
                    hits[attrib].append(getattr(hit, attrib))

                hits["KO"].append(queryresult.id)
                prev_hit_id = hit.id

    return pd.DataFrame.from_dict(hits)


def main(args):

    levels = ["Level1", "Level2", "Level3"]

    # load brite database
    brite_df = pd.read_csv(args.brite, sep="\t")

    # Loop over the HMMER tables
    counts_df = []
    for hmm_tbl in args.hmm_tbls:
        hmmer_df = hmmer_to_df(hmm_tbl, only_top_hit=False)

        # Select ids for rows with minimum e value
        idx_evalue_min = hmmer_df.groupby("id")["evalue"].idxmin()

        # Filter hmmer dataframe with these indexes
        hmmer_min_e_df = hmmer_df.loc[idx_evalue_min]
        brite_filtered = brite_df[brite_df["KO"].isin(hmmer_min_e_df.KO)]

        for level in levels:
            my_counts_df = (
                brite_filtered[level]
                .value_counts()
                .rename_axis("pathway")
                .reset_index(name="counts")
            )
            my_counts_df["level"] = level
            my_counts_df["hmm_tbl"] = hmm_tbl

            # Store in single dataframe
            counts_df = (
                my_counts_df
                if len(counts_df) == 0
                else pd.concat([counts_df, my_counts_df], ignore_index=True)
            )

    # Output the counts into text files
    for level in levels:
        output_filepath = f"{args.outprefix}_brite_{level}.tsv"

        df_for_pathways = (
            counts_df
            if args.consistent_pathways
            else counts_df[counts_df.level == level]
        )
        pathways_for_level = sorted(df_for_pathways.pathway.unique())

        logging.info(f"Writing to file {output_filepath}")
        with open(output_filepath, "w") as f:
            # Get pathways for this level so that we can have consistency in the output files even when the counts are zero
            if pathways_for_level:
                headers = ["Pathway"] + args.hmm_tbls
                f.write("\t".join(headers))
                f.write("\n")
                for pathway in pathways_for_level:
                    f.write(f"{pathway}")
                    for hmm_tbl in args.hmm_tbls:
                        filtered = counts_df[
                            (counts_df.pathway == pathway)
                            & (counts_df.level == level)
                            & (counts_df.hmm_tbl == hmm_tbl)
                        ]
                        count = filtered.counts.sum()
                        f.write(f"\t{count}")
                    f.write("\n")
            else:
                input_files = "\n\t\t".join(args.hmm_tbls)
                logging.debug(
                    f"The BRITE output for {level} appears to be empty. Please check the input files:\n{input_files}."
                )


def run(snakemake):
    # From https://stackoverflow.com/questions/16878315/what-is-the-right-way-to-treat-argparse-namespace-as-a-dictionary
    args = argparse.Namespace()
    args_dict = vars(args)

    for rule_input in ("hmm_tbls",):
        args_dict[rule_input] = snakemake.input[rule_input]

    for rule_param in ("brite", "consistent_pathways", "outprefix"):
        args_dict[rule_param] = snakemake.params[rule_param]

    main(args)


# Logging config must be at the beginning
logging.basicConfig(
    filename=str(snakemake.log),
    encoding="utf-8",
    level=logging.DEBUG,
    format="%(asctime)s %(message)s",
    datefmt="%m/%d/%Y %H:%M:%S",
)
logging.info(f"Starting script {__file__.split('/')[-1]}.")
logging.debug(f"Full script path: {__file__}")
run(snakemake)
logging.info(f"Done.")
