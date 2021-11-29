import argparse
import logging
import pandas as pd


rankedlineage_names = (
    "taxid tax_name species genus family order class phylum kingdom domain nan".split()
)
ranks = rankedlineage_names[2:-1]


def main(args):
    tax = pd.read_csv(
        args.tax_out,
        sep="\t",
        index_col=0,
    )
    rankedlineage = pd.read_csv(
        args.rankedlineage,
        sep="|",
        index_col=0,
        names=rankedlineage_names,
    )
    rankedlineage = rankedlineage.applymap(
        lambda x: x.strip() if isinstance(x, str) else x
    ).dropna(how="all", axis=1)
    tax = tax.join(rankedlineage)

    for rank in ranks:
        rank_df = tax.groupby(rank).sum()
        rank_df = pd.concat((rank_df, rank_df / rank_df.sum()), axis=1)
        rank_df.columns = "absolute", "relative"
        rank_df.to_csv(vars(args).get(rank), sep="\t")


def parse_snakemake_args(snakemake):
    args = argparse.Namespace()
    args_dict = vars(args)

    for rule_input in ("tax_out", "rankedlineage"):
        args_dict[rule_input] = snakemake.input[rule_input]

    for rank in ranks:
        # The following block prevents conflict with Python's reserved keyword 'class'
        # An old problem for Pythonistas that work with taxonomy :)
        # See the rule 'lineage_parser' output directive to see class is spelled with a 'k'
        if (rank_ := rank) == "class":
            rank_ = rank.replace("c", "k")
        else:
            pass
        args_dict[rank] = snakemake.output[rank_]

    args_dict["threads"] = snakemake.threads

    return args


def parse_args():
    # TODO: improve these arguments and add parser description
    parser = argparse.ArgumentParser(
        description="Parses output of the COG_parser script into different hierarchical levels for each tax_id."
    )
    parser.add_argument(
        "--tax_out", help="Output file of Metaphor's COG_parser script."
    )
    parser.add_argument(
        "--rankedlineage",
        help="NCBI Taxonomy rankedlineage.dmp file. Available at: https://ftp.ncbi.nih.gov/pub/taxonomy/",
    )
    for rank in ranks:
        parser.add_argument(
            f"--{rank}", help=f"Output file for rank '{rank}.' Extension is TSV."
        )
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
    logging.debug(f"Full script path: {__file__}.")
    try:
        args = parse_snakemake_args(snakemake)
        main(args)
        logging.info("Done.")
    except Exception as e:
        logging.error(e)
elif __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(message)s",
        datefmt="%m/%d/%Y %H:%M:%S",
    )
    logging.info(f"Starting script {__file__.split('/')[-1]}.")
    logging.debug(f"Full script path: {__file__}.")
    try:
        args = parse_args()
        main(args)
        logging.info("Done.")
    except Exception as e:
        logging.error(e)
