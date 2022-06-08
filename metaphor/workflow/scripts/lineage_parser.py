import argparse
import logging
import pandas as pd


rankedlineage_names = (
    "taxid tax_name species genus family order class phylum kingdom domain nan".split()
)
ranks = rankedlineage_names[2:-1]
ranks.remove("kingdom")


def main(args):
    tax = pd.read_csv(
        args.tax_out,
        sep="\t",
        index_col=0,
    )
    logging.info("Loading NCBI Taxonomy lineage data.")
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
        outfile = getattr(args, rank)
        rank_df.to_csv(outfile, sep="\t")
        logging.info(f"Wrote {len(rank_df)} records to '{outfile}'.")


def parse_snakemake_args(snakemake):
    args = argparse.Namespace()
    args_dict = vars(args)

    for directive in "input", "output", "params":
        try:
            # The following block prevents conflict with Python's reserved keyword 'class'
            # An old problem for Pythonistas that work with taxonomy :)
            # See the rule 'lineage_parser' output directive to see class is spelled with a 'k'
            for k, v in getattr(snakemake, directive).items():
                if k == "klass":
                    k = k.replace("k", "c")
                args_dict[k] = v
        except AttributeError:
            pass

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


if __name__ == "__main__":
    # The driver function is standardized across scripts in this workflow
    # Please check the workflow/scripts/utils.py module for reference
    from utils import driver

    if "snakemake" not in locals():
        snakemake = None
        parse_args_fn = parse_args
    else:
        parse_args_fn = parse_snakemake_args
    driver(main, snakemake, __file__, parse_args_fn)
