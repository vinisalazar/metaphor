import logging
import pandas as pd
from functools import lru_cache


def main():

    lin = pd.read_csv(
        "/Users/vini/Bio/MGP/databases/taxonomy/taxidlineage.dmp",
        sep="|",
        index_col=0,
        names=["NCBI_taxname", "lineage", ""],
    ).dropna(how="all", axis=1)
    lin = lin.applymap(lambda x: x.strip() if isinstance(x, str) else x)
    tax = pd.read_csv(
        "/Users/vini/Bio/MGP/databases/COG2020/H_S001_tax_out.tsv",
        sep="\t",
        index_col=0,
    )
    tax = tax.join(lin)

    nodes_names = [
        "taxid",
        "parent_tax_id",
        "rank",
    ]

    nodes = pd.read_csv(
        "/Users/vini/Bio/MGP/databases/taxonomy/nodes.dmp",
        sep="|",
        names=nodes_names,
        usecols=["taxid", "parent_tax_id", "rank"],
        index_col=0,
    )
    nodes = nodes.applymap(lambda x: x.strip() if isinstance(x, str) else x)

    @lru_cache(1024)
    def get_rank(taxid):
        try:
            return nodes.loc[int(taxid), "rank"]
        except (KeyError, TypeError):
            return None

    ranks = "superkingdom phylum class order family genus species".split()

    def parse_lineage(lineage):
        output = dict()
        if isinstance(lineage, str):
            lineage = lineage.split()

        for taxid in lineage:
            if (rank := get_rank(taxid)) in ranks:
                output[rank] = taxid

        return pd.Series(output)

    tax = tax.join(tax["lineage"].apply(parse_lineage))
    tax.to_csv("parsed_tax.csv")


if __name__ == "__main__":
    main()
