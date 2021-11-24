import logging
import pandas as pd
from functools import lru_cache


def main():

    tax = pd.read_csv(
        "/Users/vini/Bio/MGP/databases/COG2020/H_S001_tax_out.tsv",
        sep="\t",
        index_col=0,
    )
    rankedlineages_names = "taxid tax_name species genus family order class phylum kingdom domain nan".split()
    rankedlineages = pd.read_csv(
        "/Users/vini/Bio/MGP/databases/taxonomy/rankedlineage.dmp",
        sep="|",
        index_col=0,
        names=rankedlineages_names,
    )
    rankedlineages = rankedlineages.applymap(
        lambda x: x.strip() if isinstance(x, str) else x
    ).dropna(how="all", axis=1)
    tax = tax.join(rankedlineages)

    ranks = rankedlineages_names[2:-1]

    for rank in ranks:
        rank_df = tax.groupby(rank).sum()
        rank_df = pd.concat((rank_df, rank_df / rank_df.sum()), axis=1)
        rank_df.columns = "absolute", "relative"
        rank_df.to_csv(f"{rank}.tsv", sep="\t")


if __name__ == "__main__":
    main()
