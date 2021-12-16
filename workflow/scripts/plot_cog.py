#!/usr/bin/env python
import argparse
import logging
from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt


#######################################
# Category heatmap
#######################################


def create_heatmap(categories_file, filter_categories=True, cutoff=0.01):
    dataframe = pd.read_csv(categories_file, sep="\t", index_col=0)

    if filter_categories:
        dataframe = dataframe.drop("Function unknown")
        dataframe = dataframe.drop("General function prediction only")
        dataframe = dataframe[dataframe > 0.01]
        dataframe = dataframe.dropna()
        filtered = (abs(dataframe.sum() - 1)).mean()
        dataframe = dataframe / dataframe.sum()
        vmax, vmin = None, None
    else:
        nlargest = categories.sum(axis=1).nlargest().index.to_list()
        for ix in nlargest:
            if ix not in ("Function unknown", "General function prediction only"):
                vmax = dataframe.loc[ix].max()
                break

        vmin = cutoff

    fig, ax = plt.subplots(figsize=(3 + len(dataframe.columns), 6))
    sns.heatmap(dataframe, cmap="viridis", vmax=vmax, vmin=vmin, ax=ax)


#######################################
# Taxa bar plots
#######################################


def calculate_legend_width(index):
    max_len = max(len(str(i)) for i in index)
    if max_len > 20:
        return (3, 1)
    elif max_len < 10:
        return (8, 1)
    else:
        return (4, 1)


def format_index(index):
    new_index = []
    for s in index:
        try:
            s = str(s)
            if len(s.split()) >= 2 and "Undetermined/other" not in s:
                s = s.split()[0][0] + ". " + " ".join(s.split()[1:])
                s = "$" + s + "$"
        except (TypeError, ValueError):
            pass
        new_index.append(s)

    return pd.Index(new_index, name=index.name)


def create_tax_barplot(dataframe):
    figsize = (10, len(dataframe.columns))
    dataframe.index = format_index(dataframe.index)
    fig, axs = plt.subplots(
        nrows=1,
        ncols=2,
        figsize=figsize,
        gridspec_kw={"width_ratios": calculate_legend_width(dataframe.index)},
    )
    axsdataframe.T.plot(kind="barh", stacked=True, ax=axs[0], cmap="tab20c")
    handles, labels = axs[0].get_legend_handles_labels()
    axs[0].get_legend().remove()
    axs[0].set_xlabel("Relative abundance")
    axs[1].legend(handles, labels, title=dataframe.index.name.capitalize())
    axs[1].set_xticks([])
    axs[1].set_yticks([])
    axs[1].axis("off")
