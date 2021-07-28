"""
Common rules: these are special Python help functions to run other rules.

Inspired by the same one from the RNASeq workflow:
    https://github.com/snakemake-workflows/rna-seq-star-deseq2/
"""

import glob
from pathlib import Path

import pandas as pd
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], dtype={"sample_name": str})
samples = samples.set_index("sample_name", drop=False).sort_index()

validate(samples, schema="../schemas/samples.schema.yaml")

sample_IDs = samples.index.to_list()
fq1 = samples["fq1"].to_list()
fq2 = samples["fq2"].to_list()


def get_final_output():

    assembly_output = _get_assembly_output()
    mapping_output = _get_mapping_output()
    annotation_output = _get_annotation_output()

    return (
        "output/preprocess/multiqc.html",
        assembly_output,
        mapping_output,
        annotation_output,
    )


def _get_assembly_output():
    return expand(
        "output/assembly/megahit/{sample}/{sample}.contigs.fa", sample=sample_IDs
    )


def _get_mapping_output():
    return expand(
        "output/mapping/bam/{sample}.{kind}",
        sample=sample_IDs,
        kind=("map.bam", "sorted.bam"),
    )


def _get_binning_output():
    return directory("output/binning/vamb/")


def _get_annotation_output():
    return expand("output/annotation/diamond/{sample}.xml", sample=sample_IDs)


def pathfinder(path: str):
    return str(Path(__file__).joinpath(Path(path)).absolute())
