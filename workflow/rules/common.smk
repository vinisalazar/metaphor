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

# Fastq files
fq1 = samples["fq1"].to_list()
fq2 = samples["fq2"].to_list()
fq_clean = samples["clean"] = samples["sample_name"].apply(lambda s: f"output/qc/interleave/{s}-clean.fq").to_list()
fqs = fq1 + fq2 + fq_clean
fq_basenames = [str(Path(i).stem) for i in fqs]


def get_multiqc_input():
    return expand(
            "output/qc/fastqc/{sample}_fastqc.zip",
            sample=fq_basenames
        )


# Outputs
def get_final_output():

    return (
        get_qc_output(),
        get_assembly_output(),
        get_mapping_output(),
        get_annotation_output(),
    )


def get_qc_output():
    return "output/qc/multiqc.html"


def get_assembly_output():
    return expand(
        "output/assembly/megahit/{sample}/{sample}.contigs.fa", sample=sample_IDs
    )


def get_mapping_output():
    return expand(
        "output/mapping/bam/{sample}.{kind}",
        sample=sample_IDs,
        kind=("map.bam", "sorted.bam"),
    )


def get_binning_output():
    return directory("output/binning/vamb/")


def get_annotation_output():
    return expand("output/annotation/diamond/{sample}.xml", sample=sample_IDs)
