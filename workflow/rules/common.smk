"""
Common rules: these are special Python help functions to run other rules.

Inspired by the same one from the RNASeq workflow:
    https://github.com/snakemake-workflows/rna-seq-star-deseq2/
"""

from glob import glob
from pathlib import Path

import pandas as pd
from snakemake.utils import validate

validate(config, schema="../schemas/config.schema.yaml")

samples = (
    pd.read_csv(config["samples"], dtype={"sample_name": str})
    .set_index("sample_name", drop=False)
    .sort_index()
)

validate(samples, schema="../schemas/samples.schema.yaml")
sample_IDs = samples["sample_name"].drop_duplicates().to_list()

units = (
    pd.read_csv(config["units"], dtype={"sample_name": str, "unit_name": str})
    .set_index(["sample_name", "unit_name"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")
unit_names = units["unit_name"].drop_duplicates().to_list()


# Helpers
def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


def is_paired_end(sample):
    sample_units = units.loc[sample]
    fq2_null = sample_units["R2"].isnull()
    paired = ~fq2_null
    all_paired = paired.all()
    all_single = (~paired).all()
    assert (
        all_single or all_paired
    ), "invalid units for sample {}, must be all paired end or all single end".format(
        sample
    )
    return all_paired


# Inputs
def get_fastqs(wildcards):
    if config["trimming"]["activate"]:
        return expand(
            "output/qc/cutadapt/{sample}_{unit}_{read}.fq.gz",
            unit=units.loc[wildcards.sample, "unit_name"],
            sample=wildcards.sample,
            read=wildcards.read,
        )
    unit = units.loc[wildcards.sample]
    fq = "R{}".format(wildcards.read[-1])
    return units.loc[wildcards.sample, fq].tolist()


def get_cutadapt_pipe_input(wildcards):
    files = list(
        sorted(glob(units.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq]))
    )
    assert len(files) > 0, "No files were found!"
    return files


def get_cutadapt_input(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit]

    if unit["R1"].endswith("gz"):
        ending = ".gz"
    else:
        ending = ""

    if pd.isna(unit["R2"]):
        # single end local sample
        return "pipe/qc/cutadapt/{sample}_{unit}.fq{ending}".format(
            sample=unit.sample_name, unit=unit.unit_name, ending=ending
        )
    else:
        # paired end local sample
        return expand(
            "pipe/qc/cutadapt/{sample}_{unit}_{{read}}.fq{ending}".format(
                sample=unit.sample_name, unit=unit.unit_name, ending=ending
            ),
            read=["R1", "R2"],
        )


def get_fastqc_input_raw(wildcards):
    unit = units.loc[wildcards.sample].loc[wildcards.unit][wildcards.read]
    return unit


def get_fastqc_input_merged(wildcards):
    sample, read = wildcards.sample, wildcards.read
    return "output/qc/merged/{sample}_{read}.fq.gz"



def get_multiqc_input(wildcards):
    raw = expand("output/qc/fastqc/{sample}-{unit}-{read}_fastqc.zip", sample=sample_IDs, unit=unit_names, read=["R1", "R2"])
    merged = expand("output/qc/fastqc/{sample}-{read}_fastqc.zip", sample=sample_IDs, read=["R1", "R2"])
    return raw + merged


def get_map_reads_input_R1(wildcards):
    if not is_activated("merge_reads"):
        if config["trimming"]["activate"]:
            return expand(
                "output/qc/cutadapt/{sample}_{unit}_R1.fq.gz",
                unit=units.loc[wildcards.sample, "unit_name"],
                sample=wildcards.sample,
            )
        unit = units.loc[wildcards.sample]
        if all(pd.isna(unit["R1"])):
            # SRA sample (always paired-end for now)
            accession = unit["sra"]
            return expand("sra/{accession}_R1.fq", accession=accession)
        sample_units = units.loc[wildcards.sample]
        return sample_units["R1"]
    if is_paired_end(wildcards.sample):
        return "output/qc/merged/{sample}_R1.fq.gz"
    return "output/qc/merged/{sample}_single.fq.gz"


def get_map_reads_input_R2(wildcards):
    if is_paired_end(wildcards.sample):
        if not is_activated("merge_reads"):
            if config["trimming"]["activate"]:
                return expand(
                    "output/qc/cutadapt/{sample}_{unit}_R1.fq.gz",
                    unit=units.loc[wildcards.sample, "unit_name"],
                    sample=wildcards.sample,
                )
            unit = units.loc[wildcards.sample]
            if all(pd.isna(unit["R1"])):
                # SRA sample (always paired-end for now)
                accession = unit["sra"]
                return expand("sra/{accession}_R2.fq", accession=accession)
            sample_units = units.loc[wildcards.sample]
            return sample_units["R2"]
        return ("output/qc/merged/{sample}_R2.fq.gz",)
    return ""


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
