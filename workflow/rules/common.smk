"""
Common rules: these are special Python help functions to run other rules.

Inspired by the same one from the RNASeq workflow:
    https://github.com/snakemake-workflows/rna-seq-star-deseq2/

They are divided into three parts:
    - helpers: utility functions
    - input functions: take wildcards as input from Snakemake
    - output functions: do not take wildcards as inputs
"""

from glob import glob
from pathlib import Path

import pandas as pd
from snakemake.utils import validate


###############################################################
# VALIDATION
# Schema validation to check if config and sample files are OK
###############################################################

validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], dtype={"sample_name": str, "unit_name": str})
samples = samples.fillna("")
samples = samples.set_index(["sample_name", "unit_name"], drop=False).sort_index()

validate(samples, schema="../schemas/samples.schema.yaml")
sample_IDs = samples["sample_name"].drop_duplicates().to_list()
unit_names = samples["unit_name"].drop_duplicates().to_list()


###############################################################
# TOP LEVEL
# These are top level helpers for all modules
###############################################################


def get_final_output():
    final_output = (
        get_qc_output(),
        get_all_assembly_outputs(),
        get_mapping_output(),
        get_annotation_output(),
        get_binning_output(),
    )
    return final_output


def is_activated(xpath):
    c = config
    for entry in xpath.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


def get_parent(path: str) -> str:
    """Returns parent of path in string form."""
    return str(Path(path).parent)


def allow():
    """
    Allows the workflow to run even when failing certain assertions.

    Must be explicitly declared.
    """
    import sys

    allow_on = ("--lint", "--dry-run", "--dryrun", "-n")
    if any(term in sys.argv for term in allow_on):
        return True


def get_wrapper(wrapper):
    """
    Builds the string for the 'wrapper' directive
    based on 'wrapper_version' key in the config.
    """
    return str(Path(config["wrapper_version"]).joinpath(f"bio/{wrapper}"))


def get_mem_mb(wildcards, threads):
    return threads * config["resources"]["mb_per_thread"]


def is_paired_end(sample):
    sample_units = samples.loc[sample]
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


###############################################################
# QC
###############################################################


def get_fastqs(wildcards):
    if config["trimming"]["activate"]:
        return expand(
            "output/qc/cutadapt/{sample}_{unit}_{read}.fq.gz",
            unit=samples.loc[wildcards.sample, "unit_name"],
            sample=wildcards.sample,
            read=wildcards.read,
        )
    unit = samples.loc[wildcards.sample]
    fq = "R{}".format(wildcards.read[-1])
    return samples.loc[wildcards.sample, fq].tolist()


def get_cutadapt_pipe_input(wildcards):
    files = list(
        sorted(glob(samples.loc[wildcards.sample].loc[wildcards.unit, wildcards.fq]))
    )
    assert len(files) > 0, "No files were found!"
    return files


def get_cutadapt_input(wildcards):
    unit = samples.loc[wildcards.sample].loc[wildcards.unit]

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
    unit = samples.loc[wildcards.sample].loc[wildcards.unit][wildcards.read]
    return unit


def get_fastqc_input_trimmed(wildcards):
    sample, unit, read = wildcards.sample, wildcards.unit, wildcards.read
    return "output/qc/cutadapt/{sample}_{unit}_{read}.fq.gz"


def get_fastqc_input_merged(wildcards):
    sample, read = wildcards.sample, wildcards.read
    return "output/qc/merged/{sample}_{read}.fq.gz"


def get_multiqc_input(wildcards):
    raw = expand(
        "output/qc/fastqc/{sample}-{unit}-{read}-raw_fastqc.zip",
        sample=sample_IDs,
        unit=unit_names,
        read=["R1", "R2"],
    )
    trimmed = expand(
        "output/qc/fastqc/{sample}-{unit}-{read}-trimmed_fastqc.zip",
        sample=sample_IDs,
        unit=unit_names,
        read=["R1", "R2"],
    )
    merged = expand(
        "output/qc/fastqc/{sample}-{read}-merged_fastqc.zip",
        sample=sample_IDs,
        read=["R1", "R2"],
    )

    return raw + trimmed + merged


def get_qc_output():
    return "output/qc/multiqc.html"


###############################################################
# Assembly
###############################################################


def get_contigs_input(expand_=False):
    """Returns coassembly contigs if coassembly is on, else return each sample contig individually"""
    if config["coassembly"]:
        contigs = "output/assembly/megahit/coassembly/coassembly.contigs.fa"
    else:
        if expand_:
            contigs = expand(
                "output/assembly/megahit/{sample}/{sample}.contigs.fa",
                sample=sample_IDs,
            )
        else:
            contigs = "output/assembly/megahit/{sample}/{sample}.contigs.fa"
    return contigs


def get_metaquast_reference(wildcards):
    sample = wildcards.sample
    try:
        if config["coassembly"]:
            return config["metaquast"]["coassembly_reference"]
        else:
            reference = samples.loc[sample, "metaquast_reference"].unique()[0]
            assert Path(reference).is_file()
            return reference
    except (KeyError, IndexError):
        if allow():
            return ()
        else:
            print(
                f"Column 'metaquast_reference' could not be found or does not contain any files."
            )
            raise
    except AssertionError:
        if allow():
            return ()
        else:
            f"Reference file '{reference}' for sample '{sample}' could not be found."
            raise


def get_coassembly_benchmark_or_log(kind, subworkflow, rule):
    """
    This function was created to prevent formatting errors with snakefmt.

    If such errors are resolved in the future, it may be deprecated.
    """
    if not kind.endswith("s"):
        kind = f"{kind}s"
    base_path = Path(f"output/{kind}/{subworkflow}/{rule}/")
    ext = {"logs": ".log", "benchmarks": ".txt"}
    if config["coassembly"]:
        return str(base_path.joinpath("coassembly").with_suffix(ext[kind]))
    else:
        return str(base_path.joinpath("{sample}").with_suffix(ext[kind]))


def get_coassembly_or_sample_file(subworkflow, rule, suffix, add_sample_to_suffix=True):
    """
    Similarly to get_coassembly_benchmark_or_log, this function is
    used to differentiate files generated by coassembly or individual assembly.
    """
    base_path = Path(f"output/{subworkflow}/{rule}/")

    if config["coassembly"]:
        return str(base_path.joinpath(f"coassembly/coassembly_{suffix}"))
    else:
        if add_sample_to_suffix:
            suffix = f"{{sample}}_{suffix}"
        return str(base_path.joinpath(f"{{sample}}/{suffix}"))


def cleanup_megahit(contigs):
    if config["megahit"]["cleanup"]:
        intermediate_contigs = str(
            Path(contigs).parent.joinpath("intermediate_contigs")
        )
        return f"rm -rf {intermediate_contigs}"
    else:
        return ""


def get_metaquast_output():
    if config["coassembly"]:
        if Path(config["metaquast"]["coassembly_reference"]).is_file():
            return "output/assembly/metaquast_coassembly/report.html"
        else:
            return ()
    else:
        return expand(
            "output/assembly/metaquast/{sample}/report.html", sample=sample_IDs
        )


# def metaquast_cleanup(report):
#     if config["metaquast"]["cleanup"]:


def get_assembly_report(plot=None):
    base_dir = Path("output/assembly/assembly_report/")
    if not plot:
        return str(base_dir.joinpath("assembly_report.tsv"))
    else:
        return str(base_dir.joinpath(f"{plot}.png"))


def get_all_assembly_outputs():
    assemblies = [
        get_contigs_input(expand_=True),
    ]
    assemblies.append(get_assembly_report())
    if is_activated("metaquast"):
        assemblies.append(get_metaquast_output())

    return assemblies


###############################################################
# Annotation
###############################################################

ranks = "species genus family order class phylum kingdom domain".split()


def get_cog_db_file(filename):
    return str(Path(config["cog_parser"]["db"]).joinpath(filename))


def get_database_outputs():
    db_outputs = [
        get_cog_db_file("cog-20.fa.gz"),
        get_cog_db_file("cog-20.cog.csv"),
        get_cog_db_file("cog-20.org.csv"),
        get_cog_db_file("cog-20.def.tab"),
        get_cog_db_file("fun-20.tab"),
    ]

    return db_outputs


def get_diamond_output():
    return (
        "output/annotation/diamond/{sample}_dmnd.out"
        if not config["coassembly"]
        else "output/annotation/diamond/coassembly_dmnd.out"
    )


def get_all_diamond_outputs():
    return expand(get_diamond_output(), sample=sample_IDs)


def get_all_cog_parser_outputs():
    cog_outputs = ["categories", "codes", "tax", "pathways"]
    cog_valid_output_kinds = cog_outputs + ranks
    return (
        expand(
            "output/annotation/cog/{sample}/{sample}_{kind}.tsv",
            sample=sample_IDs,
            kind=cog_valid_output_kinds,
        )
        if not config["coassembly"]
        else expand(
            "output/annotation/cog/coassembly/coassembly_{kind}.tsv",
            sample=sample_IDs,
            kind=cog_valid_output_kinds,
        )
    )


def get_concatenate_cog_outputs():
    return (
        (
            "output/annotation/cog/tables/COG_categories_absolute.tsv",
            "output/annotation/cog/tables/COG_categories_relative.tsv",
            "output/annotation/cog/tables/COG_codes_absolute.tsv",
            "output/annotation/cog/tables/COG_codes_relative.tsv",
            "output/annotation/cog/tables/COG_taxs_absolute.tsv",
            "output/annotation/cog/tables/COG_taxs_relative.tsv",
            "output/annotation/cog/tables/COG_pathways_absolute.tsv",
            "output/annotation/cog/tables/COG_pathways_relative.tsv",
            "output/annotation/cog/tables/COG_species_absolute.tsv",
            "output/annotation/cog/tables/COG_species_relative.tsv",
            "output/annotation/cog/tables/COG_genus_absolute.tsv",
            "output/annotation/cog/tables/COG_genus_relative.tsv",
            "output/annotation/cog/tables/COG_family_absolute.tsv",
            "output/annotation/cog/tables/COG_family_relative.tsv",
            "output/annotation/cog/tables/COG_order_absolute.tsv",
            "output/annotation/cog/tables/COG_order_relative.tsv",
            "output/annotation/cog/tables/COG_class_absolute.tsv",
            "output/annotation/cog/tables/COG_class_relative.tsv",
            "output/annotation/cog/tables/COG_phylum_absolute.tsv",
            "output/annotation/cog/tables/COG_phylum_relative.tsv",
            "output/annotation/cog/tables/COG_kingdom_absolute.tsv",
            "output/annotation/cog/tables/COG_kingdom_relative.tsv",
            "output/annotation/cog/tables/COG_domain_absolute.tsv",
            "output/annotation/cog/tables/COG_domain_relative.tsv",
        )
        if not config["coassembly"]
        else ()
    )


def get_lineage_parser_outputs():
    return (
        get_coassembly_or_sample_file("annotation", "cog", f"{rank}.tsv")
        for rank in ranks
    )


def get_all_lineage_parser_outputs():
    return expand(get_lineage_parser_outputs(), sample=sample_IDs)


def get_prokka_output():
    return expand("output/annotation/prokka/{sample}/{sample}.faa", sample=sample_IDs)


def get_taxa_plot_outputs():
    return (
        "output/annotation/cog/plots/COG_species_relative.png",
        "output/annotation/cog/plots/COG_genus_relative.png",
        "output/annotation/cog/plots/COG_family_relative.png",
        "output/annotation/cog/plots/COG_order_relative.png",
        "output/annotation/cog/plots/COG_class_relative.png",
        "output/annotation/cog/plots/COG_phylum_relative.png",
        "output/annotation/cog/plots/COG_domain_relative.png",
    )


def get_annotation_output():
    annotations = {
        "diamond": [
            get_all_diamond_outputs(),
            config["lineage_parser"]["names"],
            config["lineage_parser"]["nodes"],
        ],
        "cog_parser": (
            get_all_cog_parser_outputs(),
            get_concatenate_cog_outputs(),
            get_database_outputs(),
        ),
        "lineage_parser": (
            get_all_lineage_parser_outputs(),
            config["lineage_parser"]["rankedlineage"],
        ),
        "plot_cog": get_taxa_plot_outputs(),
        "prokka": get_prokka_output(),
    }

    needs_activation = ("cog_parser", "lineage_parser", "plot_cog", "prokka")
    annotation_output = []

    if not is_activated("cog_parser"):
        config["lineage_parser"]["activate"] = False
        config["plot_cog"]["activate"] = False

    for k, v in annotations.items():
        if not is_activated(k) and k in needs_activation:
            continue
        else:
            annotation_output.append(v)

    return annotation_output


###############################################################
# Mapping
###############################################################
def get_map_reads_input_R1(wildcards):
    if not is_activated("merge_reads"):
        if config["trimming"]["activate"]:
            return expand(
                "output/qc/cutadapt/{sample}_{unit}_R1.fq.gz",
                unit=samples.loc[wildcards.sample, "unit_name"],
                sample=wildcards.sample,
            )
        unit = samples.loc[wildcards.sample]
        if all(pd.isna(unit["R1"])):
            # SRA sample (always paired-end for now)
            accession = unit["sra"]
            return expand("sra/{accession}_R1.fq", accession=accession)
        sample_units = samples.loc[wildcards.sample]
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
                    unit=samples.loc[wildcards.sample, "unit_name"],
                    sample=wildcards.sample,
                )
            unit = samples.loc[wildcards.sample]
            if all(pd.isna(unit["R1"])):
                # SRA sample (always paired-end for now)
                accession = unit["sra"]
                return expand("sra/{accession}_R2.fq", accession=accession)
            sample_units = samples.loc[wildcards.sample]
            return sample_units["R2"]
        return ("output/qc/merged/{sample}_R2.fq.gz",)
    return ""


def get_mapping_output():
    return expand(
        "output/mapping/bam/{sample}.{kind}",
        sample=sample_IDs,
        kind=("map.bam", "sorted.bam", "flagstat.txt"),
    )


###############################################################
# Binning
###############################################################

binners = [b for b in ("concoct", "metabat2", "vamb") if is_activated(b)]


def get_DAS_tool_input():
    scaffolds2bin = lambda binner: f"output/binning/DAS_tool/{binner}_scaffolds2bin.tsv"
    return sorted(scaffolds2bin(b) for b in binners)


def get_fasta_bins():
    binners = {
        "metabat2": "output/binning/metabat2/*.fa",
        "concoct": "output/binning/concoct/fasta_bins/*.fa",
        "vamb": "output/binning/vamb/bins/*.fna",
    }

    bins = sorted(glob(v) for k, v in binners.items() if is_activated(k))
    return bins


def get_vamb_output():
    return ("output/binning/vamb/clusters.tsv", "output/binning/vamb/log.txt")


def get_binning_output():
    binners = {
        "vamb": get_vamb_output()[0],
        "metabat2": "output/binning/metabat2/",
        "concoct": "output/binning/concoct/",
        "das_tool": "output/binning/DAS_tool/DAS_tool_proteins.faa",
    }
    return sorted(v for k, v in binners.items() if is_activated(k))


###############################################################
# Postprocessing
###############################################################
def get_postprocessing_output():
    if is_activated("postprocessing"):
        return (
            get_processing_benchmarks(),
            "output/postprocessing/runtime_barplot_sum.png",
            "output/postprocessing/runtime_barplot_errorbar.png",
            "output/postprocessing/memory_barplot_sum.png",
            "output/postprocessing/memory_barplot_errorbar.png",
        )
    else:
        return ()


def get_processing_benchmarks():
    return "output/benchmarks/postprocessing/processing_benchmarks.csv"
