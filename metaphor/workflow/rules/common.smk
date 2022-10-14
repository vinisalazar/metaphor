"""
common.smk

    Python functions that act as helpers to run other rules.

Inspired by the same one from the RNASeq workflow:
    https://github.com/snakemake-workflows/rna-seq-star-deseq2/
"""

from glob import glob
from pathlib import Path

import pandas as pd
from snakemake.utils import validate

from metaphor import wrapper_version
from metaphor.config import data_dir


###############################################################
# VALIDATION
# Schema validation to check if config and sample files are OK
###############################################################

validate(config, schema="../schemas/config.schema.yaml")

# Configure the data directory
if config["data_dir"] == "DEFAULT":
    pass
else:
    data_dir = config["data_dir"]

samples = pd.read_csv(
    config["samples"], dtype={"sample_name": str}, sep=None, engine="python"
)
if "unit_name" not in samples.columns:
    samples["unit_name"] = "unit_0"

# Coassembly/cobinning behaviour
# The only difference is really the order of statements: if groups are provided,
# they should be use as 'binning_groups' instead of sample names.
if ("group" not in samples.columns) or (samples["group"].empty):
    samples["group"] = "coassembly" if config["coassembly"] else samples["sample_name"]
    samples["binning_group"] = "cobinning" if config["cobinning"] else samples["group"]
else:
    samples["binning_group"] = "cobinning" if config["cobinning"] else samples["group"]
    samples["group"] = (
        samples["group"] if config["coassembly"] else samples["sample_name"]
    )

samples = samples.fillna("")
samples = samples.set_index(
    ["group", "sample_name", "unit_name"], drop=False
).sort_index()

validate(samples, schema="../schemas/samples.schema.yaml")
group_names = samples["group"].drop_duplicates().to_list()
sample_IDs = samples["sample_name"].drop_duplicates().to_list()
unit_names = samples["unit_name"].drop_duplicates().to_list()
binning_group_names = samples["binning_group"].drop_duplicates().to_list()

if config["host_removal"]["activate"]:
    assert (
        reference_path := Path(config["host_removal"]["reference"])
    ).exists(), f"Host removal reference path '{reference_path}' wasn't found. Please ensure it exists or deactivate host_removal setting."


###############################################################
# TOP LEVEL
# These are top level helpers for all modules
###############################################################


def get_final_output():
    """
    Requires the final output files to be generated at the end of the workflow.

    Consumed by rule 'all'.
    """
    final_output = (
        get_qc_output(),
        get_all_assembly_outputs(),
        get_annotation_output(),
        get_binning_output(),
    )
    return final_output


def is_activated(config_object):
    """
    Checks if an object in the config file is activated or not.
    Used to switch functions/rules on and off.

    config_object: a key in the config YAML file.
    """
    c = config
    for entry in config_object.split("/"):
        c = c.get(entry, {})
    return bool(c.get("activate", False))


def cleanup_rule(config_object, path):
    """
    Checks if a config_object object has the `cleanup` property set as True.
    If it does, deletes the selected `path`.

    config_object: a key in the config YAML file.
    path: the path to be deleted if cleanup is True.
    """
    if config[config_object].get("cleanup", False):
        return f"rm -rf {Path(path)}"
    else:
        return ""


def cleanup_modules():
    """
    Cleans up unnecessary files of the modules selected in config['cleanup_modules'] setting.

    The `modules` dict points to the paths to be deleted.
    """

    modules = {
        "qc": ["output/qc/cutadapt", "output/qc/filtered", "output/qc/merged"],
        "mapping": [
            "output/mapping/bam",
        ],
    }
    paths_to_delete = [
        v
        for k, v in modules.items()
        if k in config["cleanup_modules"]["modules"].split()
    ]
    paths_to_delete = [Path(p) for sublist in paths_to_delete for p in sublist]
    paths_to_delete = [str(p) for p in paths_to_delete if p.exists()]
    return paths_to_delete if is_activated("cleanup_modules") else ()


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

    wrapper: a wrapper from the Snakemake-wrappers repository, e.g. 'cutadapt-pe'
    """
    return str(Path(wrapper_version).joinpath(f"bio/{wrapper}"))


def get_threads_per_task_size(size):
    """
    Determines the number of cores to be used depending on the size of the task.

    The percentage of cores each task uses is set in the config YAML file.
    """
    assert size in (
        choices := ("small", "medium", "big")
    ), f"Size '{size}' must be one of: {choices}."
    threads_ = round(workflow.cores * config[f"cores_per_{size}_task"])
    return 1 if threads_ < 1 else threads_


def get_mb_per_cores(wildcards, threads, task_type="big"):
    """
    Calculates the amount of memory to be used based on the number of total cores
    and the config 'mb_max' object.

    wildcards: Snakemake wildcards (passed on automatically)
    threads: number of threads passed to the workflow
    """
    if config["local_execution"]:
        return threads * int(config["max_mb"] / workflow.cores)
    else:
        return get_max_mb(0.0)


def get_max_mb(margin=0.2):
    """
    Gets the config max_mb and subtracts a margin from itself.
    """
    assert 0 <= margin < 1, f"Margin '{margin}' must be between 0 and 1, zero included."
    mb = config["max_mb"] - (config["max_mb"] * margin)
    return round(mb)


def is_paired_end(sample):
    """
    Checks if a sample is paired-end or not.
    """
    sample_units = (
        samples.loc[sample] if sample in group_names else samples.xs(sample, level=1)
    )
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
    if config["cutadapt"]["activate"]:
        return expand(
            "output/qc/cutadapt/{sample}_{unit}_{read}.fq.gz",
            unit=samples.xs(wildcards.sample, level=1)["unit_name"],
            sample=wildcards.sample,
            read=wildcards.read,
        )
    unit = samples.xs(wildcards.sample, level=1).squeeze()
    fq = "R{}".format(wildcards.read[-1])
    return samples.xs(wildcards.sample, level=1)[fq].tolist().squeeze()


def get_cutadapt_pipe_input(wildcards):
    files = list(
        sorted(
            glob(
                samples.xs(wildcards.sample, level=1)
                .xs(wildcards.unit, level=1)[wildcards.fq]
                .squeeze()
            )
        )
    )
    assert len(files) > 0, "No files were found!"
    return files


def get_cutadapt_input(wildcards):
    unit = samples.xs(wildcards.sample, level=1).xs(wildcards.unit, level=1).squeeze()

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
    unit = (
        samples.xs(wildcards.sample, level=1)
        .xs(wildcards.unit, level=1)[wildcards.read]
        .squeeze()
    )
    return unit


def get_fastqc_input_trimmed(wildcards):
    sample, unit, read = wildcards.sample, wildcards.unit, wildcards.read
    return "output/qc/cutadapt/{sample}_{unit}_{read}.fq.gz"


def get_fastqc_input_merged(wildcards):
    sample, read = wildcards.sample, wildcards.read
    return "output/qc/merged/{sample}_{read}.fq.gz"


def get_fastqc_input_filtered(wildcards):
    sample, read = wildcards.sample, wildcards.read
    return "output/qc/filtered/{sample}_filtered_{read}.fq.gz"


def get_host_removal_input(wildcards):
    if not is_activated("merge_reads"):
        return get_fastqc_input_trimmed
    elif not is_activated("host_removal"):
        return get_fastqc_input_merged
    else:
        return get_fastqc_input_filtered


def get_fastq_groups(wildcards, sense, kind="filtered"):
    if is_activated("host_removal"):
        kind = "filtered"
        add = f"_{kind}_"
    elif is_activated("merge_reads"):
        kind = "merged"
        add = "_"
    elif is_activated("trimming"):
        kind = "cutadapt"
        add = getattr(wildcards, "unit", "_")

    fastq_groups = sorted(
        [
            f"output/qc/{kind}/{sample_name}{add}{sense}.fq.gz"
            for sample_name in samples.loc[wildcards.group, "sample_name"].to_list()
        ]
    )
    return fastq_groups


def get_multiqc_input():
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
    if is_activated("multiqc"):
        return "output/qc/multiqc.html"
    elif is_activated("fastqc"):
        return get_multiqc_input()
    else:
        return ()


###############################################################
# Assembly
###############################################################
def get_contigs_input(expand_=False):
    """Returns coassembly contigs if coassembly is on, else return each sample contig individually"""

    if expand_:
        contigs = expand(
            "output/assembly/megahit/{group}/{group}.contigs.fa",
            group=group_names if config["coassembly"] else sample_IDs,
        )
    else:
        contigs = "output/assembly/megahit/{group}/{group}.contigs.fa"
    return contigs


def get_metaquast_reference(wildcards):
    if config["coassembly"]:
        return config["metaquast"]["coassembly_reference"]
    try:
        reference = samples.loc[wildcards.group, "metaquast_reference"].unique()[0]
        assert Path(reference).exists()
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
            f"Reference file '{reference}' for sample/assembly group '{group}' could not be found."
            raise


def get_group_benchmark_or_log(kind, subworkflow, rule):
    """
    This function was created to prevent formatting errors with snakefmt.

    If such errors are resolved in the future, it may be deprecated.
    """
    if not kind.endswith("s"):
        kind = f"{kind}s"
    base_path = Path(f"output/{kind}/{subworkflow}/{rule}/")
    ext = {"logs": ".log", "benchmarks": ".txt"}
    return str(base_path.joinpath("{group}").with_suffix(ext[kind]))


def get_group_or_sample_file(subworkflow, rule, suffix, add_sample_to_suffix=True):
    """
    Similarly to get_group_benchmark_or_log, this function is
    used to differentiate files generated by coassembly or individual assembly.
    """
    base_path = Path(f"output/{subworkflow}/{rule}/")

    if add_sample_to_suffix:
        suffix = f"{{group}}_{suffix}"

    return str(base_path.joinpath(f"{{group}}/{suffix}"))


def get_metaquast_output():
    if config["coassembly"]:
        if Path(config["metaquast"]["coassembly_reference"]).is_file():
            return "output/assembly/metaquast/coassembly/report.html"
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

ranks = "species genus family order class phylum domain".split()


def get_cog_db_file(filename):
    dirpath = Path(data_dir).joinpath(config["cog_functional_parser"]["db"])
    return str(dirpath.joinpath(filename))


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
    return "output/annotation/diamond/{group}_dmnd.out"


def get_all_diamond_outputs():
    return expand(get_diamond_output(), group=binning_group_names)


def get_concatenate_taxonomies_outputs():
    return expand(
        "output/annotation/cog/tables/concatenated_{rank}_{kind}.tsv",
        rank=ranks
        + [
            "tax",
        ],
        kind=("absolute", "relative"),
    )


functional_kinds = ["categories", "codes", "pathways"]


def get_concatenate_cog_functional_outputs():
    return expand(
        "output/annotation/cog/tables/concatenated_{functional_kinds}_{kind}.tsv",
        functional_kinds=functional_kinds,
        kind=("absolute", "relative"),
    )


def get_lineage_parser_outputs():
    count_types = ("absolute", "relative")
    return (
        get_group_or_sample_file("annotation", "cog", f"{rank}_{count_type}.tsv")
        for rank in ranks
        for count_type in count_types
    )


def get_all_lineage_parser_outputs():
    return expand(get_lineage_parser_outputs(), group=binning_group_names)


def get_prokka_output():
    bins_dict = {}
    for group in binning_group_names:
        bins_dict[group] = glob(
            f"output/binning/DAS_tool/{group}/DAS_tool_DASTool_bins/*"
        )

    for group, list_of_bins in bins_dict.items():
        list_of_bins = [Path(bin_).stem for bin_ in list_of_bins]
        bins_dict[group] = [
            f"output/annotation/prokka/{group}/{bin_}/{bin_}.fna"
            for bin_ in list_of_bins
        ]

    return bins_dict.values()


def get_taxa_plot_outputs():
    return expand(
        "output/annotation/cog/{group}/plots/{group}_{rank}_relative.png",
        rank=ranks,
        group=binning_group_names,
    )


def get_cog_functional_plot_output(group):
    return f"output/annotation/cog/{group}/plots/{group}_cog_categories_relative.png"


def get_all_cog_functional_plot_outputs():
    return expand(
        "output/annotation/cog/{group}/plots/{group}_cog_categories_relative.png",
        group=binning_group_names,
    )


def get_annotation_output():
    annotations = {
        "diamond": [
            get_all_diamond_outputs(),
            add_data_dir(config["lineage_parser"]["names"]),
            add_data_dir(config["lineage_parser"]["nodes"]),
        ],
        "taxonomy_parser": (get_concatenate_taxonomies_outputs()),
        "plot_taxonomies": (get_taxa_plot_outputs()),
        "cog_functional_parser": (
            get_concatenate_cog_functional_outputs(),
            get_database_outputs(),
        ),
        "lineage_parser": (
            get_all_lineage_parser_outputs(),
            config["lineage_parser"]["rankedlineage"],
        ),
        "plot_cog": get_all_cog_functional_plot_outputs(),
        "prokka": get_prokka_output()
        if is_activated("das_tool")
        else (),  # Can't run Prokka without DAS Tool!
    }

    needs_activation = (
        "cog_functional_parser",
        "lineage_parser",
        "plot_cog_functional",
        "plot_taxonomies",
        "prokka",
    )
    annotation_output = []

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
        if config["cutadapt"]["activate"]:
            return expand(
                "output/qc/cutadapt/{sample}_{unit}_R1.fq.gz",
                unit=samples.xs(wildcards.sample, level=1)["unit_name"].squeeze(),
                sample=wildcards.sample,
            )
        unit = samples.xs(wildcards.sample, level=1).squeeze()
        if all(pd.isna(unit["R1"])):
            # SRA sample (always paired-end for now)
            accession = unit["sra"]
            return expand("sra/{accession}_R1.fq", accession=accession)
        sample_units = samples.xs(wildcards.sample, level=1).squeeze()
        return sample_units["R1"]
    if is_paired_end(wildcards.sample):
        return "output/qc/merged/{sample}_R1.fq.gz"
    return "output/qc/merged/{sample}_single.fq.gz"


def get_map_reads_input_R2(wildcards):
    if is_paired_end(wildcards.sample):
        if not is_activated("merge_reads"):
            if config["cutadapt"]["activate"]:
                return expand(
                    "output/qc/cutadapt/{sample}_{unit}_R1.fq.gz",
                    unit=samples.xs(wildcards.sample, level=1)["unit_name"].squeeze(),
                    sample=wildcards.sample,
                )
            unit = samples.xs(wildcards.sample, level=1).squeeze()
            if all(pd.isna(unit["R1"])):
                # SRA sample (always paired-end for now)
                accession = unit["sra"]
                return expand("sra/{accession}_R2.fq", accession=accession)
            sample_units = samples.xs(wildcards.sample, level=1).squeeze()
            return sample_units["R2"]
        return ("output/qc/merged/{sample}_R2.fq.gz",)
    return ""


def get_mapping_output():
    return expand(
        "output/mapping/bam/{binning_group}/{group}.{kind}",
        binning_group=binning_group_names,
        group=group_names,
        kind=("map.bam", "sorted.bam", "flagstat.txt"),
    )


###############################################################
# Binning
###############################################################

binners = [b for b in ("concoct", "metabat2", "vamb") if is_activated(b)]


def get_DAS_tool_input():
    scaffolds2bin = (
        lambda binner: f"output/binning/{binner}/{{binning_group}}/{binner}_scaffolds2bin.tsv"
    )
    return sorted(scaffolds2bin(b) for b in binners)


def get_fasta_bins():
    binners = {
        "metabat2": "output/binning/metabat2/*/*.fa",
        "concoct": "output/binning/concoct/*/fasta_bins/*.fa",
        "vamb": "output/binning/vamb/*/bins/*.fna",
    }

    bins = sorted(glob(v) for k, v in binners.items() if is_activated(k))
    return bins


def get_vamb_output():
    return "output/binning/vamb/{binning_group}/clusters.tsv"


def get_all_vamb_output():
    return tuple(
        expand(
            "output/binning/vamb/{binning_group}/clusters.tsv",
            binning_group=binning_group_names,
        )
    )


def get_binning_report_output(binning_group):
    plots = [f"bin_{i}" for i in "quality scores quantity sizes N50".split()]
    plots_dict = {
        plot: report(
            f"output/binning/plots/{binning_group}/{plot}.png", category="Binning"
        )
        for plot in plots
    }
    return plots_dict


def get_binning_output():
    binners = {
        "vamb": get_all_vamb_output(),
        "metabat2": expand(
            "output/binning/metabat2/{binning_group}", binning_group=binning_group_names
        ),
        "concoct": expand(
            "output/binning/concoct/{binning_group}", binning_group=binning_group_names
        ),
        "das_tool": [
            expand(
                "output/binning/DAS_tool/{binning_group}/{binning_group}_DASTool_summary.tsv",
                binning_group=binning_group_names,
            ),
            [
                get_binning_report_output(binning_group).values()
                for binning_group in binning_group_names
            ]
            if (config["das_tool"]["bins_report"]) and (is_activated("das_tool"))
            else [],
        ],
    }
    return (v for k, v in binners.items() if is_activated(k))


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
            "output/postprocessing/memory_barplot_errorbar.png",
            "output/postprocessing/cleanup.txt"
            if is_activated("cleanup_modules")
            else (),
        )
    else:
        return ()


def get_processing_benchmarks():
    return "output/benchmarks/postprocessing/processing_benchmarks.csv"


def add_data_dir(path) -> str:
    """Add the 'data_dir' config variable as a prefix to path."""
    return str(Path(data_dir).joinpath(path))
