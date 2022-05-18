"""
qc.smk

    Quality Control module. Filter reads by quality and trim with CutAdapt. Run FastQC at each stage
    and join results in MultiQC report.

QC rules:
    - cutadapt_pipe: get pipe input for cutadapt
    - cutadapt_pe: trim paired end reads with cutadapt
    - merge_fastqs: merge files from different lanes in the same sample with cat
    - fastqc_raw: check quality of raw reads with FastQC
    - fastqc_trimmed: check quality of trimmed reads with FastQC
    - fastqc_merged: check quality of trimmed and merged reads with FastQC
    - multiqc: combine reports of rules 'fastqc_raw', '_trimmed' and '_merged' with MultiQC
"""
from pathlib import Path


rule cutadapt_pipe:
    """Pipe reads into cutadapt."""
    input:
        get_cutadapt_pipe_input,
    output:
        pipe("pipe/qc/cutadapt/{sample}_{unit}_{fq}.{ext}"),
    resources:
        mem_mb=get_max_mb(0.5),
    log:
        "output/logs/qc/cutadapt/{sample}_{unit}_{fq}_pipe.{ext}.log",
    benchmark:
        "output/benchmarks/qc/cutadapt/{sample}_{unit}_{fq}_pipe.{ext}.txt"
    wildcard_constraints:
        sample="|".join(sample_IDs),
        unit="|".join(unit_names),
        ext=r"fq|fq\.gz|fastq|fastq\.gz",
    threads: 1
    conda:
        "../envs/utils.yaml"
    shell:
        """cat {input} > {output} 2> {log}"""


rule cutadapt_pe:
    """
    Trim paired end reads with cutadapt.
    """
    input:
        get_cutadapt_input,
    output:
        fastq1="output/qc/cutadapt/{sample}_{unit}_R1.fq.gz",
        fastq2="output/qc/cutadapt/{sample}_{unit}_R2.fq.gz",
        qc="output/qc/cutadapt/{sample}_{unit}.paired.qc.txt",
    resources:
        mem_mb=get_max_mb(0.5),
    log:
        "output/logs/qc/cutadapt/{sample}-{unit}.log",
    benchmark:
        "output/benchmarks/qc/cutadapt/{sample}-{unit}.txt"
    threads: round(workflow.cores * config["cores_per_small_task"])
    params:
        # adapters=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
        others="",
        adapters="-a AGAGCACACGTCTGAACTCCAGTCAC -g AGATCGGAAGAGCACACGT -A AGAGCACACGTCTGAACTCCAGTCAC -G AGATCGGAAGAGCACACGT",
        extra=(
            f"--minimum-length {config['cutadapt']['minimum_length']} "
            + f"--quality-cutoff {config['cutadapt']['quality_cutoff']} "
            + f"--quality-base {config['cutadapt']['phred']} "
            + f"-u {config['cutadapt']['clip_r5']} "
            + f"-u -{config['cutadapt']['clip_r3']} "
            + f"-U {config['cutadapt']['clip_r5']} "
            + f"-U -{config['cutadapt']['clip_r3']} "
        ),
    wildcard_constraints:
        sample="|".join(sample_IDs),
    wrapper:
        get_wrapper("cutadapt/pe")


rule merge_fastqs:
    """
    Concatenate paired-end reads from different units.
    """
    input:
        get_fastqs,
    output:
        "output/qc/merged/{sample}_{read}.fq.gz",
    threads: 1
    resources:
        mem_mb=get_max_mb(0.5),
    log:
        "output/logs/qc/merge_fastqs/{sample}.{read}.log",
    benchmark:
        "output/benchmarks/qc/merge_fastqs/{sample}.{read}.txt"
    wildcard_constraints:
        sample="|".join(sample_IDs),
        read="single|R1|R2",
    conda:
        "../envs/utils.yaml"
    shell:
        """cat {input} > {output} 2> {log}"""


rule fastqc_raw:  # qc on raw, unmerged reads
    input:
        get_fastqc_input_raw,
    output:
        zip="output/qc/fastqc/{sample}-{unit}-{read}-raw_fastqc.zip",
        html="output/qc/fastqc/{sample}-{unit}-{read}-raw.html",
    params:
        "--quiet",
    log:
        "output/logs/qc/fastqc_raw/{sample}-{unit}-{read}.log",
    benchmark:
        "output/benchmarks/qc/fastqc_raw/{sample}-{unit}-{read}.txt"
    threads: round(workflow.cores * config["cores_per_medium_task"])
    wrapper:
        get_wrapper("fastqc")


rule fastqc_trimmed:  # qc on trimmed reads
    input:
        get_fastqc_input_trimmed,
    output:
        zip="output/qc/fastqc/{sample}-{unit}-{read}-trimmed_fastqc.zip",
        html="output/qc/fastqc/{sample}-{unit}-{read}-trimmed.html",
    params:
        "--quiet",
    log:
        "output/logs/qc/fastqc_trimmed/{sample}-{unit}-{read}.log",
    benchmark:
        "output/benchmarks/qc/fastqc_trimmed/{sample}-{unit}-{read}.txt"
    threads: round(workflow.cores * config["cores_per_medium_task"])
    wrapper:
        get_wrapper("fastqc")


rule fastqc_merged:  # qc on trimmed, merged reads
    input:
        get_fastqc_input_merged,
    output:
        zip="output/qc/fastqc/{sample}-{read}-merged_fastqc.zip",
        html="output/qc/fastqc/{sample}-{read}-merged.html",
    params:
        "--quiet",
    log:
        "output/logs/qc/fastqc_merged/{sample}-{read}.log",
    benchmark:
        "output/benchmarks/qc/fastqc_merged/{sample}-{read}.txt"
    threads: round(workflow.cores * config["cores_per_medium_task"])
    wrapper:
        get_wrapper("fastqc")


rule multiqc:
    input:
        get_multiqc_input(),
    output:
        report="output/qc/multiqc.html",
    log:
        "output/logs/qc/multiqc.log",
    benchmark:
        "output/benchmarks/qc/multiqc.txt"
    wrapper:
        get_wrapper("multiqc")
