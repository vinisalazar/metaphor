"""
qc.smk

    Quality Control module. Filter reads by quality and trim with fastp. Run FastQC at each stage
    and join results in MultiQC report.

QC rules:
    - fastp_pipe: get pipe input for fastp
    - fastp_pe: trim paired end reads with fastp
    - merge_fastqs: merge files from different lanes in the same sample with cat
    - fastqc_raw: check quality of raw reads with FastQC
    - fastqc_trimmed: check quality of trimmed reads with FastQC
    - fastqc_merged: check quality of trimmed and merged reads with FastQC
    - multiqc: combine reports of rules 'fastqc_raw', '_trimmed' and '_merged' with MultiQC
"""
from pathlib import Path


rule fastp_pipe:
    """Pipe reads into fastp."""
    input:
        get_fastp_pipe_input,
    output:
        pipe("pipe/qc/fastp/{sample}_{unit}_{fq}.{ext}"),
    resources:
        mem_mb=get_max_mb(0.5),
    log:
        "output/logs/qc/fastp/{sample}_{unit}_{fq}_pipe.{ext}.log",
    benchmark:
        "output/benchmarks/qc/fastp/{sample}_{unit}_{fq}_pipe.{ext}.txt"
    wildcard_constraints:
        sample="|".join(sample_IDs),
        unit="|".join(unit_names),
        ext=r"fq|fq\.gz|fastq|fastq\.gz",
    threads: 1
    conda:
        "../envs/utils.yaml"
    shell:
        """cat {input} > {output} 2> {log}"""


rule fastp_pe:
    """
    Trim paired end reads with fastp.
    """
    input:
        sample=get_fastp_pe_input,
    output:
        trimmed=[
            temp("output/qc/fastp/{sample}_{unit}_R1.fq.gz"),
            temp("output/qc/fastp/{sample}_{unit}_R2.fq.gz"),
        ],
        html="output/qc/fastp/{sample}_{unit}.html",
        json="output/qc/fastp/{sample}_{unit}.json",
    resources:
        mem_mb=get_max_mb(0.5),
    log:
        "output/logs/qc/fastp/{sample}-{unit}.log",
    benchmark:
        "output/benchmarks/qc/fastp/{sample}-{unit}.txt"
    threads: get_threads_per_task_size("small")
    params:
        adapters="",
        extra=" ".join(
            [
                "-M",
                str(config["fastp"]["cut_mean_quality"]),
                "-l",
                str(config["fastp"]["length_required"]),
                str(config["fastp"]["extra"]),
            ]
        ),
    wildcard_constraints:
        sample="|".join(sample_IDs),
    wrapper:
        get_wrapper("fastp")


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
    resources:
        mem_mb=10000 if get_max_mb(0.5) > 10000 else get_max_mb(0.5),
    wildcard_constraints:
        sample="|".join(sample_IDs),
        unit="|".join(unit_names),
        read="single|R1|R2",
    log:
        "output/logs/qc/fastqc_raw/{sample}-{unit}-{read}.log",
    benchmark:
        "output/benchmarks/qc/fastqc_raw/{sample}-{unit}-{read}.txt"
    threads: get_threads_per_task_size("medium")
    wrapper:
        get_wrapper("fastqc")


rule fastqc_trimmed:  # qc on trimmed reads
    input:
        get_fastqc_input_trimmed,
    output:
        zip="output/qc/fastqc/{sample}-{unit}-{read}-trimmed_fastqc.zip",
        html="output/qc/fastqc/{sample}-{unit}-{read}-trimmed.html",
    resources:
        mem_mb=10000 if get_max_mb(0.5) > 10000 else get_max_mb(0.5),
    params:
        "--quiet",
    wildcard_constraints:
        sample="|".join(sample_IDs),
        unit="|".join(unit_names),
        read="single|R1|R2",
    log:
        "output/logs/qc/fastqc_trimmed/{sample}-{unit}-{read}.log",
    benchmark:
        "output/benchmarks/qc/fastqc_trimmed/{sample}-{unit}-{read}.txt"
    threads: get_threads_per_task_size("medium")
    wrapper:
        get_wrapper("fastqc")


# This rule is deprecated. It will not run by default.
# There is no point in running fastqc on the merged reads,
# as they are merely the concatenation of multiple units.
rule fastqc_merged:  # qc on trimmed, merged reads
    input:
        get_fastqc_input_merged,
    output:
        zip="output/qc/fastqc/{sample}-{read}-merged_fastqc.zip",
        html="output/qc/fastqc/{sample}-{read}-merged.html",
    resources:
        mem_mb=10000 if get_max_mb(0.5) > 10000 else get_max_mb(0.5),
    params:
        "--quiet",
    wildcard_constraints:
        sample="|".join(sample_IDs),
        read="single|R1|R2",
    log:
        "output/logs/qc/fastqc_merged/{sample}-{read}.log",
    benchmark:
        "output/benchmarks/qc/fastqc_merged/{sample}-{read}.txt"
    threads: get_threads_per_task_size("medium")
    wrapper:
        get_wrapper("fastqc")


rule multiqc:
    input:
        get_multiqc_input(),
    output:
        report=report("output/qc/multiqc.html", category="QC"),
    log:
        "output/logs/qc/multiqc.log",
    benchmark:
        "output/benchmarks/qc/multiqc.txt"
    wrapper:
        get_wrapper("multiqc")


rule host_removal_create_index:
    input:
        host_removal_reference=config["host_removal"]["reference"],
    output:
        host_removal_index="output/qc/host_removal_reference_db.mmi",
    resources:
        mem_mb=get_max_mb(),
    log:
        "output/logs/qc/host_removal_create_index.log",
    benchmark:
        "output/benchmarks/qc/host_removal_create_index.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        minimap2 -d {output} {input} &> {log}
        """


# Would be nice to have some sort of report after this one.
rule host_removal:
    input:
        fastqs=get_fastqc_input_merged
        if is_activated("merge_reads")
        else get_fastqc_input_trimmed,
        reference="output/qc/host_removal_reference_db.mmi",
    output:
        unpaired=temp("output/qc/filtered/{sample}_unpaired_{read}.fq"),
    params:
        preset="sr",
    wildcard_constraints:
        sample="|".join(sample_IDs),
        read="single|R1|R2",
    threads: get_threads_per_task_size("big")
    resources:
        mem_mb=get_mb_per_cores,
    log:
        "output/logs/qc/host_removal/{sample}-{read}.log",
    benchmark:
        "output/benchmarks/qc/host_removal/{sample}-{read}.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        {{ minimap2 -t {threads}           \
                 -ax {params.preset}       \
                 {input.reference}         \
                 {input.fastqs} ; }}                                2>> {log}   | 
        {{ samtools view -buSh -f 4 ; }}                            2>> {log}   |
        {{ samtools fastq - > {output.unpaired} ; }}                2>> {log}
        """


rule fastq_pair:
    input:
        unpaired=expand(
            "output/qc/filtered/{{sample}}_unpaired_{read}.fq", read=["R1", "R2"]
        ),
    output:
        paired=temp(
            expand(
                ("output/qc/filtered/{{sample}}_unpaired_{read}.fq.paired.fq"),
                read=["R1", "R2"],
            )
        ),
    resources:
        mem_mb=get_max_mb(),
    log:
        "output/logs/qc/host_removal/{sample}-pairing.log",
    benchmark:
        "output/benchmarks/qc/host_removal/{sample}-pairing.txt"
    conda:
        "../envs/fastq-pair.yaml"
    shell:
        """
        fastq_pair {input.unpaired} 2>> {log}      
        """


rule compress_paired:
    input:
        "output/qc/filtered/{sample}_unpaired_{read}.fq.paired.fq",
    output:
        "output/qc/filtered/{sample}_filtered_{read}.fq.gz",
    resources:
        mem_mb=get_mb_per_cores,
    threads: get_threads_per_task_size("medium")
    log:
        "output/logs/qc/host_removal/{sample}-compress-{read}.log",
    benchmark:
        "output/benchmarks/qc/host_removal/{sample}-compress-{read}.txt"
    conda:
        "../envs/utils.yaml"
    shell:
        """
        {{ pigz -p {threads} -f -c {input} > {output} ; }} &> {log}
        """
