"""
Preprocess rules:
    - cutadapt_pipe: get pipe input for cutadapt
    - cutadapt_pe: trim paired end reads with cutadapt
    - merge_fastqs: merge files from different lanes in the same sample with cat
    - fastqc_raw: check quality of raw reads with FastQC
    - fastqc_merged: check quality of trimmed and merged reads with FastQC
    - multiqc: combine reports of rules 'fastqc_raw' and 'fastqc_clean' with MultiQC
"""


rule cutadapt_pipe:
    input:
        get_cutadapt_pipe_input,
    output:
        pipe("pipe/qc/cutadapt/{sample}_{unit}_{fq}.{ext}"),
    log:
        "output/logs/qc/cutadapt/{sample}_{unit}_{fq}_pipe.{ext}.log",
    benchmark:
        "output/benchmarks/qc/cutadapt/{sample}_{unit}_{fq}_pipe.{ext}.txt",
    wildcard_constraints:
        ext=r"fq|fq\.gz|fastq|fastq\.gz",
    threads: 1
    conda:
        "../envs/bash.yaml"
    shell:
        "cat {input} > {output} 2> {log}"


rule cutadapt_pe:
    input:
        get_cutadapt_input,
    output:
        fastq1="output/qc/cutadapt/{sample}_{unit}_R1.fq.gz",
        fastq2="output/qc/cutadapt/{sample}_{unit}_R2.fq.gz",
        qc="output/qc/cutadapt/{sample}_{unit}.paired.qc.txt",
    log:
        "output/logs/qc/cutadapt/{sample}-{unit}.log",
    benchmark:
        "output/benchmarks/qc/cutadapt/{sample}-{unit}.txt",
    params:
        others="",  # config["params"]["cutadapt-pe"],
        adapters="-a AGAGCACACGTCTGAACTCCAGTCAC -g AGATCGGAAGAGCACACGT -A AGAGCACACGTCTGAACTCCAGTCAC -G AGATCGGAAGAGCACACGT",
        extra="--minimum-length 1 -q 20",
        # adapters=lambda w: str(units.loc[w.sample].loc[w.unit, "adapters"]),
    threads: 1
    wrapper:
        "0.59.2/bio/cutadapt/pe"


rule merge_fastqs:
    input:
        get_fastqs,
    output:
        "output/qc/merged/{sample}_{read}.fq.gz",
    log:
        "output/logs/qc/merge_fastqs/{sample}.{read}.log",
    benchmark:
        "output/benchmarks/qc/merge_fastqs/{sample}.{read}.txt",
    wildcard_constraints:
        read="single|R1|R2",
    conda:
        "../envs/bash.yaml"
    shell:
        "cat {input} > {output} 2> {log}"


rule fastqc_raw:  # qc on raw, unmerged reads
    input:
        get_fastqc_input_raw,
    output:
        zip="{output}/qc/fastqc/{sample}-{unit}-{read}_fastqc.zip",
        html="{output}/qc/fastqc/{sample}-{unit}-{read}.html",
    params:
        "--quiet",
    log:
        "{output}/logs/qc/fastqc/{sample}-{unit}-{read}-fastqc.log",
    benchmark:
        "{output}/benchmarks/qc/fastqc/{sample}-{unit}-{read}-fastqc.txt"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"


rule fastqc_merged:  # qc on trimmed, merged reads
    input:
        get_fastqc_input_merged,
    output:
        zip="{output}/qc/fastqc/{sample}-{read}_fastqc.zip",
        html="{output}/qc/fastqc/{sample}-{read}.html",
    params:
        "--quiet",
    log:
        "{output}/logs/qc/fastqc/{sample}-{read}-fastqc.log",
    benchmark:
        "{output}/benchmarks/qc/fastqc/{sample}-{read}-fastqc.txt"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"


rule multiqc:
    input:
        get_multiqc_input,
    output:
        report="output/qc/multiqc.html",
    log:
        "output/logs/qc/multiqc.log",
    benchmark:
        "output/benchmarks/qc/multiqc.txt"
    wrapper:
        "0.77.0/bio/multiqc"
