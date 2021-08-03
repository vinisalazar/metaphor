"""
Preprocess rules:

    - fastqc_raw: check quality of raw reads with FastQC
    - flash: extend and combine reads with Flash
    - interleave: interleave clean reads and extended fragments with custom script
    - hostremoval: map reads against reference database with ???
    - fastqc_clean: check quality after previous steps with FastQC
    - multiqc: combine reports of rules 'fastqc_raw' and 'fastqc_clean' with MultiQC
"""

rule trim_pipe:
    input:
        get_trim_pipe_input,
    output:
        pipe("pipe/qc/trim/{sample}_{unit}_{fq}.{ext}"),
    log:
        "output/logs/qc/trim_pipe/{sample}_{unit}_{fq}.{ext}.log",
    wildcard_constraints:
        ext=r"fq|fq\.gz|fastq|fastq\.gz",
    threads: 1
    shell:
        "cat {input} > {output} 2> {log}"


rule trim_galore_pe:
    input:
        get_trim_input,
    output:
        fastq1="output/qc/trim/{sample}_{unit}_R1_trim.fq.gz",
        qc_fq1="output/qc/trim/{sample}_{unit}_R1.gz_trimming_report.txt",
        fastq2="output/qc/trim/{sample}_{unit}_R2_trim.fq.gz",
        qc_fq2="output/qc/trim/{sample}_{unit}_R2.gz_trimming_report.txt"
    params:
        extra="--illumina -q 20"
    log:
        "output/logs/qc/trim_galore/{sample}_{unit}.log"
    benchmark:
        "output/benchmarks/qc/trim_galore/{sample}_{unit}.log"
    wrapper:
        "0.77.0/bio/trim_galore/pe"


# rule trim_galore_se:
#     input:
#         get_trim_input,
#     output:
#         fastq="output/qc/trim/{sample}_{unit}.trim.fq.gz",
#         qc="output/qc/trim/{sample}_{unit}.fq.gz_trimming_report.txt",
#     params:
#         extra="--illumina -q 20"
#     log:
#         "output/logs/qc/trim_galore/{sample}_{unit}.log"
#     benchmark:
#         "output/benchmarks//qc/trim_galore/{sample}_{unit}.log"
#     wrapper:
#         "0.77.0/bio/trim_galore/se"
    

rule merge_fastqs:
    input:
        get_fastqs,
    output:
        "output/qc/merged/{sample}_{read}.fq.gz",
    log:
        "output/logs/qc/merge_fastqs/{sample}.{read}.log",
    wildcard_constraints:
        read="single|R1|R2",
    shell:
        "cat {input} > {output} 2> {log}"


rule fastqc:  # qc on raw reads
    input:
        "output/qc/{kind}/{basename}.fq.gz"
    output:
        zip="{output}/qc/fastqc/{basename}-{kind}_fastqc.zip",
        html="{output}/qc/fastqc/{basename}-{kind}.html",
    params:
        "--quiet",
    log:
        "{output}/logs/qc/fastqc/{basename}-{kind}-fastqc.log",
    benchmark:
        "{output}/benchmarks/qc/fastqc/{basename}-{kind}_fastqc.txt"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"


rule multiqc:
    input:
        get_multiqc_input,
    output:
        report="{output}/qc/multiqc.html",
    log:
        "{output}/logs/qc/multiqc.log",
    benchmark:
        "{output}/benchmarks/qc/multiqc.txt"
    wrapper:
        "0.77.0/bio/multiqc"
