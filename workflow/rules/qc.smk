"""
Preprocess rules:

    - fastqc_raw: check quality of raw reads with FastQC
    - flash: extend and combine reads with Flash
    - interleave: interleave clean reads and extended fragments with custom script
    - hostremoval: map reads against reference database with ???
    - fastqc_clean: check quality after previous steps with FastQC
    - multiqc: combine reports of rules 'fastqc_raw' and 'fastqc_clean' with MultiQC
"""


rule fastqc_raw:  # qc on raw reads
    input:
        "data/{sample}.fq"
    output:
        zip="{output}/qc/fastqc/{sample}_fastqc.zip",
        html="{output}/qc/fastqc/{sample}.html",
    params:
        "--quiet",
    log:
        "{output}/logs/qc/fastqc/{sample}-fastqc.log",
    benchmark:
        "{output}/benchmarks/qc/fastqc/{sample}_fastqc.txt"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"


rule fastqc_clean:  # qc on clean reads
    input:
        "{output}/qc/interleave/{sample}-clean.fq"
    output:
        zip="{output}/qc/fastqc/{sample}-clean_fastqc.zip",
        html="{output}/qc/fastqc/{sample}-clean.html",
    params:
        "--quiet",
    log:
        "{output}/logs/qc/fastqc/{sample}-fastqc.log",
    benchmark:
        "{output}/benchmarks/qc/fastqc/{sample}_fastqc.txt"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"


rule flash:
    input:
        fqforward=fq1,
        fqreverse=fq2,
    output:
        flash_notcombined1="{output}/qc/flash/{sample}.notCombined_1.fastq",
        flash_notcombined2="{output}/qc/flash/{sample}.notCombined_2.fastq",
        flash_extended="{output}/qc/flash/{sample}.extendedFrags.fastq",
    params:
        max_overlap=120,
    log:
        "{output}/logs/qc/flash/{sample}-flash.log",
    benchmark:
        "{output}/benchmarks/qc/flash/{sample}.txt"
    conda:
        "../envs/flash.yaml"
    shell:
        """
        flash -d {wildcards.output}/qc/flash -o {wildcards.sample} \
              -M {params.max_overlap} {input} &> {log}
        """


rule interleave:
    input:
        flash_notcombined1="{output}/qc/flash/{sample}.notCombined_1.fastq",
        flash_notcombined2="{output}/qc/flash/{sample}.notCombined_2.fastq",
        flash_extended="{output}/qc/flash/{sample}.extendedFrags.fastq",
    output:
        interleaved="{output}/qc/interleave/{sample}-interleaved.fq",
        clean="{output}/qc/interleave/{sample}-clean.fq",
    log:
        "{output}/logs/qc/interleave/{sample}-interleave.log",
    benchmark:
        "{output}/benchmarks/qc/interleave/{sample}-interleave.txt"
    conda:
        "../envs/bash.yaml"
    shell:
        """
        {{ bash workflow/scripts/interleave_fastq.sh {input.flash_notcombined1} {input.flash_notcombined2} > {output.interleaved} ; }} &> {log}

        cat {input.flash_extended} {output.interleaved} > {output.clean}
        """


rule multiqc:
    input:
        get_multiqc_input(),
    output:
        report="{output}/qc/multiqc.html",
    log:
        "{output}/logs/qc/multiqc.log",
    benchmark:
        "{output}/benchmarks/qc/multiqc.txt"
    wrapper:
        "0.77.0/bio/multiqc"
