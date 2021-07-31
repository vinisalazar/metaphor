"""
Preprocess rules:

    - fastqc_raw: check quality of raw reads with FastQC
    - flash: extend and combine reads with Flash
    - interleave: interleave merged reads and extended fragments with custom script
    - hostremoval: map reads against reference database with ???
    - fastqc_clean: check quality after previous steps with FastQC
    - multiqc: combine reports of rules 'fastqc_raw' and 'fastqc_clean' with MultiQC
"""


rule fastqc_raw:  # qc on raw reads
    input:
        "data/{sample}.fq"
    output:
        zip="{output}/preprocess/fastqc/{sample}_fastqc.zip",
        html="{output}/preprocess/fastqc/{sample}.html",
    params:
        "--quiet",
    log:
        "{output}/logs/preprocess/fastqc/{sample}-fastqc.log",
    benchmark:
        "{output}/benchmarks/preprocess/fastqc/{sample}_fastqc.txt"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"


rule fastqc_clean:  # qc on clean reads
    input:
        "{output}/preprocess/interleave/{sample}-clean.fq"
    output:
        zip="{output}/preprocess/fastqc/{sample}-clean_fastqc.zip",
        html="{output}/preprocess/fastqc/{sample}-clean.html",
    params:
        "--quiet",
    log:
        "{output}/logs/preprocess/fastqc/{sample}-fastqc.log",
    benchmark:
        "{output}/benchmarks/preprocess/fastqc/{sample}_fastqc.txt"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"


rule flash:
    input:
        fqforward=fq1,
        fqreverse=fq2,
    output:
        flash_notcombined1="{output}/preprocess/flash/{sample}.notCombined_1.fastq",
        flash_notcombined2="{output}/preprocess/flash/{sample}.notCombined_2.fastq",
        flash_extended="{output}/preprocess/flash/{sample}.extendedFrags.fastq",
    params:
        max_overlap=120,
    log:
        "{output}/logs/preprocess/flash/{sample}-flash.log",
    benchmark:
        "{output}/benchmarks/preprocess/flash/{sample}.txt"
    conda:
        "../envs/flash.yaml"
    shell:
        """
        flash -d {wildcards.output}/preprocess/flash -o {wildcards.sample} \
              -M {params.max_overlap} {input} &> {log}
        """


rule interleave:
    input:
        flash_notcombined1="{output}/preprocess/flash/{sample}.notCombined_1.fastq",
        flash_notcombined2="{output}/preprocess/flash/{sample}.notCombined_2.fastq",
        flash_extended="{output}/preprocess/flash/{sample}.extendedFrags.fastq",
    output:
        interleaved="{output}/preprocess/interleave/{sample}-interleaved.fq",
        merged="{output}/preprocess/interleave/{sample}-merged.fq",
    log:
        "{output}/logs/preprocess/interleave/{sample}-interleave.log",
    benchmark:
        "{output}/benchmarks/preprocess/interleave/{sample}-interleave.txt"
    conda:
        "../envs/bash.yaml"
    shell:
        """
        {{ bash workflow/scripts/interleave_fastq.sh {input.flash_notcombined1} {input.flash_notcombined2} > {output.interleaved} ; }} &> {log}

        cat {input.flash_extended} {output.interleaved} > {output.merged}
        """


rule hostremoval:
    input:
        flash_merged="{output}/preprocess/interleave/{sample}-merged.fq",
    params:
        alignment_id_threshold=70,
        alignment_coverage_threshold=70,
        dbs="mm1,mm2,mm3,mm4,mm5,mm6",
    output:
        host_removal_output="{output}/preprocess/interleave/{sample}-clean.fq",
    log:
        "{output}/logs/preprocess/interleave/{sample}-hostremoval.log",
    benchmark:
        "{output}/benchmarks/preprocess/interleave/{sample}-hostremoval.txt"
    conda:
        "../envs/bash.yaml"
    shell:
        # This rule is off for now due to the lack of databases
        # "perl bin/dqc/deconseq.pl -dbs {params.dbs}" \
        # " -f {input}" \
        # " -i {params.alignment_id_threshold}" \
        # " -c {params.alignment_coverage_threshold} " \
        # " -id {wildcards.sample}"
        "cat {input} > {output}"


rule multiqc:
    input:
        get_multiqc_input(),
    output:
        report="{output}/preprocess/multiqc.html",
    log:
        "{output}/logs/preprocess/multiqc.log",
    benchmark:
        "{output}/benchmarks/preprocess/multiqc.txt"
    wrapper:
        "0.77.0/bio/multiqc"
