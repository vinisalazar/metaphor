"""
Preprocessing rules:

    - fastqc_raw: check quality of raw reads with FastQC
    - flash: extend and combine reads with Flash
    - interleave: interleave merged reads and extended fragments with custom script
    - hostremoval: map reads against reference database with ???
    - fastqc_clean: check quality after previous steps with FastQC
    - multiqc: combine reports of rules 'fastqc_raw' and 'fastqc_clean' with MultiQC
"""


rule fastqc_raw:  # qc on raw reads
    input:
        expand("{sample}", sample=raw_fqs),
    output:
        html="{output}/qc/{sample}.html",
        zip="{output}/qc/{sample}_fastqc.zip",
    params:
        "--quiet",
    log:
        "{output}/logs/qc/{sample}-fastqc.log",
    benchmark:
        "{output}/benchmarks/qc/{sample}_fastqc.txt"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"


rule flash:
    input:
        fqforward="data/{sample}-1.fq",
        fqreverse="data/{sample}-2.fq",
    output:
        flash_notcombined1="{output}/flash/{sample}.notCombined_1.fastq",
        flash_notcombined2="{output}/flash/{sample}.notCombined_2.fastq",
        flash_extended="{output}/flash/{sample}.extendedFrags.fastq",
    params:
        max_overlap=120,
    log:
        "{output}/logs/qc/{sample}-flash.log",
    benchmark:
        "{output}/benchmarks/flash/{sample}.txt"
    conda:
        "../envs/flash.yaml"
    shell:
        """
        flash -d {wildcards.output}/flash -o {wildcards.sample} \
              -M {params.max_overlap} {input} &> {log}
        """


rule interleave:
    input:
        flash_notcombined1="{output}/flash/{sample}.notCombined_1.fastq",
        flash_notcombined2="{output}/flash/{sample}.notCombined_2.fastq",
        flash_extended="{output}/flash/{sample}.extendedFrags.fastq",
    output:
        interleaved="{output}/interleave/{sample}-interleaved.fq",
        merged="{output}/interleave/{sample}-merged.fq",
    log:
        "{output}/logs/interleave/{sample}-interleave.log",
    benchmark:
        "{output}/benchmarks/interleave/{sample}-interleave.txt"
    conda:
        "../envs/bash.yaml"
    shell:
        """
        {{ bash ../scripts/interleave_fastq.sh {input.flash_notcombined1} {input.flash_notcombined2} > {output.interleaved} ; }} &> {log}

        cat {input.flash_extended} {output.interleaved} > {output.merged}
        """


rule hostremoval:
    input:
        flash_merged="{output}/interleave/{sample}-merged.fq",
    params:
        alignment_id_threshold=70,
        alignment_coverage_threshold=70,
        dbs="mm1,mm2,mm3,mm4,mm5,mm6",
    output:
        host_removal_output="{output}/interleave/{sample}-clean.fq",
    log:
        "{output}/logs/interleave/{sample}-hostremoval.log",
    benchmark:
        "{output}/benchmarks/interleave/{sample}-hostremoval.txt"
    conda:
        "../envs/bash.yaml"
    shell:
        "cat {input} > {output}"
        # This rule is off for now due to the lack of databases
        # "perl bin/dqc/deconseq.pl -dbs {params.dbs}" \
        # " -f {input}" \
        # " -i {params.alignment_id_threshold}" \
        # " -c {params.alignment_coverage_threshold} " \
        # " -id {wildcards.sample}"


rule fastqc_clean:  # qc on merged reads, after rules 'flash', 'interleave', and 'hostremoval'
    input:
        "{output}/interleave/{sample}-clean.fq",
    output:
        html="{output}/qc/{sample}-clean.html",
        zip="{output}/qc/{sample}-clean.zip",
    params:
        "--quiet",
    log:
        "{output}/logs/qc/{sample}-fastqc_clean.log",
    benchmark:
        "{output}/benchmarks/qc/{sample}_fastqc_clean.txt"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"


rule multiqc:
    input:
        expand(
            "output/qc/{sample}-{kind}_fastqc.zip",
            sample=[
                "readsa",
            ],
            kind=["1", "2", "merged"],
        ),
    output:
        report="{output}/qc/multiqc.html",
    log:
        "{output}/logs/qc/multiqc.log",
    benchmark:
        "{output}/benchmarks/qc/multiqc.txt"
    wrapper:
        "0.77.0/bio/multiqc"