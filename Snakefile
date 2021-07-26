# MetaSnakePipe
# Original MetaGenePipe workflow by Bobbie Shaban
# Snakemake port by Vini Salazar

rule all:
    input: 
        expand("data/{sample}-{readsno}.fq",
                sample=["readsa",],
                readsno=["1", "2"]) 


rule fastqc_raw:  # qc on raw reads
    input:
        "data/{sample}-{readsno}.fq"
    output:
        html="{output}/qc/{sample}-{readsno}.html",
        zip="{output}/qc/{sample}-{readsno}_fastqc.zip"
    params: "--quiet"
    log:
        "{output}/logs/qc/{sample}-{readsno}.log"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"


rule flash:
    input:
        fqforward="data/{sample}-1.fq",
        fqreverse="data/{sample}-2.fq" 
    output: 
        flash_notcombined1="{output}/flash/{sample}.notCombined_1.fastq",
        flash_notcombined2="{output}/flash/{sample}.notCombined_2.fastq",
        flash_extended="{output}/flash/{sample}.extendedFrags.fastq"
    shell: 
        "flash -d {wildcards.output}/flash -o {wildcards.sample} {input}"


rule interleave:
    input: 
        flash_notcombined1="{output}/flash/{sample}.notCombined_1.fastq",
        flash_notcombined2="{output}/flash/{sample}.notCombined_2.fastq",
        flash_extended="{output}/flash/{sample}.extendedFrags.fastq"
    output: 
        interleaved="{output}/interleave/{sample}-interleaved.fq",
        merged="{output}/interleave/{sample}-merged.fq"
    shell: 
        "bash scripts/interleave_fastqc.sh {input.flash_notcombined1} {input.flash_notcombined2} " \
        " > {output.interleaved} ; cat {output.interleaved} {input.flash_extended} > {output.merged}"


rule hostremoval:
    input:
        flash_merged="{output}/interleave/{sample}-merged.fq"
    params:
        alignment_id_threshold=70,
        alignment_coverage_threshold=70,
        dbs="mm1,mm2,mm3,mm4,mm5,mm6"
    output:
        host_removal_output="{output}/interleave/{sample}-clean.fq"
    shell:
        "cat {input} > {output}"
        # This rule is off for now due to the lack of databases
        # "perl bin/dqc/deconseq.pl -dbs {params.dbs}" \
        # " -f {input}" \
        # " -i {params.alignment_id_threshold}" \
        # " -c {params.alignment_coverage_threshold} " \
        # " -id {wildcards.sample}"


rule megahit:
    input:
        host_removal_output="{output}/interleave/{sample}-clean.fq"
    output:
        megahit_out=directory("{output}/megahit/{sample}")
    params:
        out_preffix="{sample}",
        min_contig_len=200,
        memory=0.5
    shell:
        "megahit -r {input} -o {output} --out-prefix {params.out_preffix} " \
        "--min-contig-len {params.min_contig_len} " \
        "-t {workflow.cores} " \
        "-m {params.memory} " \
