# MetaSnakePipe
# Original MetaGenePipe workflow by Bobbie Shaban
# Snakemake port by Vini Salazar

rule all:
    input: 
        expand("data/{sample}-{readsno}.fq",
                sample=["readsa",],
                readsno=["1", "2"]) 


rule fastqc:
    input:
        "data/{sample}-{readsno}.fq"
    output:
        html="{output}/fastqc/{sample}-{readsno}.html",
        zip="{output}/fastqc/{sample}-{readsno}_fastqc.zip"
    params: "--quiet"
    log:
        "{output}/logs/fastqc/{sample}-{readsno}.log"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"


rule flash:
    input:
        fqforward="data/{sample}-1.fq",
        fqreverse="data/{sample}-2.fq" 
    output: 
        flash_notcombined1="{output}/flash/{sample}_notcombined_1.fastq",
        flash_notcombined2="{output}/flash/{sample}_notcombined_2.fastq",
        flash_extended="{output}/flash/{sample}_extendedfrags.fastq"
    shell: 
        "flash -d {wildcards.output}/flash -o {wildcards.sample} {input}"


rule interleave:
    input: 
        flash_notcombined1="{output}/flash/{sample}_notcombined_1.fastq",
        flash_notcombined2="{output}/flash/{sample}_notcombined_2.fastq",
        flash_extended="{output}/flash/{sample}_extendedfrags.fastq"
    output: 
        interleaved="{output}/interleave/{sample}-interleaved.fq",
        merged="{output}/interleave/{sample}-merged.fq"
    shell: 
        "bash scripts/interleave_fastqc.sh {input.flash_notcombined1} {input.flash_notcombined2} > {output.interleaved} ; cat {output.interleaved} {input.flash_extended} > {output.merged}"
