# MetaSnakePipe
# Original MetaGenePipe workflow by Bobbie Shaban
# Snakemake port by Vini Salazar

rule fastqc:
    input:
        "data/{sample}.fq"
    output:
        html="{output}/fastqc/{sample}.html",
        zip="{output}/fastqc/{sample}_fastqc.zip"
    params: "--quiet"
    log:
        "{output}/logs/fastqc/{sample}.log"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"


rule flash:
    input:
        fqforward="data/{sample}-1.fq",
        fqreverse="data/{sample}-2.fq" 
    output: 
        flashNotCombined1="{output}/flash/{sample}.notCombined_1.fastq",
        flashNotCombined2="{output}/flash/{sample}.notCombined_2.fastq",
        flashExtended="{output}/flash/{sample}.extendedFrags.fastq"
    log:

    shell: 
        "mkdir -p {output.output_dir} ; flash -o {output.output_dir}/readsa {input}"


rule interleave:
    input: 
        flashNotCombined1="output/flash/readsa.notCombined_1.fastq",
        flashNotCombined2="output/flash/readsa.notCombined_2.fastq",
        flashExtended="output/flash/readsa.extendedFrags.fastq"
    output: 
        output_dir=directory("output/interleave"),
        interleaved="output/interleave/readsa-interleaved.fq",
        merged="output/interleave/readsa-merged.fq"
    shell: 
        "mkdir {output.output_dir} ; bash scripts/interleave_fastqc {input.flashNotCombined1} {input.flashNotCombined2} > {output.interleaved} ; cat {output.interleaved} {input.flashExtended} > {output.merged}"