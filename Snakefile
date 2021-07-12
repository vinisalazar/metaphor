# MetaSnakePipe
# Original MetaGenePipe workflow by Bobbie Shaban
# Snakemake port by Vini Salazar

rule fastqc:
    input: 
        fqforward="data/readsa-1.fq",
        fqreverse="data/readsa-2.fq"
    output: 
        output_dir=directory("output/fastqc"),
    shell: 
        "mkdir -p {output.output_dir} ; fastqc {input} -o {output.output_dir}"


rule flash:
    input:
        fqforward="data/readsa-1.fq",
        fqreverse="data/readsa-2.fq" 
    output: 
        output_dir=directory("output/flash"),
        flashNotCombined1="output/flash/readsa.notCombined_1.fastq",
        flashNotCombined2="output/flash/readsa.notCombined_2.fastq",
        flashExtended="output/flash/readsa.extendedFrags.fastq"
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