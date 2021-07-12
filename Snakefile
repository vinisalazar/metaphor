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
        output_dir=directory("output/flash")
    shell: 
        "mkdir -p {output.output_dir} ; flash -o {output.output_dir}/readsa {input}"
