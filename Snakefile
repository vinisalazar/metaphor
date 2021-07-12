# MetaSnakePipe
# Original MetaGenePipe workflow by Bobbie Shaban
# Snakemake port by Vini Salazar

rule fastqc:
    input: 
        "data/readsa-1.fq",
        "data/readsa-2.fq"
    output: 
        "readsa_fastqc"
    shell: 
        "fastqc {input} -o {output}"
