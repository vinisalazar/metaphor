# MetaSnakePipe
# Original MetaGenePipe workflow by Bobbie Shaban
# Snakemake port by Vini Salazar

rule fastqc:
    input: 
        forwardfq="data/readsa-1.fq",
        reversefq="data/readsa-2.fq",
    output: 
        output_dir="reads_fastqc"
    shell: 
        "fastqc {input} -o {output.output_dir}"
