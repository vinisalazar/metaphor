"""
Assembly rules:

    - megahit: assemble preprocessed reads with Megahit
"""

from pathlib import Path


rule megahit:
    input:
        fastq1="{output}/qc/merged/{sample}_R1.fq.gz",
        fastq2="{output}/qc/merged/{sample}_R2.fq.gz",
    output:
        contigs="{output}/assembly/megahit/{sample}/{sample}.contigs.fa",
    params:
        out_dir=lambda w, output: str(Path(output.contigs).parent.parent),  # this is equivalent to "{output}/megahit"
        min_contig_len=200,
        k_list="21,29,39,59,79,99,119,141",
    threads: workflow.cores
    log:
        "{output}/logs/assembly/megahit/{sample}-megahit.log",
    benchmark:
        "{output}/benchmarks/assembly/megahit/{sample}.txt"
    conda:
        "../envs/megahit.yaml"
    # Using the '--12' flag yielded slightly better results than the '-r' flag, but also eventually presented errors
    shell:
        """
        # MegaHit has no --force flag, so we must remove the created directory prior to running
        rm -rf {params.out_dir}/{wildcards.sample}

        megahit -1 {input.fastq1} -2 {input.fastq2} \
                -o {params.out_dir}/{wildcards.sample} \
                --out-prefix {wildcards.sample} \
                --min-contig-len {params.min_contig_len}  \
                -t {threads}  \
                --k-list {params.k_list} &> {log}
        """
