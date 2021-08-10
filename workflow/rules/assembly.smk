"""
Assembly rules:

    - megahit: assemble preprocessed reads with Megahit
"""


rule megahit:
    input:
        fastq1="{output}/qc/merged/{sample}_R1.fq.gz",
        fastq2="{output}/qc/merged/{sample}_R2.fq.gz",
    output:
        contigs="{output}/assembly/megahit/{sample}/{sample}.contigs.fa",
    params:
        out_dir=lambda w, output: get_parent(get_parent(output.contigs)),  # this is equivalent to "{output}/megahit"
        min_contig_len=200,
        k_list="21,29,39,59,79,99,119,141",
        preset=config["megahit"]["preset"],
    threads: round(workflow.cores * 0.75)
    log:
        "{output}/logs/assembly/megahit/{sample}-megahit.log",
    benchmark:
        "{output}/benchmarks/assembly/megahit/{sample}.txt"
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        # MegaHit has no --force flag, so we must remove the created directory prior to running
        rm -rf {params.out_dir}/{wildcards.sample}

        megahit -1 {input.fastq1} -2 {input.fastq2}         \
                -o {params.out_dir}/{wildcards.sample}      \
                --presets {params.preset}                   \
                --out-prefix {wildcards.sample}             \
                --min-contig-len {params.min_contig_len}    \
                -t {threads}                                \
                --k-list {params.k_list} &> {log}
        """
