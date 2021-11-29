"""
Assembly rules:
    - concatenate_merged_reads: prepare MegaHIT input if there coassembly is True in config file
    - megahit: assemble preprocessed reads with Megahit
    - metaquast: evaluate assembly results with MetaQuast
"""


rule concatenate_merged_reads:
    input:
        R1=expand("output/qc/merged/{sample}_R1.fq.gz", sample=sample_IDs),
        R2=expand("output/qc/merged/{sample}_R2.fq.gz", sample=sample_IDs),
    output:
        R1_concat="output/qc/merged/all_samples_R1.fq.gz",
        R2_concat="output/qc/merged/all_samples_R2.fq.gz",
    log:
        "output/logs/assembly/concatenate_merged_reads.log",
    benchmark:
        "output/benchmarks/assembly/concatenate_merged_reads.txt"
    conda:
        "../envs/bash.yaml"
    shell:
        """
        {{ cat {input.R1} > {output.R1_concat} ; }} > {log}
        {{ cat {input.R2} > {output.R2_concat} ; }} >> {log}
        """


rule megahit:
    input:
        fastq1="output/qc/merged/{sample}_R1.fq.gz",
        fastq2="output/qc/merged/{sample}_R2.fq.gz",
    output:
        contigs="output/assembly/megahit/{sample}/{sample}.contigs.fa",
    params:
        out_dir=lambda w, output: get_parent(get_parent(output.contigs)),  # this is equivalent to "{output}/megahit"
        min_contig_len=200,
        k_list="21,29,39,59,79,99,119,141",
        preset=config["megahit"]["preset"],
    threads: round(workflow.cores * 0.75)
    log:
        "output/logs/assembly/megahit/{sample}.log",
    benchmark:
        "output/benchmarks/assembly/megahit/{sample}.txt"
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


rule megahit_coassembly:
    input:
        fastq1="output/qc/merged/all_samples_R1.fq.gz",
        fastq2="output/qc/merged/all_samples_R2.fq.gz",
    output:
        contigs="output/assembly/megahit/coassembly.contigs.fa",
    params:
        out_dir=lambda w, output: get_parent(output.contigs),
        min_contig_len=200,
        k_list="21,29,39,59,79,99,119,141",
        preset=config["megahit"]["preset"],
    threads: round(workflow.cores * 0.75)
    log:
        "output/logs/assembly/megahit/coassembly.log",
    benchmark:
        "output/benchmarks/assembly/megahit/coassembly.txt"
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        # MegaHit has no --force flag, so we must remove the created directory prior to running
        rm -rf {params.out_dir}

        megahit -1 {input.fastq1} -2 {input.fastq2}         \
                -o {params.out_dir}                         \
                --presets {params.preset}                   \
                --out-prefix coassembly                     \
                --min-contig-len {params.min_contig_len}    \
                -t {threads}                                \
                --k-list {params.k_list} &> {log}
        """


rule metaquast:
    input:
        assemblies=get_contigs_input()
        if config["coassembly"]
        else expand(
            "output/assembly/megahit/{sample}/{sample}.contigs.fa", sample=sample_IDs
        ),
    output:
        outfile=get_metaquast_output(),
    params:
        mincontig=500,
        reference=get_metaquast_reference(),
        labels=",".join(sample_IDs),
        outdir=lambda w, output: str(Path(output.outfile).parent.parent),
        extra_params="--no-icarus",
    threads: round(workflow.cores * 0.75)
    log:
        "output/logs/assembly/metaquast.log",
    benchmark:
        "output/benchmarks/assembly/metaquast.txt"
    conda:
        "../envs/quast.yaml"
    shell:
        """
        metaquast.py -t {threads}               \
                     -o {params.outdir}         \
                     -l {params.labels}         \
                     -m {params.mincontig}      \
                     --no-icarus                \
                     {params.reference}         \
                     {input.assemblies}
        """
