"""
Assembly rules:

    - megahit: assemble preprocessed reads with Megahit
"""

from pathlib import Path


rule megahit:
    input:
        host_removal_output="{output}/preprocess/interleave/{sample}-clean.fq",
    output:
        contigs="{output}/assembly/megahit/{sample}/{sample}.contigs.fa",
    params:
        out_dir=lambda w, output: str(Path(output.contigs).parent.parent),  # this is equivalent to "{output}/megahit"
        min_contig_len=200,
        k_list="21,29",
        memory=0.5,
    threads: int(workflow.cores * 0.75)
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

        megahit -r {input} -o {params.out_dir}/{wildcards.sample} \
                --out-prefix {wildcards.sample} \
                --min-contig-len {params.min_contig_len}  \
                -t {threads}  \
                -m {params.memory} \
                --k-list {params.k_list} &> {log}
        """
