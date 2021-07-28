"""
Assembly rules:

    - megahit: assemble preprocessed reads with Megahit
    - vamb: refine contigs with vamb
"""

from pathlib import Path


rule megahit:
    input:
        host_removal_output="{output}/preprocessing/interleave/{sample}-clean.fq",
    output:
        contigs="{output}/assembly/megahit/{sample}/{sample}.contigs.fa",
    params:
        out_dir=lambda w, output: str(Path(output.contigs).parent.parent),  # this is equivalent to "{output}/megahit"
        min_contig_len=200,
        k_list="21,29",
        memory=0.5,
        cpus=workflow.cores,
    log:
        "{output}/logs/assembly/megahit/{sample}-megahit.log",
    benchmark:
        "{output}/benchmarks/assembly/megahit/{sample}.txt"
    conda:
        "../envs/megahit.yaml"
    # Using the '--12' flag yielded slightly better results than the '-r' flag
    shell:
        """
        # MegaHit has no --force flag, so we must remove the created directory prior to running
        rm -rf {params.out_dir}/{wildcards.sample}

        megahit --12 {input} -o {params.out_dir}/{wildcards.sample} \
                --out-prefix {wildcards.sample} \
                --min-contig-len {params.min_contig_len}  \
                -t {params.cpus}  \
                -m {params.memory} \
                --k-list {params.k_list} &> {log}
        """


rule vamb:
    input:
        bamfiles=expand(
            "output/mapping/bam/{sample}",
            sample=[
                "readsa",
            ],
        ),
        catalogue="{output}/mapping/catalogue.fna.gz",
    output:
        outdir=directory("output/mapping/vamb"),
    params:  # defaults in vamb's README
        binsplit_sep="C",
        minfasta=200000,
    log:
        "output/log/mapping/vamb.log",
    benchmark:
        "output/benchmarks/mapping/vamb.txt"
    conda:
        "../envs/vamb.yaml"
    shell:
        """
        rm -rf {output}

        vamb --outdir {output} \
             --fasta {input.catalogue} \
             --bamfiles {input.bamfiles} \
             -o {params.binsplit_sep} \
             --minfasta {params.minfasta} &> {log}
        """
