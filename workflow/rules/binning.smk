"""
Binning rules:
    - vamb: Bin contigs with vamb
    - MetaBAT2: Bin contigs with MetaBat2
"""
from pathlib import Path


rule vamb:
    input:
        bam_contig_depths="output/mapping/bam_contig_depths.txt",
        catalogue="output/mapping/catalogue.fna.gz",
    output:
        get_vamb_output(),
    params:  # defaults in vamb's README
        outdir=lambda w, output: get_parent(output[0]),
        binsplit_sep="C",
        minfasta=200000,
        batchsize=256,
    log:
        "output/logs/binning/vamb.log",
    benchmark:
        "output/benchmarks/binning/vamb.txt"
    conda:
        "../envs/vamb.yaml"
    shell:
        """
        rm -rf {params.outdir}

        vamb --outdir {params.outdir}           \
             --fasta {input.catalogue}          \
             --jgi {input.bam_contig_depths}    \
             -o {params.binsplit_sep}           \
             -t {params.batchsize}              \
             --minfasta {params.minfasta} &> {log}
        """


rule metabat2:
    input:
        contigs="output/mapping/catalogue.fna.gz",
        depths="output/mapping/bam_contig_depths.txt",
    output:
        outdir=directory("output/binning/metabat2")
    params:
        minContig=2500,
        seed=config["metabat2"]["seed"],
        outfile=lambda w, output: output.outdir + "/bin"
    threads: round(workflow.cores * 0.75)
    log:
        "output/logs/binning/metabat2.log",
    benchmark:
        "output/benchmarks/binning/metabat2.txt"
    conda:
        "../envs/metabat2.yaml"
    shell:
        """
        rm -rf {output} && mkdir {output}

        metabat2 -i {input.contigs}             \
                 -a {input.depths}              \
                 -m {params.minContig}          \
                 -t {threads}                   \
                 --seed {params.seed}           \
                 -o {params.outfile} &> {log}
        """
