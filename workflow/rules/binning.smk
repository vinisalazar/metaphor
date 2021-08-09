"""
Binning rules:
    - vamb: Bin contigs with vamb
    - maxbin2: Bin contigs with MaxBin2
"""


rule vamb:
    input:
        bam_contig_depths="output/mapping/bam_contig_depths.txt",
        catalogue="output/mapping/catalogue.fna.gz",
    output:
        clusters="output/binning/vamb/clusters.tsv",
    params:  # defaults in vamb's README
        outdir=lambda w, output: get_parent(output.clusters),
        binsplit_sep="C",
        minfasta=200000,
    log:
        "output/log/binning/vamb.log",
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
             --minfasta {params.minfasta} &> {log}
        """
