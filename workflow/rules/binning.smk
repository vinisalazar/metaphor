"""
Binning rules:
    - vamb: Bin contigs with vamb
    - maxbin2: Bin contigs with MaxBin2
"""


rule vamb:
    input:
        bamfiles=expand(
            "output/mapping/bam/{sample}",
            sample=sample_IDs,
        ),
        catalogue="{output}/binning/catalogue.fna.gz",
    output:
        outdir=directory("output/binning/vamb"),
    params:  # defaults in vamb's README
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
        rm -rf {output}

        vamb --outdir {output} \
             --fasta {input.catalogue} \
             --bamfiles {input.bamfiles} \
             -o {params.binsplit_sep} \
             --minfasta {params.minfasta} &> {log}
        """
