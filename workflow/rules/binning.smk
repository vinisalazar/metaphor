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
        outdir=directory("output/binning/metabat2"),
    params:
        minContig=2500,
        seed=config["metabat2"]["seed"],
        outfile=lambda w, output: output.outdir + "/" + config["metabat2"]["preffix"]
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


rule concoct:
    input:
        catalogue="output/mapping/catalogue.fna.gz",
        bams=expand("output/mapping/bam/{sample}.sorted.bam", sample=sample_IDs)
    output:
        outdir=directory("output/binning/concoct/"),
    params:
        contig_size=10000,
        bed=lambda w, output: output.outdir + "/contigs.bed",
        contigs=lambda w, output: output.outdir + "/contigs.fa",
        coverage_table=lambda w, output: output.outdir + "/coverage_table.tsv",
        fasta_bins=lambda w, output: output.outdir + "/fasta_bins",
        clustering_gt=lambda w, output: output.outdir + "/clustering_gt1000.csv",
        clustering_merged=lambda w, output: output.outdir + "/clustering_merged.csv",
        uncompressed_catalogue=lambda w, input: input.catalogue.replace(".gz", "")
    threads: round(workflow.cores * 0.75)
    log:
        "output/logs/binning/concoct.log"
    benchmark:
        "output/benchmarks/binning/concoct.log"
    conda:
        "../envs/concoct.yaml"
    shell:
        """ 
        rm -rf {output.outdir}
        mkdir {output.outdir} 

        {{ pigz -d -p {threads} -k {input.catalogue} ; }} 2>> {log}

        {{ cut_up_fasta.py {params.uncompressed_catalogue}          \
                        -c {params.contig_size}                     \
                        -o 0                                        \
                        -b {params.bed}                             \
                        --merge_last                                \
                        > {params.contigs}  ; }} 2>> {log}

        {{ concoct_coverage_table.py {params.bed}                   \
                                    {input.bams}                    \
                                    > {params.coverage_table} ; }} 2>> {log}

        {{ concoct --composition_file {params.contigs}              \
                  --coverage_file {params.coverage_table}           \
                  -b {output.outdir}                                \
                  -t {threads}  ; }} 2>> {log}

        {{ merge_cutup_clustering.py {params.clustering_gt}          \
            > {params.clustering_merged}  ; }} 2>> {log}

        mkdir {params.fasta_bins}

        {{ extract_fasta_bins.py {params.uncompressed_catalogue}    \
                              {params.clustering_merged}            \
                              --output_path {params.fasta_bins} ; }} 2>> {log}
        
        rm {params.uncompressed_catalogue}
        mv {output.outdir}/../concoct_* {output.outdir}
        """
