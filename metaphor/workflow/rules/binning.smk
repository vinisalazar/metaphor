"""
binning.smk

    Cluster contigs into genome bins and refine said bins.

Binning rules:
    - vamb: Bin contigs with vamb
    - MetaBAT2: Bin contigs with MetaBat2
    - concoct: Bin contigs with CONCOCT
    - DAS_tool: Refine bins with DAS_tool
"""
from pathlib import Path


rule vamb:
    input:
        bam_contig_depths="output/mapping/{group}/bam_contig_depths.txt",
        catalogue="output/mapping/{group}/catalogue.fna.gz",
    output:
        clusters=get_vamb_output(),
        scaffolds2bin="output/binning/DAS_tool/{group}/vamb_scaffolds2bin.tsv",
    params:  # defaults in vamb's README
        outdir=lambda w, output: get_parent(output.clusters),
        binsplit_sep="C",
        minfasta=config["vamb"]["minfasta"],
        batchsize=256,
    threads: round(workflow.cores * config["cores_per_big_task"])
    resources:
        mem_mb=get_mb_per_cores,
    log:
        "output/logs/binning/{group}/vamb.log",
    benchmark:
        "output/benchmarks/binning/{group}/vamb.txt"
    conda:
        "../envs/vamb.yaml"
    shell:
        """
        rm -rf {params.outdir}

        vamb --outdir {params.outdir}           \
             --fasta {input.catalogue}          \
             --jgi {input.bam_contig_depths}    \
             -p {threads}                       \
             -o {params.binsplit_sep}           \
             -t {params.batchsize}              \
             --minfasta {params.minfasta} &> {log}

        {{ awk -v OFS='\t' '{{ print $2, $1 }}' {output.clusters} |  \
        sed "s/$(echo '\t')/$(echo '\t')vamb./g" >          \
        {output.scaffolds2bin} ; }} >> {log} 2>&1
        """


rule metabat2:
    input:
        contigs="output/mapping/{group}/catalogue.fna.gz",
        depths="output/mapping/{group}/bam_contig_depths.txt",
    output:
        outdir=directory("output/binning/metabat2/{group}/"),
        scaffolds2bin="output/binning/DAS_tool/{group}/metabat2_scaffolds2bin.tsv",
    params:
        minContig=2500,
        seed=config["metabat2"]["seed"],
        outfile=lambda w, output: str(
            Path(output.outdir).joinpath(config["metabat2"]["preffix"])
        ),
    threads: round(workflow.cores * config["cores_per_big_task"])
    resources:
        mem_mb=get_mb_per_cores,
    log:
        "output/logs/binning/{group}/metabat2.log",
    benchmark:
        "output/benchmarks/binning/{group}/metabat2.txt"
    conda:
        "../envs/metabat2.yaml"
    shell:
        """
        rm -rf {output.outdir} && mkdir {output.outdir}

        metabat2 -i {input.contigs}             \
                 -a {input.depths}              \
                 -m {params.minContig}          \
                 -t {threads}                   \
                 --seed {params.seed}           \
                 --saveCls                      \
                 -o {params.outfile} &> {log}

        sed "s/$(echo '\t')/$(echo '\t')metabat2./g" {params.outfile} > {output.scaffolds2bin}
        """


rule concoct:
    input:
        catalogue="output/mapping/{group}/catalogue.fna",
        bams=expand("output/mapping/bam/{{group}}/{sample}.sorted.bam", sample=sample_IDs),
        bais=expand("output/mapping/bam/{{group}}/{sample}.sorted.bam.bai", sample=sample_IDs),
    output:
        outdir=directory("output/binning/concoct/{group}/"),
        scaffolds2bin="output/binning/DAS_tool/{group}/concoct_scaffolds2bin.tsv",
    params:
        contig_size=10000,
        bed=lambda w, output: str(Path(output.outdir).joinpath("contigs.bed")),
        contigs=lambda w, output: str(Path(output.outdir).joinpath("contigs.fa")),
        coverage_table=lambda w, output: str(
            Path(output.outdir).joinpath("coverage_table.tsv")
        ),
        fasta_bins=lambda w, output: str(Path(output.outdir).joinpath("fasta_bins")),
        clustering_gt=lambda w, output: str(
            Path(output.outdir).joinpath("clustering_gt1000.csv")
        ),
        clustering_merged=lambda w, output: str(
            Path(output.outdir).joinpath("clustering_merged.csv")
        ),
    threads: round(workflow.cores * config["cores_per_big_task"])
    resources:
        mem_mb=get_mb_per_cores,
    log:
        "output/logs/binning/{group}/concoct.log",
    benchmark:
        "output/benchmarks/binning/{group}/concoct.txt"
    conda:
        "../envs/concoct.yaml"
    shell:
        """ 
        rm -rf {output.outdir}
        mkdir {output.outdir} 

        {{ cut_up_fasta.py {input.catalogue}                        \
                           -c {params.contig_size}                  \
                           -o 0                                     \
                           -b {params.bed}                          \
                           --merge_last                             \
                           > {params.contigs}  ; }} 2>> {log}

        {{ concoct_coverage_table.py {params.bed}                   \
                                     {input.bams}                   \
                                     > {params.coverage_table} ; }} 2>> {log}

        {{ concoct --composition_file {params.contigs}              \
                   --coverage_file {params.coverage_table}          \
                   -b {output.outdir}                               \
                   -t {threads}  ; }} 2>> {log}

        {{ merge_cutup_clustering.py {params.clustering_gt}         \
                                     > {params.clustering_merged}  ; }} 2>> {log}

        mkdir {params.fasta_bins}

        {{ extract_fasta_bins.py {output.uncompressed_catalogue}    \
                                 {params.clustering_merged}         \
                                 --output_path {params.fasta_bins} ; }} 2>> {log}

        sed "s/,/$(echo '\t')concoct./g" {params.clustering_merged} | tail -n +2 > {output.scaffolds2bin}
        """


rule DAS_tool:
    """
    Refine bins assembled with one or more binners.
    """
    input:
        contigs="output/mapping/{group}/catalogue.fna",
        scaffolds2bin=get_DAS_tool_input(),
    output:
        proteins="output/binning/DAS_tool/{group}/DAS_tool_proteins.faa",
    params:
        fmt_scaffolds2bin=lambda w, input: ",".join(input.scaffolds2bin),
        binners=",".join(binners),
        outpreffix=lambda w, output: str(
            Path(output.proteins).parent.joinpath("DAS_tool")
        ),
        score_threshold=config["das_tool"]["score_threshold"],
    threads: round(workflow.cores * config["cores_per_big_task"])
    resources:
        mem_mb=get_mb_per_cores,
    log:
        "output/logs/binning/{group}/DAS_tool.log",
    benchmark:
        "output/benchmarks/binning/{group}/DAS_tool.txt"
    conda:
        "../envs/das_tool.yaml"
    shell:
        """
        DAS_Tool -i {params.fmt_scaffolds2bin}                      \
                 -l {params.binners}                                \
                 -c {input.contigs}                                 \
                 -o {params.outpreffix}                             \
                 --score_threshold {params.score_threshold}         \
                 --search_engine diamond                            \
                 --write_bins                                       \
                 --threads {threads} &> {log}

        rm -rf {input.contigs}
        """
