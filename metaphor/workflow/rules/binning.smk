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
        bam_contig_depths="output/mapping/{binning_group}/bam_contigs_depths.txt",
        catalogue="output/mapping/{binning_group}/{binning_group}_contigs_catalogue.fna.gz",
    output:
        clusters=get_vamb_output(),
        scaffolds2bin="output/binning/vamb/{binning_group}/vamb_scaffolds2bin.tsv",
    params:  # defaults in vamb's README
        outdir=lambda w, output: get_parent(output.clusters),
        binsplit_sep="C",
        minfasta=config["vamb"]["minfasta"],
        batchsize=config["vamb"]["batchsize"],
    threads: get_threads_per_task_size("big")
    resources:
        mem_mb=get_mb_per_cores,
    log:
        "output/logs/binning/vamb/{binning_group}.log",
    benchmark:
        "output/benchmarks/binning/vamb/{binning_group}.txt"
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
        contigs="output/mapping/{binning_group}/{binning_group}_contigs_catalogue.fna.gz",
        depths="output/mapping/{binning_group}/bam_contigs_depths.txt",
    output:
        outdir=directory("output/binning/metabat2/{binning_group}/"),
        scaffolds2bin="output/binning/metabat2/{binning_group}/metabat2_scaffolds2bin.tsv",
    params:
        minContig=2500,
        seed=config["metabat2"]["seed"],
        outfile=lambda w, output: str(
            Path(output.outdir).joinpath(config["metabat2"]["preffix"])
        ),
    threads: get_threads_per_task_size("big")
    resources:
        mem_mb=get_mb_per_cores,
    log:
        "output/logs/binning/metabat2/{binning_group}.log",
    benchmark:
        "output/benchmarks/binning/metabat2/{binning_group}.txt"
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
        catalogue="output/mapping/{binning_group}/{binning_group}_contigs_catalogue.fna",
        bams=lambda wildcards: expand(
            "output/mapping/bam/{{binning_group}}/{sample}-to-contigs.sorted.bam",
        sample=list(
            samples.query(f"binning_group == '{wildcards.binning_group}'")[
        "sample_name"
                ].unique()
            ),
        ),
        bais=lambda wildcards: expand(
            "output/mapping/bam/{{binning_group}}/{sample}-to-contigs.sorted.bam.bai",
        sample=list(
            samples.query(f"binning_group == '{wildcards.binning_group}'")[
        "sample_name"
                ].unique()
            ),
        ),
    output:
        outdir=directory("output/binning/concoct/{binning_group}/"),
        scaffolds2bin="output/binning/concoct/{binning_group}/concoct_scaffolds2bin.tsv",
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
    threads: get_threads_per_task_size("big")
    resources:
        mem_mb=get_mb_per_cores,
    log:
        "output/logs/binning/concoct/{binning_group}.log",
    benchmark:
        "output/benchmarks/binning/concoct/{binning_group}.txt"
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

        {{ extract_fasta_bins.py {input.catalogue}                  \
                                 {params.clustering_merged}         \
                                 --output_path {params.fasta_bins} ; }} 2>> {log}

        sed "s/,/$(echo '\t')concoct./g" {params.clustering_merged} | tail -n +2 > {output.scaffolds2bin}
        """


rule DAS_tool:
    """
    Refine bins assembled with one or more binners.
    """
    input:
        contigs="output/mapping/{binning_group}/{binning_group}_contigs_catalogue.fna",
        scaffolds2bin=get_DAS_tool_input(),
        proteins="output/annotation/prodigal/{binning_group}/{binning_group}_proteins.faa"
    output:
        bin_evals="output/binning/DAS_tool/{binning_group}/{binning_group}_DASTool_summary.tsv",
    params:
        fmt_scaffolds2bin=lambda w, input: ",".join(input.scaffolds2bin),
        binners=",".join(binners),
        outpreffix=lambda w, output: str(
            Path(output.bin_evals).parent.joinpath(w.binning_group)
        ),
        score_threshold=config["das_tool"]["score_threshold"],
        extra=""
    threads: get_threads_per_task_size("big")
    resources:
        mem_mb=get_mb_per_cores,
    log:
        "output/logs/binning/DAS_tool/{binning_group}.log",
    benchmark:
        "output/benchmarks/binning/DAS_tool/{binning_group}.txt"
    conda:
        "../envs/das_tool.yaml"
    shell:
        """
        DAS_Tool -i {params.fmt_scaffolds2bin}                      \
                 -l {params.binners}                                \
                 -c {input.contigs}                                 \
                 -o {params.outpreffix}                             \
                 -p {input.proteins}                                \
                 --score_threshold {params.score_threshold}         \
                 --search_engine diamond                            \
                 --write_bins                                       \
                 {params.extra}                                     \
                 --threads {threads} &> {log}
        """


rule bins_report:
    """
    Plots binning metrics generated by DAS Tool.
    """
    input:
        bins_eval="output/binning/DAS_tool/{binning_group}/{binning_group}_DASTool_summary.tsv",
    output:
        **get_binning_report_output("{binning_group}"),
    params:
        score_threshold=config["das_tool"]["score_threshold"],
        binning_group=lambda w: w.binning_group,
    log:
        "output/logs/binning/plots/{binning_group}.log",
    benchmark:
        "output/benchmarks/binning/plots/{binning_group}.txt"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/bins_report.py"
