"""
annotation.smk

    Generates annotation data from contigs. First, coding sequences are predicted from contigs using Prodigal.
    Then, annotation (both taxonomic and functional) are done with the NCBI COG database
    (https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/). A Diamond database is built from the reference amino acid sequences
    and taxonomy is added with the NCBI Taxonomy database. It is possible to use a custom database for taxonomic annotation
    as long as a taxonmap (mapping of sequence accession number to TaxIDs) is provided.

Annotation rules:

    - prodigal: gene prediction with Prodigal
    - prokka: annotate bins with Prokka (WIP, may be deprecated)
    - download_COG_database: download COG database from NCBI URL
    - generate_COG_taxonmap: generate taxon map for Diamond db with NCBI COG and Taxonomy DBs
    - download_taxonomy_database: download NCBI taxonomy database
    - diamond: protein annotation with Diamond
    - cog_functional_parser: parses output of Diamond to get COG functional counts
    - taxonomy_parser: parses staxids columns from Diamond output
    - concatenate_cog_functional: concatenates cog_functional_parser output of all samples
    - concatenate_taxonomies: concatenates taxonomy_parser output of all samples
    - lineage_parser: gets complete lineage for each TaxID
    - plot_cog_functional: plots COG_functional counts
    - plot_cog_taxonomy: plots COG_taxonomy counts
"""

from pathlib import Path


rule prodigal:
    input:
        contigs=get_contigs_input(),
    output:
        proteins=get_group_or_sample_file("annotation", "prodigal", "proteins.faa"),
        genbank=get_group_or_sample_file("annotation", "prodigal", "genbank.gbk"),
        genes=(get_group_or_sample_file("annotation", "prodigal", "genes.fna")),
        # if config["prodigal"]["genes"]
        # else (),
        scores=(get_group_or_sample_file("annotation", "prodigal", "scores.cds"))
        if config["prodigal"]["scores"]
        else [],
    params:
        mode=config["prodigal"]["mode"],
        genes=lambda w, output: f"-d {output.genes}",
        # if config["prodigal"]["genes"]
        # else "",
        scores=lambda w, output: f"-s {output.scores}"
        if config["prodigal"]["scores"]
        else "",
        quiet="-q" if config["prodigal"]["quiet"] else "",
    wildcard_constraints:
        group="|".join(group_names),
    resources:
        mem_mb=get_max_mb(),
    log:
        get_group_benchmark_or_log("log", "annotation", "prodigal"),
    benchmark:
        get_group_benchmark_or_log("benchmark", "annotation", "prodigal")
    conda:
        "../envs/prodigal.yaml"
    shell:
        """
        prodigal {params.quiet}         \
                 -p {params.mode}       \
                 -i {input}             \
                 -a {output.proteins}   \
                 -o {output.genbank}    \
                 {params.genes}         \
                 {params.scores} &> {log}
        """


rule prokka:
    input:
        genome_bin="output/binning/DAS_tool/{binning_group}/DAS_tool_DASTool_bins/{bin}.fa",
    output:
        outfile="output/annotation/prokka/{binning_group}/{bin}/{bin}.fna",
    params:
        outdir=lambda w, output: str(Path(output.outfile).parent),
        kingdom=config["prokka"]["kingdom"],
        args=config["prokka"]["args"],
    wildcard_constraints:
        binning_group="|".join(binning_group_names),
    threads: get_threads_per_task_size("small")
    resources:
        mem_mb=get_max_mb(),
    log:
        "output/logs/annotation/prokka/{binning_group}/{bin}.log",
    benchmark:
        "output/benchmarks/annotation/prokka/{binning_group}/{bin}.txt"
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        prokka --outdir {params.outdir}     \
               --kingdom {params.kingdom}   \
               --cpus {threads}             \
               --prefix {wildcards.bin}     \
               {params.args}                \
               {input.genome_bin}          
        """


rule download_COG_database:
    output:
        get_database_outputs(),
    log:
        "output/logs/annotation/cog/download_COG_database.log",
    benchmark:
        "output/benchmarks/annotation/cog/download_COG_database.txt"
    conda:
        "../envs/utils.yaml"
    shell:
        """
        for file in {output}; do
            wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/$(basename $file) -O $file 2>> {log};
        done
        """


rule generate_COG_taxonmap:
    input:
        cog_csv=get_cog_db_file("cog-20.cog.csv"),
        org_csv=get_cog_db_file("cog-20.org.csv"),
    output:
        taxonmap=get_cog_db_file("cog-20.taxonmap.tsv"),
    resources:
        mem_mb=get_max_mb(),
    log:
        "output/logs/annotation/generate_COG_taxonmap.log",
    benchmark:
        "output/benchmarks/annotation/generate_COG_taxonmap.txt"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/create_cog_taxonmap.py"


rule diamond_makedb:
    input:
        fname=config["diamond"]["db_source"],
        taxonmap=config["lineage_parser"]["taxonmap"],
        taxonnodes=config["lineage_parser"]["nodes"],
        taxonnames=config["lineage_parser"]["names"],
    output:
        fname=config["diamond"]["db"],
    params:
        extra=lambda w, input: f"--taxonmap {input.taxonmap} --taxonnames {input.taxonnames} --taxonnodes {input.taxonnodes}",
    log:
        "output/logs/annotation/diamond/diamond_makedb.log",
    threads: get_threads_per_task_size("small")
    resources:
        mem_mb=get_max_mb(),
    wrapper:
        get_wrapper("diamond/makedb")


rule download_taxonomy_database:
    output:
        # Replace the .test/ directory with data/ directory if it is set like so
        rankedlineage=config["lineage_parser"]["rankedlineage"],
        names=config["lineage_parser"]["names"],
        nodes=config["lineage_parser"]["nodes"],
    params:
        download_url=config["lineage_parser"]["download_url"],
        output_dir=lambda w, output: str(Path(output.rankedlineage).parent),
    log:
        "output/logs/annotation/cog/download_taxonomy_database.log",
    conda:
        "../envs/utils.yaml"
    shell:
        """
        mkdir -p {params.output_dir} 2>> {log}

        DST={params.output_dir}/$(basename {params.download_url})

        # No verbose (nv) option to avoid progress bar
        wget -nv {params.download_url} -O $DST 2>> {log}

        tar zxvf $DST -C {params.output_dir} 2>> {log}

        ls {params.output_dir} 2>> {log}
        """


rule diamond:
    input:
        fname_fasta=get_group_or_sample_file("annotation", "prodigal", "proteins.faa"),
        fname_db=config["diamond"]["db"],
    output:
        fname=get_diamond_output(),
    params:
        output_type=config["diamond"]["output_type"],
        output_format=config["diamond"]["output_format"],
        extra="--iterate --top 0",
    wildcard_constraints:
        group="|".join(group_names),
    threads: get_threads_per_task_size("big")
    resources:
        mem_mb=get_mb_per_cores,
    log:
        get_group_benchmark_or_log("log", "annotation", "diamond"),
    benchmark:
        get_group_benchmark_or_log("benchmark", "annotation", "diamond")
    conda:
        "../envs/diamond.yaml"
    shell:
        """
        echo {params.output_format} | sed -e 's/ /\t/g' > {output.fname}
        {{ diamond blastp -q {input.fname_fasta}                \
                   -p {threads}                                 \
                   -d {input.fname_db}                          \
                   -f {params.output_type}                      \
                   {params.output_format}                       \
                   {params.extra}                               \
                   >> {output.fname} ; }} &> {log}
        """


rule cog_functional_parser:
    input:
        dmnd_out=get_diamond_output(),
        cog_csv=get_cog_db_file("cog-20.cog.csv"),
        def_tab=get_cog_db_file("cog-20.def.tab"),
        fun_tab=get_cog_db_file("fun-20.tab"),
    output:
        categories_out=get_group_or_sample_file("annotation", "cog", "categories.tsv"),
        codes_out=get_group_or_sample_file("annotation", "cog", "codes.tsv"),
        pathways_out=get_group_or_sample_file("annotation", "cog", "pathways.tsv"),
    log:
        get_group_benchmark_or_log("log", "annotation", "cog_functional_parser"),
    benchmark:
        get_group_benchmark_or_log("benchmark", "annotation", "cog_functional_parser")
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/cog_functional_parser.py"


rule taxonomy_parser:
    input:
        dmnd_out=get_diamond_output(),
    output:
        tax_out=get_group_or_sample_file("annotation", "cog", "tax.tsv"),
    log:
        get_group_benchmark_or_log("log", "annotation", "cog_parser"),
    benchmark:
        get_group_benchmark_or_log("benchmark", "annotation", "cog_parser")
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/taxonomy_parser.py"


rule concatenate_cog_functional:
    input:
        functional_counts=lambda wildcards: expand(
            "output/annotation/cog/{group}/{group}_{kind}.tsv",
            group=group_names,
            kind=wildcards.kind,
        ),
    output:
        functional_absolute_counts="output/annotation/cog/tables/COG_{kind}_absolute.tsv",
        functional_relative_counts="output/annotation/cog/tables/COG_{kind}_relative.tsv",
    wildcard_constraints:
        kind="|".join(functional_kinds),
    log:
        "output/logs/annotation/concatenate_cog_{kind}.log",
    benchmark:
        "output/benchmarks/annotation/concatenate_cog_{kind}.txt"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/concatenate_cog_functional.py"


rule concatenate_taxonomies:
    input:
        files=lambda wildcards: expand(
            "output/annotation/cog/{group}/{group}_{rank}.tsv",
            group=group_names,
            rank=wildcards.rank,
        ),
    output:
        absolute_counts="output/annotation/cog/tables/COG_{rank}_absolute.tsv",
        relative_counts="output/annotation/cog/tables/COG_{rank}_relative.tsv",
    wildcard_constraints:
        rank="|".join(
            ranks
            + [
                "tax",
            ]
        ),
    log:
        "output/logs/annotation/concatenate_cog_{rank}.log",
    benchmark:
        "output/benchmarks/annotation/concatenate_cog_{rank}.txt"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/concatenate_taxonomies.py"


rule lineage_parser:
    input:
        tax_out=get_group_or_sample_file("annotation", "cog", "tax.tsv"),
        rankedlineage=config["lineage_parser"]["rankedlineage"],
    output:
        # Class must be spelled with a 'k' to prevent conflicts with the Python keyword
        species=get_group_or_sample_file("annotation", "cog", "species.tsv"),
        genus=get_group_or_sample_file("annotation", "cog", "genus.tsv"),
        family=get_group_or_sample_file("annotation", "cog", "family.tsv"),
        order=get_group_or_sample_file("annotation", "cog", "order.tsv"),
        klass=get_group_or_sample_file("annotation", "cog", "class.tsv"),
        phylum=get_group_or_sample_file("annotation", "cog", "phylum.tsv"),
        domain=get_group_or_sample_file("annotation", "cog", "domain.tsv"),
    resources:
        mem_mb=get_max_mb(0.5),
    log:
        get_group_benchmark_or_log("log", "annotation", "lineage_parser"),
    benchmark:
        get_group_benchmark_or_log("benchmark", "annotation", "lineage_parser")
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/lineage_parser.py"


rule plot_cog_functional:
    input:
        categories_file="output/annotation/cog/tables/COG_categories_relative.tsv",
    output:
        categories_plot=report(
            get_cog_functional_plot_outputs(),
            category="Annotation",
        ),
    params:
        filter_categories=config["plot_cog_functional"]["filter_categories"],
        categories_cutoff=config["plot_cog_functional"]["categories_cutoff"],
    log:
        "output/logs/annotation/plot_cog_functional.log",
    benchmark:
        "output/benchmarks/annotation/plot_cog_functional.txt"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/plot_cog_functional.py"


rule plot_cog_taxonomy:
    input:
        taxonomy_relative_counts="output/annotation/cog/tables/COG_{rank}_relative.tsv",
    output:
        taxonomy_barplot=report(
            "output/annotation/cog/plots/COG_{rank}_relative.png",
            category="Annotation",
        ),
    params:
        tax_cutoff=config["plot_taxonomies"]["tax_cutoff"],
        colormap=config["plot_taxonomies"]["colormap"],
    log:
        "output/logs/annotation/plot_taxonomies_{rank}.log",
    benchmark:
        "output/benchmarks/annotation/plot_taxonomies_{rank}.txt"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/plot_taxonomies.py"
