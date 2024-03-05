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
        contigs="output/mapping/{group}/{group}_contigs_catalogue.fna",
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
        group="|".join(binning_group_names),
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
        genome_bin="output/binning/DAS_tool/{binning_group}/{binning_group}_DASTool_bins/{bin}.fa",
        bin_evals="output/binning/DAS_tool/{binning_group}/{binning_group}_DASTool_summary.tsv",
    output:
        outfile="output/annotation/prokka/{binning_group}/{bin}/{bin}.fna",
    params:
        outdir=lambda w, output: str(Path(output.outfile).parent),
        args=config["prokka"]["args"] + " --centre X --compliant",
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
        # Get kingdom from bin eval file
        bin_clean=$(echo {wildcards.bin} | sed 's/_sub$//g')  # Remove '_sub' from corrected bins
        kingdom=$(grep $bin_clean {input.bin_evals} | cut -f 5)
        kingdom=$(echo $kingdom | cut -f 1 -d ' ')
        kingdom=$(echo $kingdom | head -c 1 | tr '[a-z]' '[A-Z]'; echo $kingdom | tail -c +2)

        prokka --outdir {params.outdir}     \
               --kingdom $kingdom           \
               --cpus {threads}             \
               --prefix {wildcards.bin}     \
               {params.args}                \
               {input.genome_bin} &> {log}
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
        "../envs/python-utils.yaml"
    script:
        "../scripts/create_cog_taxonmap.py"


rule diamond_makedb:
    input:
        fname=add_data_dir(config["diamond"]["db_source"]),
        taxonmap=add_data_dir(config["lineage_parser"]["taxonmap"]),
        taxonnodes=add_data_dir(config["lineage_parser"]["nodes"]),
        taxonnames=add_data_dir(config["lineage_parser"]["names"]),
    output:
        fname=add_data_dir(config["diamond"]["db"]),
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
        rankedlineage=add_data_dir(config["lineage_parser"]["rankedlineage"]),
        names=add_data_dir(config["lineage_parser"]["names"]),
        nodes=add_data_dir(config["lineage_parser"]["nodes"]),
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
        fname_db=add_data_dir(config["diamond"]["db"]),
    output:
        fname=get_diamond_output(),
    params:
        output_type=config["diamond"]["output_type"],
        output_format=config["diamond"]["output_format"],
        extra="--iterate --top 0",
    wildcard_constraints:
        group="|".join(binning_group_names),
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
        coverage_depths="output/mapping/{group}/bam_contigs_depths.txt",
    output:
        categories_out_absolute=get_group_or_sample_file(
            "annotation", "cog", "categories_absolute.tsv"
        ),
        categories_out_relative=get_group_or_sample_file(
            "annotation", "cog", "categories_relative.tsv"
        ),
        codes_out_absolute=get_group_or_sample_file(
            "annotation", "cog", "codes_absolute.tsv"
        ),
        codes_out_relative=get_group_or_sample_file(
            "annotation", "cog", "codes_relative.tsv"
        ),
        pathways_out_absolute=get_group_or_sample_file(
            "annotation", "cog", "pathways_absolute.tsv"
        ),
        pathways_out_relative=get_group_or_sample_file(
            "annotation", "cog", "pathways_relative.tsv"
        ),
    wildcard_constraints:
        group="|".join(binning_group_names),
    log:
        get_group_benchmark_or_log("log", "annotation", "cog_functional_parser"),
    benchmark:
        get_group_benchmark_or_log("benchmark", "annotation", "cog_functional_parser")
    conda:
        "../envs/python-utils.yaml"
    script:
        "../scripts/cog_functional_parser.py"


rule taxonomy_parser:
    input:
        dmnd_out=get_diamond_output(),
        coverage_depths="output/mapping/{group}/bam_contigs_depths.txt",
    output:
        tax_out_absolute=get_group_or_sample_file(
            "annotation", "cog", "tax_absolute.tsv"
        ),
        tax_out_relative=get_group_or_sample_file(
            "annotation", "cog", "tax_relative.tsv"
        ),
    log:
        get_group_benchmark_or_log("log", "annotation", "cog_parser"),
    wildcard_constraints:
        group="|".join(binning_group_names),
    benchmark:
        get_group_benchmark_or_log("benchmark", "annotation", "cog_parser")
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/taxonomy_parser.py"


rule concatenate_cog_functional:
    input:
        functional_counts=lambda wildcards: expand(
            "output/annotation/cog/{group}/{group}_{kind}_{count_type}.tsv",
            group=binning_group_names,
            kind=wildcards.kind,
            count_type=wildcards.count_type,
        ),
    output:
        concatenated_functional_counts="output/annotation/cog/tables/concatenated_{kind}_{count_type}.tsv",
    wildcard_constraints:
        kind="|".join(functional_kinds),
        count_type="absolute|relative",
    log:
        "output/logs/annotation/concatenate_{kind}_{count_type}.log",
    benchmark:
        "output/benchmarks/annotation/concatenate_{kind}_{count_type}.txt"
    conda:
        "../envs/python-utils.yaml"
    script:
        "../scripts/concatenate_cog_functional.py"


rule concatenate_taxonomies:
    input:
        taxonomy_counts=lambda wildcards: expand(
            "output/annotation/cog/{group}/{group}_{rank}_{count_type}.tsv",
            group=binning_group_names,
            rank=wildcards.rank,
            count_type=wildcards.count_type,
        ),
    output:
        concatenated_taxonomy_counts="output/annotation/cog/tables/concatenated_{rank}_{count_type}.tsv",
    wildcard_constraints:
        rank="|".join(
            ranks
            + [
                "tax",
            ]
        ),
        count_type="absolute|relative",
    log:
        "output/logs/annotation/concatenate_cog_{rank}_{count_type}.log",
    benchmark:
        "output/benchmarks/annotation/concatenate_cog_{rank}_{count_type}.txt"
    conda:
        "../envs/python-utils.yaml"
    script:
        "../scripts/concatenate_taxonomies.py"


rule lineage_parser:
    input:
        tax_out=get_group_or_sample_file("annotation", "cog", "tax_{count_type}.tsv"),
        rankedlineage=add_data_dir(config["lineage_parser"]["rankedlineage"]),
    output:
        # Class must be spelled with a 'k' to prevent conflicts with the Python keyword
        species=get_group_or_sample_file(
            "annotation", "cog", "species_{count_type}.tsv"
        ),
        genus=get_group_or_sample_file("annotation", "cog", "genus_{count_type}.tsv"),
        family=get_group_or_sample_file("annotation", "cog", "family_{count_type}.tsv"),
        order=get_group_or_sample_file("annotation", "cog", "order_{count_type}.tsv"),
        klass=get_group_or_sample_file("annotation", "cog", "class_{count_type}.tsv"),
        phylum=get_group_or_sample_file("annotation", "cog", "phylum_{count_type}.tsv"),
        domain=get_group_or_sample_file("annotation", "cog", "domain_{count_type}.tsv"),
    resources:
        mem_mb=get_max_mb(0.5),
    wildcard_constraints:
        group="|".join(binning_group_names),
    log:
        get_group_benchmark_or_log("log", "annotation", "lineage_parser_{count_type}"),
    benchmark:
        get_group_benchmark_or_log(
            "benchmark", "annotation", "lineage_parser_{count_type}"
        )
    conda:
        "../envs/python-utils.yaml"
    script:
        "../scripts/lineage_parser.py"


rule plot_cog_functional:
    input:
        categories_file="output/annotation/cog/{group}/{group}_categories_relative.tsv",
    output:
        categories_plot=report(
            get_cog_functional_plot_output("{group}"),
            category="Annotation",
        ),
    params:
        filter_categories=config["plot_cog_functional"]["filter_categories"],
        categories_cutoff=config["plot_cog_functional"]["categories_cutoff"],
        white_background=not config["transparent_background"],
        dpi=config["dpi"],
    wildcard_constraints:
        group="|".join(binning_group_names),
    log:
        "output/logs/annotation/plot_cog_functional/{group}.log",
    benchmark:
        "output/benchmarks/annotation/plot_cog_functional/{group}.txt"
    conda:
        "../envs/python-utils.yaml"
    script:
        "../scripts/plot_cog_functional.py"


rule plot_cog_taxonomy:
    input:
        taxonomy_relative_counts="output/annotation/cog/{group}/{group}_{rank}_relative.tsv",
    output:
        taxonomy_barplot=report(
            "output/annotation/cog/{group}/plots/{group}_{rank}_relative.png",
            category="Annotation",
        ),
    params:
        tax_cutoff=config["plot_taxonomies"]["tax_cutoff"],
        colormap=config["plot_taxonomies"]["colormap"],
        white_background=not config["transparent_background"],
        dpi=config["dpi"],
        output_format=config["output_format"],
    log:
        "output/logs/annotation/plot_taxonomies_{rank}/{group}.log",
    benchmark:
        "output/benchmarks/annotation/plot_taxonomies_{rank}/{group}.txt"
    conda:
        "../envs/python-utils.yaml"
    script:
        "../scripts/plot_taxonomies.py"
