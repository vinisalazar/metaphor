"""
Annotation rules:

    - prodigal: gene prediction with Prodigal
    - diamond: protein annotation with Diamond
    - cog_parser: parse collated outputs with custom Python script
"""

from pathlib import Path


rule prodigal:
    input:
        contigs=get_contigs_input(),
    output:
        proteins=get_coassembly_or_sample_file("annotation", "prodigal", "proteins.faa"),
        genbank=get_coassembly_or_sample_file("annotation", "prodigal", "genbank.gbk"),
        genes=(get_coassembly_or_sample_file("annotation", "prodigal", "genes.fna")),
        # if config["prodigal"]["genes"]
        # else (),
        scores=(get_coassembly_or_sample_file("annotation", "prodigal", "scores.cds"))
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
    log:
        get_coassembly_benchmark_or_log("log", "annotation", "prodigal"),
    benchmark:
        get_coassembly_benchmark_or_log("benchmark", "annotation", "prodigal")
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


# TODO: refactor for coassembly
rule prokka:
    input:
        contigs=get_contigs_input(),
    output:
        outfile="output/annotation/prokka/{sample}/{sample}.faa",
    params:
        sample=lambda w: w.sample,
        outdir=lambda w, output: str(Path(output.outfile).parent),
        kingdom=config["prokka"]["kingdom"],
        args=config["prokka"]["args"],
    threads: round(workflow.cores * 0.25)
    log:
        "output/logs/annotation/prokka/{sample}.log",
    benchmark:
        "output/benchmarks/annotation/prokka/{sample}.txt"
    conda:
        "../envs/prokka.yaml"
    shell:
        """
        prokka --outdir {params.outdir}     \
               --kingdom {params.kingdom}   \
               --cpus {threads}             \
               --prefix {params.sample}     \
               {params.args}                \
               {input.contigs}          
        """


rule download_COG_database:
    output:
        get_database_outputs(),
    log:
        "output/logs/annotation/cog/download_COG_database.log",
    benchmark:
        "output/benchmarks/annotation/cog/download_COG_database.txt"
    conda:
        "../envs/bash.yaml"
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
    log:
        "output/logs/annotation/generate_COG_taxonmap.log",
    benchmark:
        "output/benchmarks/annotation/generate_COG_taxonmap.txt"
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/create_cog_taxonmap.py"


rule diamond_makedb:
    input:
        fname=get_cog_db_file("cog-20.fa.gz"),
        taxonmap=config["lineage_parser"]["taxonmap"],
        taxonnodes=config["lineage_parser"]["nodes"],
        taxonnames=config["lineage_parser"]["names"],
    output:
        fname=config["diamond"]["db"],
    params:
        extra=lambda w, input: f"--taxonmap {input.taxonmap} --taxonnames {input.taxonnames} --taxonnodes {input.taxonnodes}",
    log:
        "output/logs/annotation/diamond/diamond_makedb.log",
    threads: round(workflow.cores * 0.25)
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
        "../envs/bash.yaml"
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
        fname_fasta=get_coassembly_or_sample_file(
            "annotation", "prodigal", "proteins.faa"
        ),
        fname_db=config["diamond"]["db"],
    output:
        fname=get_diamond_output(),
    params:
        output_type=config["diamond"]["output_type"],
        output_format=config["diamond"]["output_format"],
        extra="--iterate --top 0",
    threads: round(workflow.cores * 0.75)
    log:
        get_coassembly_benchmark_or_log("log", "annotation", "diamond"),
    benchmark:
        get_coassembly_benchmark_or_log("benchmark", "annotation", "diamond")
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
        categories_out=get_coassembly_or_sample_file(
            "annotation", "cog", "categories.tsv"
        ),
        codes_out=get_coassembly_or_sample_file("annotation", "cog", "codes.tsv"),
        pathways_out=get_coassembly_or_sample_file("annotation", "cog", "pathways.tsv"),
    log:
        get_coassembly_benchmark_or_log("log", "annotation", "cog_functional_parser"),
    benchmark:
        get_coassembly_benchmark_or_log(
            "benchmark", "annotation", "cog_functional_parser"
        )
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/cog_functional_parser.py"


rule taxonomy_parser:
    input:
        dmnd_out=get_diamond_output(),
    output:
        tax_out=get_coassembly_or_sample_file("annotation", "cog", "tax.tsv"),
    log:
        get_coassembly_benchmark_or_log("log", "annotation", "cog_parser"),
    benchmark:
        get_coassembly_benchmark_or_log("benchmark", "annotation", "cog_parser")
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/taxonomy_parser.py"


rule concatenate_cog_functional:
    input:
        functional_counts=expand(
            "output/annotation/cog/{sample}/{sample}_{kind}.tsv",
            sample=sample_IDs,
            kind=functional_kinds,
        )
        if not config["coassembly"]
        else expand(
            get_coassembly_or_sample_file("annotation", "cog", "{kind}.tsv"),
            kind=functional_kinds,
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
        "../envs/bash.yaml"
    script:
        "../scripts/concatenate_cog_functional.py"


rule concatenate_taxonomies:
    input:
        files=lambda wildcards: expand(
            "output/annotation/cog/{sample}/{sample}_{rank}.tsv",
            sample=sample_IDs,
            rank=wildcards.rank,
        )
        if not config["coassembly"]
        else lambda wildcards: expand(
            get_coassembly_or_sample_file("annotation", "cog", "{rank}.tsv"),
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
        "../envs/bash.yaml"
    script:
        "../scripts/concatenate_taxonomies.py"


rule lineage_parser:
    input:
        tax_out=get_coassembly_or_sample_file("annotation", "cog", "tax.tsv"),
        rankedlineage=config["lineage_parser"]["rankedlineage"],
    output:
        # Class must be spelled with a 'k' to prevent conflicts with the Python keyword
        species=get_coassembly_or_sample_file("annotation", "cog", "species.tsv"),
        genus=get_coassembly_or_sample_file("annotation", "cog", "genus.tsv"),
        family=get_coassembly_or_sample_file("annotation", "cog", "family.tsv"),
        order=get_coassembly_or_sample_file("annotation", "cog", "order.tsv"),
        klass=get_coassembly_or_sample_file("annotation", "cog", "class.tsv"),
        phylum=get_coassembly_or_sample_file("annotation", "cog", "phylum.tsv"),
        kingdom=get_coassembly_or_sample_file("annotation", "cog", "kingdom.tsv"),
        domain=get_coassembly_or_sample_file("annotation", "cog", "domain.tsv"),
    log:
        get_coassembly_benchmark_or_log("log", "annotation", "lineage_parser"),
    benchmark:
        get_coassembly_benchmark_or_log("benchmark", "annotation", "lineage_parser")
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/lineage_parser.py"


rule plot_cog_functional:
    input:
        categories_file="output/annotation/cog/tables/COG_categories_relative.tsv",
    output:
        categories_plt=report(
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
        "../envs/bash.yaml"
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
    log:
        "output/logs/annotation/plot_taxonomies_{rank}.log",
    benchmark:
        "output/benchmarks/annotation/plot_taxonomies_{rank}.txt"
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/plot_taxonomies.py"