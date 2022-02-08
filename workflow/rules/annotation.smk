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
        proteins=get_coassembly_or_sample_file("annotation", "prodigal", "proteins.faa")
        genbank=get_coassembly_or_sample_file("annotation", "prodigal", "genbank.gbk")
        genes=(
            get_coassembly_or_sample_file("annotation", "prodigal", "genes.fna")
        ),
        # if config["prodigal"]["genes"]
        # else (),
        scores=(
            get_coassembly_or_sample_file("annotation", "prodigal", "scores.cds")
        )
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
        "output/logs/annotation/download_COG_database.log",
    benchmark:
        "output/benchmarks/annotation/download_COG_database.txt"
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
        taxonmap=get_cog_db_file("cog-20.taxonmap.tsv"),
        taxonnodes=config["lineage_parser"]["names"],
        taxonnames=config["lineage_parser"]["nodes"],
    output:
        fname=config["diamond"]["db"],
    params:
        extra=lambda w, input: "--taxonmap {input.taxonmap} --taxonnames {input.taxonnames} --taxonnodes {input.taxonnodes}"
    log:
        "output/logs/annotation/diamond/diamond_makedb.log",
    threads: round(workflow.cores * 0.25)
    wrapper:
        get_wrapper("diamond/makedb")


rule download_taxonomy_database:
    output:
        # Replace the .test/ directory with data/ directory if it is set like so
        rankedlineage=config["lineage_parser"]["lineage_parser"],
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

        wget {params.download_url} -O $DST 2>> {log}

        tar zxvf $DST -C {params.output_dir} 2>> {log}

        ls {params.output_dir} 2>> {log}
        """


rule diamond:
    input:
        fname_fasta="output/annotation/prodigal/{sample}/{sample}_proteins.faa"
        if not config["coassembly"]
        else "output/annotation/prodigal/coassembly_proteins.faa",
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


rule cog_parser:
    input:
        dmnd_out=get_diamond_output(),
        cog_csv=get_cog_db_file("cog-20.cog.csv"),
        def_tab=get_cog_db_file("cog-20.def.tab"),
        fun_tab=get_cog_db_file("fun-20.tab"),
    output:
        categories_out="output/annotation/cog/{sample}/{sample}_categories.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_categories.tsv",
        codes_out="output/annotation/cog/{sample}/{sample}_codes.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_codes.tsv",
        tax_out="output/annotation/cog/{sample}/{sample}_tax.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_tax.tsv",
        pathways_out="output/annotation/cog/{sample}/{sample}_pathways.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_pathways.tsv",
    log:
        get_coassembly_benchmark_or_log("log", "annotation", "cog_parser"),
    benchmark:
        get_coassembly_benchmark_or_log("benchmark", "annotation", "cog_parser")
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/cog_parser.py"


rule concatenate_cog:
    input:
        categories=expand(
            "output/annotation/cog/{sample}/{sample}_categories.tsv", sample=sample_IDs
        ),
        codes=expand(
            "output/annotation/cog/{sample}/{sample}_codes.tsv", sample=sample_IDs
        ),
        taxs=expand(
            "output/annotation/cog/{sample}/{sample}_tax.tsv", sample=sample_IDs
        ),
        pathways=expand(
            "output/annotation/cog/{sample}/{sample}_pathways.tsv", sample=sample_IDs
        ),
        species=expand(
            "output/annotation/cog/{sample}/{sample}_species.tsv", sample=sample_IDs
        ),
        genus=expand(
            "output/annotation/cog/{sample}/{sample}_genus.tsv", sample=sample_IDs
        ),
        family=expand(
            "output/annotation/cog/{sample}/{sample}_family.tsv", sample=sample_IDs
        ),
        order=expand(
            "output/annotation/cog/{sample}/{sample}_order.tsv", sample=sample_IDs
        ),
        klass=expand(
            "output/annotation/cog/{sample}/{sample}_class.tsv", sample=sample_IDs
        ),
        phylum=expand(
            "output/annotation/cog/{sample}/{sample}_phylum.tsv", sample=sample_IDs
        ),
        kingdom=expand(
            "output/annotation/cog/{sample}/{sample}_kingdom.tsv", sample=sample_IDs
        ),
        domain=expand(
            "output/annotation/cog/{sample}/{sample}_domain.tsv", sample=sample_IDs
        ),
    output:
        # Unfortunately this ugly block of code is required due to standardization of argument parsing across the workflow
        categories_absolute="output/annotation/cog/COG_categories_absolute.tsv",
        categories_relative="output/annotation/cog/COG_categories_relative.tsv",
        codes_absolute="output/annotation/cog/COG_codes_absolute.tsv",
        codes_relative="output/annotation/cog/COG_codes_relative.tsv",
        taxs_absolute="output/annotation/cog/COG_taxs_absolute.tsv",
        taxs_relative="output/annotation/cog/COG_taxs_relative.tsv",
        pathways_absolute="output/annotation/cog/COG_pathways_absolute.tsv",
        pathways_relative="output/annotation/cog/COG_pathways_relative.tsv",
        species_absolute="output/annotation/cog/COG_species_absolute.tsv",
        species_relative="output/annotation/cog/COG_species_relative.tsv",
        genus_absolute="output/annotation/cog/COG_genus_absolute.tsv",
        genus_relative="output/annotation/cog/COG_genus_relative.tsv",
        family_absolute="output/annotation/cog/COG_family_absolute.tsv",
        family_relative="output/annotation/cog/COG_family_relative.tsv",
        order_absolute="output/annotation/cog/COG_order_absolute.tsv",
        order_relative="output/annotation/cog/COG_order_relative.tsv",
        klass_absolute="output/annotation/cog/COG_class_absolute.tsv",
        klass_relative="output/annotation/cog/COG_class_relative.tsv",
        phylum_absolute="output/annotation/cog/COG_phylum_absolute.tsv",
        phylum_relative="output/annotation/cog/COG_phylum_relative.tsv",
        kingdom_absolute="output/annotation/cog/COG_kingdom_absolute.tsv",
        kingdom_relative="output/annotation/cog/COG_kingdom_relative.tsv",
        domain_absolute="output/annotation/cog/COG_domain_absolute.tsv",
        domain_relative="output/annotation/cog/COG_domain_relative.tsv",
    log:
        "output/logs/annotation/concatenate_cog.log",
    benchmark:
        "output/benchmarks/annotation/concatenate_cog.txt"
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/concatenate_cog.py"


rule lineage_parser:
    input:
        tax_out="output/annotation/cog/{sample}/{sample}_tax.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_tax.tsv",
        rankedlineage=config["lineage_parser"]["db"],
    output:
        # Class must be spelled with a 'k' to prevent conflicts with the Python keyword
        species="output/annotation/cog/{sample}/{sample}_species.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_species.tsv",
        genus="output/annotation/cog/{sample}/{sample}_genus.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_genus.tsv",
        family="output/annotation/cog/{sample}/{sample}_family.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_family.tsv",
        order="output/annotation/cog/{sample}/{sample}_order.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_order.tsv",
        klass="output/annotation/cog/{sample}/{sample}_class.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_class.tsv",
        phylum="output/annotation/cog/{sample}/{sample}_phylum.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_phylum.tsv",
        kingdom="output/annotation/cog/{sample}/{sample}_kingdom.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_kingdom.tsv",
        domain="output/annotation/cog/{sample}/{sample}_domain.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_domain.tsv",
    log:
        get_coassembly_benchmark_or_log("log", "annotation", "lineage_parser"),
    benchmark:
        get_coassembly_benchmark_or_log("benchmark", "annotation", "lineage_parser")
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/lineage_parser.py"


rule plot_cog:
    input:
        categories_file="output/annotation/cog/COG_categories_relative.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_categories.tsv",
        species="output/annotation/cog/COG_species_relative.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_species.tsv",
        genus="output/annotation/cog/COG_genus_relative.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_genus.tsv",
        family="output/annotation/cog/COG_family_relative.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_family.tsv",
        order="output/annotation/cog/COG_order_relative.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_order.tsv",
        klass="output/annotation/cog/COG_class_relative.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_class.tsv",
        phylum="output/annotation/cog/COG_phylum_relative.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_phylum.tsv",
        kingdom="output/annotation/cog/COG_kingdom_relative.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_kingdom.tsv",
        domain="output/annotation/cog/COG_domain_relative.tsv"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_domain.tsv",
    output:
        categories_plt="output/annotation/cog/COG_categories_relative.png"
        if not config["coassembly"]
        else "output/annotation/cog/coassembly_categories.png",
        taxa_barplots=get_taxa_plot_outputs(),
    params:
        filter_categories=config["plot_cog"]["filter_categories"],
        categories_cutoff=config["plot_cog"]["categories_cutoff"],
        tax_cutoff=config["plot_cog"]["tax_cutoff"],
        coassembly=config["coassembly"],
    log:
        "output/logs/annotation/plot_cog.log",
    benchmark:
        "output/benchmarks/annotation/plot_cog.txt"
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/plot_cog.py"
