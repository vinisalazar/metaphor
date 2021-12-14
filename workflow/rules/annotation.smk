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
        genes="output/annotation/prodigal/{sample}/{sample}_genes.fna"
        if not config["coassembly"]
        else "output/annotation/prodigal/coassembly_genes.fna",
        proteins="output/annotation/prodigal/{sample}/{sample}_proteins.faa"
        if not config["coassembly"]
        else "output/annotation/prodigal/coassembly_proteins.faa",
        scores="output/annotation/prodigal/{sample}/{sample}_scores.cds"
        if not config["coassembly"]
        else "output/annotation/prodigal/coassembly_scores.cds",
        genbank="output/annotation/prodigal/{sample}/{sample}_genbank.gbk"
        if not config["coassembly"]
        else "output/annotation/prodigal/coassembly_genbank.gbk",
    params:
        mode=config["prodigal"]["mode"],
        quiet="-q" if config["prodigal"]["quiet"] else "",
    log:
        "output/logs/annotation/prodigal/{sample}.log"
        if not config["coassembly"]
        else "output/logs/annotation/prodigal/coassembly.log",
    benchmark:
        "output/benchmarks/annotation/prodigal/{sample}.txt" if not config[
        "coassembly"
        ] else "output/benchmarks/annotation/prodigal/coassembly.txt"
    conda:
        "../envs/prodigal.yaml"
    shell:
        """
        prodigal {params.quiet}         \
                 -p {params.mode}       \
                 -i {input}             \
                 -d {output.genes}      \
                 -a {output.proteins}   \
                 -s {output.scores}     \
                 -o {output.genbank} &> {log}
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
        cog_fasta=get_cog_db_file("cog-20.fa.gz"),
        cog_csv=get_cog_db_file("cog-20.cog.csv"),
        def_tab=get_cog_db_file("cog-20.def.tab"),
        fun_tab=get_cog_db_file("fun-20.tab"),
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


rule diamond_makedb:
    input:
        fname=get_cog_db_file("cog-20.fa.gz"),
    output:
        fname=config["diamond"]["db"],
    log:
        "output/logs/annotation/diamond/diamond_makedb.log",
    threads: round(workflow.cores * 0.25)
    wrapper:
        str(Path(config["wrapper_version"]).joinpath("bio/diamond/makedb"))


rule download_taxonomy_database:
    output:
        # Replace the .test/ directory with data/ directory if it is set like so
        rankedlineage=config["lineage_parser"]["db"],
    params:
        download_url=config["lineage_parser"]["download_url"],
        output_dir=lambda w, input_: str(Path(input_.rankedlineage).parent),
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

        ls {output.rankedlineage} 2>> {log}
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
        max_target_seqs=25,
        output_type=config["diamond"]["output_type"],
        output_format=config["diamond"]["output_format"],
        extra_params="--iterate --top 0",
    threads: round(workflow.cores * 0.75)
    log:
        "output/logs/annotation/diamond/{sample}.log"
        if not config["coassembly"]
        else "output/logs/annotation/diamond/coassembly.log",
    benchmark:
        "output/benchmarks/annotation/diamond/{sample}.txt" if not config[
        "coassembly"
        ] else "output/benchmarks/annotation/diamond/coassembly.txt"
    conda:
        "../envs/diamond.yaml"
    shell:
        """
        echo {params.output_format} | sed -e 's/ /\t/g' > {output.fname}
        {{ diamond blastp -q {input.fname_fasta}                \
                   --max-target-seqs {params.max_target_seqs}   \
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
        "output/logs/annotation/cog_parser/{sample}.log"
        if not config["coassembly"]
        else "output/logs/annotation/cog_parser/coassembly.log",
    benchmark:
        "output/benchmarks/annotation/cog_parser/{sample}.txt" if not config[
        "coassembly"
        ] else "output/benchmarks/annotation/cog_parser/coassembly.txt"
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
    output:
        concat_categories_absolute="output/annotation/cog/COG_categories_absolute.tsv",
        concat_categories_relative="output/annotation/cog/COG_categories_relative.tsv",
        concat_codes_absolute="output/annotation/cog/COG_codes_absolute.tsv",
        concat_codes_relative="output/annotation/cog/COG_codes_relative.tsv",
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
        "output/logs/annotation/lineage_parser/{sample}.log"
        if not config["coassembly"]
        else "output/logs/annotation/lineage_parser/coassembly.log",
    benchmark:
        "output/benchmarks/annotation/lineage_parser/{sample}.txt" if not config[
        "coassembly"
        ] else "output/benchmarks/annotation/lineage_parser/coassembly.txt"
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/lineage_parser.py"
