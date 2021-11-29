"""
Annotation rules:

    - prodigal: gene prediction with Prodigal
    - diamond: protein annotation with Diamond
    - cog_parser: parse collated outputs with custom Python script
"""

from pathlib import Path


rule prodigal:
    input:
        contigs="{output}/assembly/megahit/{sample}/{sample}.contigs.fa",
    output:
        genes="{output}/annotation/prodigal/{sample}/{sample}_genes.fna",
        proteins="{output}/annotation/prodigal/{sample}/{sample}_proteins.faa",
        scores="{output}/annotation/prodigal/{sample}/{sample}_scores.cds",
        genbank="{output}/annotation/prodigal/{sample}/{sample}_genbank.gbk",
    params:
        mode=config["prodigal"]["mode"],
        quiet="-q" if config["prodigal"]["quiet"] else "",
    log:
        "{output}/logs/annotation/prodigal/{sample}.log",
    benchmark:
        "{output}/benchmarks/annotation/prodigal/{sample}.txt"
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


rule prokka:
    input:
        contigs="output/assembly/megahit/{sample}/{sample}.contigs.fa",
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


rule diamond:
    input:
        proteins="output/annotation/prodigal/{sample}/{sample}_proteins.faa",
    output:
        dmnd_out=get_diamond_output(),
    params:
        db=config["diamond"]["db"],
        max_target_seqs=1,
        output_type=config["diamond"]["output_type"],
        output_format=config["diamond"]["output_format"],
    threads: round(workflow.cores * 0.75)
    log:
        "output/logs/annotation/diamond/{sample}.log",
    benchmark:
        "output/benchmarks/annotation/diamond/{sample}.txt"
    conda:
        "../envs/diamond.yaml"
    shell:
        """
        echo {params.output_format} | sed -e 's/ /\t/g' > {output}
        {{ diamond blastp -q {input}                            \
                   --max-target-seqs {params.max_target_seqs}   \
                   -p {threads}                                 \
                   -d {params.db}                               \
                   -f {params.output_type}                      \
                   {params.output_format}                       \
                   >> {output} ; }} &> {log}
        """


rule cog_parser:
    input:
        dmnd_out=get_diamond_output(),
    output:
        categories_out="output/annotation/cog/{sample}/{sample}_categories.tsv",
        codes_out="output/annotation/cog/{sample}/{sample}_codes.tsv",
        tax_out="output/annotation/cog/{sample}/{sample}_tax.tsv",
        pathways_out="output/annotation/cog/{sample}/{sample}_pathways.tsv",
    params:
        cog_csv=get_cog_db_file("cog-20.cog.csv*"),
        def_tab=get_cog_db_file("cog-20.def.tab*"),
        fun_tab=get_cog_db_file("fun-20.tab*"),
    log:
        "output/logs/annotation/cog_parser/{sample}.log",
    benchmark:
        "output/benchmarks/annotation/cog_parser/{sample}.txt"
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
        tax_out="output/annotation/cog/{sample}/{sample}_tax.tsv",
        rankedlineage=config["lineage_parser"]["db"],
    output:
        # Class must be spelled with a 'k' to prevent conflicts with the Python keyword
        species="output/annotation/cog/{sample}/{sample}_species.tsv",
        genus="output/annotation/cog/{sample}/{sample}_genus.tsv",
        family="output/annotation/cog/{sample}/{sample}_family.tsv",
        order="output/annotation/cog/{sample}/{sample}_order.tsv",
        klass="output/annotation/cog/{sample}/{sample}_class.tsv",
        phylum="output/annotation/cog/{sample}/{sample}_phylum.tsv",
        kingdom="output/annotation/cog/{sample}/{sample}_kingdom.tsv",
        domain="output/annotation/cog/{sample}/{sample}_domain.tsv",
    log:
        "output/logs/annotation/lineage_parser/{sample}.log",
    benchmark:
        "output/benchmarks/annotation/lineage_parser/{sample}.txt"
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/lineage_parser.py"
