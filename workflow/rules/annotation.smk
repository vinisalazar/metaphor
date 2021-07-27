"""
Annotation rules:

    - prodigal: gene prediction with Prodigal
    - diamond: protein annotation with Diamond
    - collation: format Diamond XML outputs with sed
    - xmlparser: parse collated outputs with custom scripts
"""

from pathlib import Path


rule prodigal:
    input:
        contigs="{output}/megahit/{sample}/{sample}.contigs.fa"
    output:
        genes="{output}/prodigal/{sample}/{sample}_genes.fna",
        proteins="{output}/prodigal/{sample}/{sample}_proteins.faa",
        scores="{output}/prodigal/{sample}/{sample}_scores.cds",
        genbank="{output}/prodigal/{sample}/{sample}_genbank.gbk"
    log:
        "{output}/logs/prodigal/{sample}"
    benchmark:
        "{output}/benchmarks/prodigal/{sample}.txt"
    conda:
        "../envs/prodigal.yaml"
    shell:
        """
        prodigal -i {input} \
                 -d {output.genes} \
                 -a {output.proteins} \
                 -s {output.scores} \
                 -o {output.genbank} &> {log}
        """


rule diamond:
    input:
        proteins="{output}/prodigal/{sample}/{sample}_proteins.faa"
    output:
        xmlout="{output}/diamond/{sample}.xml"
    params:
        db=config["diamond_db"],
        max_target_seqs=1,
        format=5,
        cpus=workflow.cores
    log:
        "{output}/logs/diamond/{sample}.log"
    benchmark:
        "{output}/benchmarks/diamond/{sample}.txt"
    shell:
        """
        diamond blastp -q {input} \
                --max-target-seqs {params.max_target_seqs} \
                -p {params.cpus} \
                -f {params.format} \
                -d {params.db} \
                -o {output} &> {log}
        """


rule collation:
    input: 
        xmlout="{output}/diamond/{sample}.xml"
    output:
        collationout="{output}/collation/{sample}-collated.xml"
    log:
        "{output}/logs/collation/{sample}-collation.log"
    benchmark:
        "{output}/benchmarks/collation/{sample}-collation.txt"
    conda:
        "../envs/bash.yaml"
    shell: 
        """
        sed 's/\&quot;//g' '{input}' | sed 's/\&//g' > {output}
        """


rule xmlparser:
    input: 
        collationout="{output}/collation/{sample}-collated.xml"
    output:
        gene_count_table="{output}/collation/{sample}_gene_count_table.txt",
        otu_table="{output}/collation/{sample}_OTU_out.txt"
    params:
        output_dir=str(Path("{input}").parent),
        ko_formatted_file="bin/db/formatted.xml.out",
        kegg_species_file="bin/db/species_prokaryotes.dat",
        tax_rank_file="bin/db/tax_rank",
        full_lineage_file="fullnamelineage.dmp"
    log:
        "{output}/logs/collation/{sample}-xmlparser.log"
    benchmark:
        "{output}/benchmarks/collation/{sample}-xmlparser.txt"
    conda:
        "../envs/bash.yaml"
    shell:
        """
        perl scripts/xml_parser.function.pl \
             {output.gene_count_table} 1 \
             {params.ko_formatted_file} \
             {params.kegg_species_file} \
             {input}

        perl scripts/orgID_2_name.pl \
             {params.tax_rank_file} \
             {params.full_lineage_file} \
             {params.output_dir} > {output.otu_table}
        """  # must add 'sep' next to input (see L30) of metaGenePipe/tasks/xml_parser.wdl
