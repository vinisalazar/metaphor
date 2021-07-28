"""
Annotation rules:

    - prodigal: gene prediction with Prodigal
    - diamond: protein annotation with Diamond
    - xmlparser: parse collated outputs with custom scripts
"""

from pathlib import Path


rule prodigal:
    input:
        contigs="{output}/megahit/{sample}/{sample}.contigs.fa",
    output:
        genes="{output}/prodigal/{sample}/{sample}_genes.fna",
        proteins="{output}/prodigal/{sample}/{sample}_proteins.faa",
        scores="{output}/prodigal/{sample}/{sample}_scores.cds",
        genbank="{output}/prodigal/{sample}/{sample}_genbank.gbk",
    params:
        mode=config["prodigal"]["mode"],
    log:
        "{output}/logs/prodigal/{sample}",
    benchmark:
        "{output}/benchmarks/prodigal/{sample}.txt"
    conda:
        "../envs/prodigal.yaml"
    shell:
        """
        prodigal -p {params.mode} \
                 -i {input} \
                 -d {output.genes} \
                 -a {output.proteins} \
                 -s {output.scores} \
                 -o {output.genbank} &> {log}
        """


rule diamond:
    input:
        proteins="{output}/prodigal/{sample}/{sample}_proteins.faa",
    output:
        xmlout="{output}/diamond/{sample}.xml",
    params:
        db=config["diamond"]["db"],
        max_target_seqs=1,
        format=5,
        cpus=workflow.cores,
    log:
        "{output}/logs/diamond/{sample}.log",
    benchmark:
        "{output}/benchmarks/diamond/{sample}.txt"
    conda:
        "../envs/diamond.yaml"
    shell:
        """
        {{ diamond blastp -q {input} \
                   --max-target-seqs {params.max_target_seqs} \
                   -p {params.cpus} \
                   -f {params.format} \
                   -d {params.db} \
                   | sed 's/\&quot;//g' \
                   | sed 's/\&//g' > {output} ; }} &> {log}
        """


rule xmlparser:
    input:
        xmlout="{output}/diamond/{sample}.xml",
    output:
        gene_count_table="{output}/diamond/{sample}_gene_count_table.txt",
        otu_table="{output}/diamond/{sample}_OTU_out.txt",
    params:
        output_dir=str(Path("{input}").parent),
        ko_formatted_file="bin/db/formatted.xml.out",
        kegg_species_file="bin/db/species_prokaryotes.dat",
        tax_rank_file="bin/db/tax_rank",
        full_lineage_file="fullnamelineage.dmp",
        xml_parser=pathfinder("../scripts/xml_parser.function.pl"),
        orgID_2_name=pathfinder("../scripts/orgID_2_name.pl")
    log:
        "{output}/logs/diamond/{sample}-xmlparser.log",
    benchmark:
        "{output}/benchmarks/diamond/{sample}-xmlparser.txt"
    conda:
        "../envs/bash.yaml"
    shell:
        """
        perl {params.xml_parser} \
             {output.gene_count_table} 1 \
             {params.ko_formatted_file} \
             {params.kegg_species_file} \
             {input}

        perl {params.orgID_2_name} \
             {params.tax_rank_file} \
             {params.full_lineage_file} \
             {params.output_dir} > {output.otu_table}
        """  # must add 'sep' next to input (see L30) of metaGenePipe/tasks/xml_parser.wdl
