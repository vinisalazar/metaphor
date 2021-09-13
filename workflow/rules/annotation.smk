"""
Annotation rules:

    - prodigal: gene prediction with Prodigal
    - hmmsearch: find KEGG categories with hmmsearch
    - diamond: protein annotation with Diamond
    - diamond_parser: parse collated outputs with custom Python script
    - hmmer_parser: parse the output of hmmsearch with custom Python script
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
    log:
        "{output}/logs/annotation/prodigal/{sample}.log",
    benchmark:
        "{output}/benchmarks/annotation/prodigal/{sample}.txt"
    conda:
        "../envs/prodigal.yaml"
    shell:
        """
        prodigal -p {params.mode}       \
                 -i {input}             \
                 -d {output.genes}      \
                 -a {output.proteins}   \
                 -s {output.scores}     \
                 -o {output.genbank} &> {log}
        """


rule hmmsearch:
    input:
        fasta="{output}/annotation/prodigal/{sample}/{sample}_proteins.faa",
        profile=config["hmmsearch"]["db"],
    output:
        # only one of these is required
        tblout="{output}/annotation/hmmsearch/{sample}_hmmer.tblout",  # save parseable table of per-sequence hits to file <f>
        outfile="{output}/annotation/hmmsearch/{sample}_hmmer.out",  # Direct the main human-readable output to a file <f> instead of the default stdout.
        # domtblout="{output}/annotation/hmmsearch/{sample}_hmmer.domtblout", # save parseable table of per-domain hits to file <f>
        # alignment is disabled because it's too big
        # alignment_hits="{output}/annotation/hmmsearch/{sample}_hmmer.aln", # Save a multiple alignment of all significant hits (those satisfying inclusion thresholds) to the file <f>
    log:
        "{output}/logs/annotation/hmmsearch/{sample}.log",
    benchmark:
        "{output}/benchmarks/annotation/hmmsearch/{sample}.txt"
    params:
        evalue_threshold=0.00001,
        # if bitscore threshold provided, hmmsearch will use that instead
        # score_threshold=50,
        extra="",
    threads: round(workflow.cores * 0.75)
    wrapper:
        "0.77.0/bio/hmmer/hmmsearch"


rule diamond:
    input:
        proteins="{output}/annotation/prodigal/{sample}/{sample}_proteins.faa",
    output:
        dmnd_out="{output}/annotation/diamond/{sample}_dmnd.out",
    params:
        db=config["diamond"]["db"],
        max_target_seqs=1,
        output_type=config["diamond"]["output_type"],
        output_format=config["diamond"]["output_format"],
    threads: round(workflow.cores * 0.75)
    log:
        "{output}/logs/annotation/diamond/{sample}.log",
    benchmark:
        "{output}/benchmarks/annotation/diamond/{sample}.txt"
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


rule diamond_parser:
    input:
        dmnd_out=get_diamond_output(),
    output:
        outfile=get_diamond_parser_output(),
    params:
        db=config["diamond_parser"]["db"],
    log:
        "output/logs/annotation/diamond_parser/diamond_parser.log",
    benchmark:
        "output/benchmarks/annotation/diamond_parser/diamond_parser.txt"
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/diamond_parser.py"


rule hmmer_parser:
    input:
        hmm_tbls=get_hmmsearch_output(),
    output:
        brite_level1="{output}/annotation/brite/{sample}_brite_Level1.tsv",
        brite_level2="{output}/annotation/brite/{sample}_brite_Level2.tsv",
        brite_level3="{output}/annotation/brite/{sample}_brite_Level3.tsv",
    params:
        brite=config["hmmer_parser"]["db"],
        consistent_pathways=config["hmmer_parser"]["consistent_pathways"],
        outprefix=lambda w, output: output.brite_level1.split("_brite_Level")[0],  # Cannot provide output file prefix, must infer instead.
    log:
        "{output}/logs/annotation/hmmer_parser/{sample}.log",
    benchmark:
        "{output}/benchmarks/annotation/hmmer_parser/{sample}.txt"
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/hmmer_parser.py"
