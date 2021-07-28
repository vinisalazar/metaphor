"""
Mapping rules:

    - concatenate_contigs: concatenate contigs into a compressed file with vamb script
    - create_mapping: create mapping file from contig catalogue with minimap2
    - map_reads: map short reads against contigs with minimap2 and samtools
    - sort_reads: sort BAM file with samtools
"""


rule concatenate_contigs:
    input:
        contigs=expand(
            "output/assembly/megahit/{sample}/{sample}.contigs.fa",
            sample=[
                "readsa",
            ],
        ),
    output:
        catalogue="{output}/mapping/catalogue.fna.gz",
    params:
        sequence_length_cutoff=config["concatenate_contigs"]["sequence_length_cutoff"],
    log:
        "{output}/logs/mapping/concatenate_contigs.log",
    benchmark:
        "{output}/benchmarks/mapping/concatenate_contigs.txt"
    conda:
        "../envs/vamb.yaml"
    shell:
        """
        concatenate.py -m {params.sequence_length_cutoff} {output} {input} &> {log}
        """


rule create_mapping:
    input:
        catalogue_fna="{output}/mapping/catalogue.fna.gz",
    output:
        catalogue_idx="{output}/mapping/catalogue.mmi",
    log:
        "{output}/logs/mapping/create_mapping.log",
    benchmark:
        "{output}/benchmarks/mapping/create_mapping.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        minimap2 -d {output} {input} &> {log}
        """


rule map_reads:
    input:
        catalogue_idx="{output}/mapping/catalogue.mmi",
        reads="{output}/preprocess/interleave/{sample}-clean.fq",
    output:
        bam="{output}/mapping/bam/{sample}.bam",
    params:
        threads=workflow.cores,
        N=50,
        preset="sr",
        flags=3584,
    log:
        "{output}/logs/mapping/{sample}_map_reads.log",
    benchmark:
        "{output}/benchmarks/mapping/{sample}_map_reads.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        {{ minimap2 -t {params.threads} -N {params.N} -a -x {params.preset} \
                 {input.catalogue_idx} {input.reads} | samtools view \
                 -F {params.flags} -b --threads {params.threads} > {output.bam} ; }} &> {log}
        """


rule sort_reads:
    input:
        bam="{output}/mapping/bam/{sample}.bam",
    output:
        sort="{output}/mapping/bam/{sample}.sorted.bam",
    params:
        threads=workflow.cores,
    log:
        "{output}/logs/mapping/{sample}_sort_reads.log",
    benchmark:
        "{output}/benchmarks/mapping/{sample}_sort_reads.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        {{ samtools sort -@ {params.threads} -o {output.sort} {input.bam} ; }} &> {log}
        """


rule jgi_summarize_bam_contig_depths:
    input:
        expand("/output/mapping/bam/{sample}.sorted.bam", sample=sample_IDs)
    output: 
        output("{output}/mapping/bam_contig_depths.txt")
    log:
        "/{output}/logs/mapping/jgi_summarize_bam_contig_depths.log"
    benchmark:
        "/{output}/benchmarks/mapping/jgi_summarize_bam_contig_depths.txt"
    conda:
        "../envs/metabat2.yaml"
    shell: 