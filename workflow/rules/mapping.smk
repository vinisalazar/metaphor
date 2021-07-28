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
            "output/megahit/{sample}/{sample}.contigs.fa",
            sample=[
                "readsa",
            ],
        ),
    output:
        catalogue="{output}/vamb/catalogue.fna.gz",
    params:
        sequence_length_cutoff=config["concatenate_contigs"]["sequence_length_cutoff"],
    log:
        "{output}/logs/vamb/concatenate_contigs.log",
    benchmark:
        "{output}/benchmarks/vamb/concatenate_contigs.txt"
    conda:
        "../envs/vamb.yaml"
    shell:
        """
        concatenate.py -m {params.sequence_length_cutoff} {output} {input} &> {log}
        """


rule create_mapping:
    input:
        catalogue_fna="{output}/vamb/catalogue.fna.gz",
    output:
        catalogue_idx="{output}/vamb/catalogue.mmi",
    log:
        "{output}/logs/vamb/create_mapping.log",
    benchmark:
        "{output}/benchmarks/vamb/create_mapping.txt"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        minimap2 -d {output} {input} &> {log}
        """


rule map_reads:
    input:
        catalogue_idx="{output}/vamb/catalogue.mmi",
        reads="{output}/interleave/{sample}-clean.fq",
    output:
        bam="{output}/vamb/bam/{sample}.bam",
    params:
        threads=workflow.cores,
        N=50,
        preset="sr",
        flags=3584,
    log:
        "{output}/logs/vamb/{sample}_map_reads.log",
    benchmark:
        "{output}/benchmarks/vamb/{sample}_map_reads.txt"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        {{ minimap2 -t {params.threads} -N {params.N} -a -x {params.preset} \
                 {input.catalogue_idx} {input.reads} | samtools view \
                 -F {params.flags} -b --threads {params.threads} > {output.bam} ; }} &> {log}
        """


rule sort_reads:
    input:
        bam="{output}/vamb/bam/{sample}.bam",
    output:
        sort="{output}/vamb/bam/{sample}.sorted.bam",
    params:
        threads=workflow.cores,
    log:
        "{output}/logs/vamb/{sample}_sort_reads.log",
    benchmark:
        "{output}/benchmarks/vamb/{sample}_sort_reads.txt"
    conda:
        "../envs/bwa.yaml"
    shell:
        """
        {{ samtools sort -@ {params.threads} -o {output.sort} {input.bam} ; }} &> {log}
        """
