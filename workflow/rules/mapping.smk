"""
Mapping rules:

    - concatenate_contigs: concatenate contigs into a compressed file with vamb script
    - create_mapping: create mapping file from contig catalogue with minimap2
    - map_reads: map short reads against contigs with minimap2 and samtools
    - sort_reads: sort BAM file with samtools
    - flagstat: calculate flagstat with samtools
"""


rule concatenate_contigs:
    input:
        contigs=expand(
            "output/assembly/megahit/{sample}/{sample}.contigs.fa", sample=sample_IDs,
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
        fastq1=get_map_reads_input_R1,
        fastq2=get_map_reads_input_R2
    output:
        bam="{output}/mapping/bam/{sample}.map.bam",
    params:
        N=50,
        preset="sr",
        flags=3584,
    threads: workflow.cores
    log:
        "{output}/logs/mapping/{sample}_map_reads.log",
    benchmark:
        "{output}/benchmarks/mapping/{sample}_map_reads.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        {{ minimap2 -t {threads} -N {params.N} -a -x {params.preset} \
                 {input.catalogue_idx} {input.fastq1} {input.fastq2} \
                 | samtools view -F {params.flags} -b --threads \
                   {threads} > {output.bam} ; }} &> {log}
        """


rule sort_reads:
    input:
        bam="{output}/mapping/bam/{sample}.map.bam",
    output:
        sort="{output}/mapping/bam/{sample}.sorted.bam",
    threads: int(workflow.cores * 0.25)
    log:
        "{output}/logs/mapping/{sample}_sort_reads.log",
    benchmark:
        "{output}/benchmarks/mapping/{sample}_sort_reads.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        {{ samtools sort -@ {threads} -o {output.sort} {input.bam} ; }} &> {log}
        """


rule flagstat:
    input: 
        sort="{output}/mapping/bam/{sample}.sorted.bam",
    output:
        flagstat="{output}/mapping/bam/{sample}.flagstat.txt",
    threads: int(workflow.cores * 0.25) 
    log:
        "{output}/logs/mapping/{sample}_flagstat.log",
    benchmark:
        "{output}/benchmarks/mapping/{sample}_flagstat.txt"
    conda:
        "../envs/samtools.yaml"
    shell: 
        "{{ samtools flagstat -@ {threads} {input.sort} > {output.flagstat} ; }} &> {log}"


# WIP - for vamb step
rule jgi_summarize_bam_contig_depths:
    input:
        expand("output/mapping/bam/{sample}.sorted.bam", sample=sample_IDs),
    output:
        contig_depths="{output}/mapping/bam_contig_depths.txt",
    log:
        "{output}/logs/mapping/jgi_summarize_bam_contig_depths.log",
    benchmark:
        "{output}/benchmarks/mapping/jgi_summarize_bam_contig_depths.txt"
    conda:
        "../envs/metabat2.yaml"
    shell:
        "jgi_summarize_bam_contig_depths {input} --outputDepth {output}"
