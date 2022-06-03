"""
mapping.smk

    Map reads to contigs. Files generated in this module are required for binning.

Mapping rules:
    - concatenate_contigs: concatenate contigs into a compressed file with vamb script
    - concatenate_contigs: concatenate proteins into a single file to be used by CONCOCT
    - create_index: create index file from contig catalogue with minimap2
    - map_reads: map short reads against contigs with minimap2 and samtools
    - sort_reads: sort BAM file with samtools
    - index_reads: index sorte BAM file with samtools
    - flagstat: calculate flagstat with samtools
    - jgi_summarize_bam_contig_depths: calculate coverage depth of each contig
"""
from pathlib import Path


rule concatenate_contigs:
    input:
        contigs=get_contigs_input(expand_=True),
    output:
        catalogue="output/mapping/{binning_group}/{binning_group}_contig_catalogue.fna.gz",
    params:
        sequence_length_cutoff=config["megahit"]["min_contig_len"],
    resources:
        mem_mb=get_max_mb(),
    wildcard_constraints:
        binning_group="cobinning"
        if config["cobinning"]
        else "|".join(binning_group_names),
    log:
        "output/logs/mapping/concatenate_contigs/{binning_group}.log",
    benchmark:
        "output/benchmarks/mapping/concatenate_contigs/{binning_group}.txt"
    conda:
        "../envs/vamb.yaml"
    shell:
        """
        concatenate.py -m {params.sequence_length_cutoff} {output} {input} &> {log}
        """


rule decompress_catalogue:
    input:
        catalogue_gz="output/mapping/{binning_group}/{binning_group}_contig_catalogue.fna.gz",
    output:
        catalogue="output/mapping/{binning_group}/{binning_group}_contig_catalogue.fna",
    threads: get_threads_per_task_size("small")
    resources:
        mem_mb=get_max_mb(),
    wildcard_constraints:
        binning_group="cobinning"
        if config["cobinning"]
        else "|".join(binning_group_names),
    log:
        "output/logs/mapping/decompress_catalogue/{binning_group}.log",
    benchmark:
        "output/benchmarks/mapping/decompress_catalogue/{binning_group}.txt"
    conda:
        "../envs/utils.yaml"
    shell:
        """
        pigz -d -f -p {threads} -k {input.catalogue_gz} &> {log}
        """


rule concatenate_proteins:
    """
    Used by DAS_Tool (skips the Prodigal run). This is deactivated for now.
    """
    input:
        proteins=expand(
            "output/annotation/prodigal/{group}/{group}_proteins.faa",
            group=group_names,
        ),
    output:
        prot_catalogue="output/mapping/{binning_group}_proteins_catalogue.faa",
    resources:
        mem_mb=get_max_mb(),
    wildcard_constraints:
        binning_group="cobinning"
        if config["cobinning"]
        else "|".join(binning_group_names),
    log:
        "output/logs/mapping/{binning_group}_concatenate_proteins.log",
    benchmark:
        "output/benchmarks/mapping/{binning_group}_concatenate_proteins.txt"
    conda:
        "../envs/utils.yaml"
    shell:
        """cat {input} > {output}"""


rule create_index:
    input:
        catalogue_fna="output/mapping/{binning_group}/{binning_group}_contig_catalogue.fna.gz",
    output:
        catalogue_idx="output/mapping/{binning_group}/{binning_group}_contig_catalogue.mmi",
    resources:
        mem_mb=get_max_mb(),
    wildcard_constraints:
        binning_group="cobinning"
        if config["cobinning"]
        else "|".join(binning_group_names),
    log:
        "output/logs/mapping/create_index/{binning_group}.log",
    benchmark:
        "output/benchmarks/mapping/create_index/{binning_group}.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        minimap2 -d {output} {input} &> {log}
        """


rule map_reads:
    input:
        catalogue_idx="output/mapping/{binning_group}/{binning_group}_contig_catalogue.mmi",
        fastq1=get_map_reads_input_R1,
        fastq2=get_map_reads_input_R2,
    output:
        # This one needs wildcard named 'sample' because of the input functions.
        bam="output/mapping/bam/{binning_group}/{sample}.map.bam",
    params:
        N=50,
        preset="sr",
        flags=3584,
        split_prefix="output/mapping/bam/{binning_group}/{sample}",
    threads: get_threads_per_task_size("big")
    resources:
        mem_mb=get_mb_per_cores,
    wildcard_constraints:
        binning_group="cobinning"
        if config["cobinning"]
        else "|".join(binning_group_names),
        sample="|".join(sample_IDs),
    log:
        "output/logs/mapping/map_reads/{binning_group}-{sample}.log",
    benchmark:
        "output/benchmarks/mapping/map_reads/{binning_group}-{sample}.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        {{ minimap2 -t {threads}                        \
                    -N {params.N}                       \
                    -ax {params.preset}                 \
                    --split-prefix {params.split_prefix}\
                    {input.catalogue_idx}               \
                    {input.fastq1}                      \
                    {input.fastq2} ; }} 2>> {log}       |
        {{ samtools view                                \
                    -F {params.flags}                   \
                    -b --threads                        \
                    {threads} > {output.bam} ; }} 2>> {log}
        """


rule sort_reads:
    input:
        bam="output/mapping/bam/{binning_group}/{sample}.map.bam",
    output:
        sort="output/mapping/bam/{binning_group}/{sample}.sorted.bam",
    threads: get_threads_per_task_size("big")
    resources:
        mem_mb=get_mb_per_cores,
    wildcard_constraints:
        binning_group="cobinning"
        if config["cobinning"]
        else "|".join(binning_group_names),
        sample="|".join(sample_IDs),
    log:
        "output/logs/mapping/sort_reads/{binning_group}-{sample}.log",
    benchmark:
        "output/benchmarks/mapping/sort_reads/{binning_group}-{sample}.txt"
    conda:
        "../envs/samtools.yaml"
    wrapper:
        get_wrapper("samtools/sort")


rule index_reads:
    input:
        sort="output/mapping/bam/{binning_group}/{sample}.sorted.bam",
    output:
        index="output/mapping/bam/{binning_group}/{sample}.sorted.bam.bai",
    threads: get_threads_per_task_size("big")
    resources:
        mem_mb=get_mb_per_cores,
    wildcard_constraints:
        binning_group="cobinning"
        if config["cobinning"]
        else "|".join(binning_group_names),
        sample="|".join(sample_IDs),
    log:
        "output/logs/mapping/index_reads/{binning_group}-{sample}.log",
    benchmark:
        "output/benchmarks/mapping/index_reads/{binning_group}-{sample}.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        "samtools index -@ {threads} {input} {output} &> {log}"


rule flagstat:
    input:
        sort="output/mapping/bam/{binning_group}/{sample}.sorted.bam",
    output:
        flagstat="output/mapping/bam/{binning_group}/{sample}.flagstat.txt",
    threads: get_threads_per_task_size("big")
    resources:
        mem_mb=get_mb_per_cores,
    wildcard_constraints:
        binning_group="cobinning"
        if config["cobinning"]
        else "|".join(binning_group_names),
        sample="|".join(sample_IDs),
    log:
        "output/logs/mapping/flagstat/{binning_group}/{sample}.log",
    benchmark:
        "output/benchmarks/mapping/flagstat/{binning_group}/{sample}.txt"
    conda:
        "../envs/samtools.yaml"
    shell:
        "{{ samtools flagstat -@ {threads} {input.sort} > {output.flagstat} ; }} &> {log}"


rule jgi_summarize_bam_contig_depths:
    input:
        lambda wildcards: expand(
            "output/mapping/bam/{binning_group}/{sample}.sorted.bam",
        binning_group=wildcards.binning_group,
        sample=samples.query(f"binning_group == '{wildcards.binning_group}'")[
        "sample_name"
            ].unique(),
        ),
    output:
        contig_depths="output/mapping/{binning_group}/bam_contig_depths.txt",
    wildcard_constraints:
        binning_group="cobinning"
        if config["cobinning"]
        else "|".join(binning_group_names),
    log:
        "output/logs/mapping/jgi_summarize_bam_contig_depths/{binning_group}.log",
    benchmark:
        "output/benchmarks/mapping/jgi_summarize_bam_contig_depths/{binning_group}.txt"
    conda:
        "../envs/metabat2.yaml"
    shell:
        "jgi_summarize_bam_contig_depths {input} --outputDepth {output} 2> {log}"
