# MetaSnakePipe
# Original MetaGenePipe workflow by Bobbie Shaban
# Snakemake port by Vini Salazar
from pathlib import Path


configfile: "config.yaml"


rule all:
    input: 
        expand("data/{sample}-{readsno}.fq",
                sample=["readsa",],
                readsno=["1", "2"]) 


rule fastqc_raw:  # qc on raw reads
    input:
        "data/{sample}.fq"
    output:
        html="{output}/qc/{sample}.html",
        zip="{output}/qc/{sample}_fastqc.zip"
    params: "--quiet"
    log:
        "{output}/logs/qc/{sample}-fastqc.log"
    benchmark:
        "{output}/benchmarks/qc/{sample}_fastqc.txt"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"


rule flash:
    input:
        fqforward="data/{sample}-1.fq",
        fqreverse="data/{sample}-2.fq" 
    output: 
        flash_notcombined1="{output}/flash/{sample}.notCombined_1.fastq",
        flash_notcombined2="{output}/flash/{sample}.notCombined_2.fastq",
        flash_extended="{output}/flash/{sample}.extendedFrags.fastq"
    params:
        max_overlap=120
    log:
        "{output}/logs/qc/{sample}-flash.log"
    benchmark:
        "{output}/benchmarks/flash/{sample}.txt"
    conda:
        "envs/flash.yaml"
    shell: 
        """
        flash -d {wildcards.output}/flash -o {wildcards.sample} \
              -M {params.max_overlap} {input} &> {log}
        """


rule interleave:
    input: 
        flash_notcombined1="{output}/flash/{sample}.notCombined_1.fastq",
        flash_notcombined2="{output}/flash/{sample}.notCombined_2.fastq",
        flash_extended="{output}/flash/{sample}.extendedFrags.fastq"
    output: 
        interleaved="{output}/interleave/{sample}-interleaved.fq",
        merged="{output}/interleave/{sample}-merged.fq"
    log:
        "{output}/logs/interleave/{sample}-interleave.log"
    benchmark:
        "{output}/benchmarks/interleave/{sample}-interleave.txt"
    conda:
        "envs/bash.yaml"
    shell: 
        """
        bash scripts/interleave_fastq.sh {input.flash_notcombined1} {input.flash_notcombined2} > {output.interleaved}

        cat {output.interleaved} {input.flash_extended} > {output.merged}
        """

rule fastqc_merged:  # qc on merged reads, after rules 'flash' and 'interleave'
    input:
        "{output}/interleave/{sample}-merged.fq"
    output:
        html="{output}/qc/{sample}-merged.html",
        zip="{output}/qc/{sample}-merged_fastqc.zip"
    params: "--quiet"
    log:
        "{output}/logs/qc/{sample}-fastqc_merged.log"
    benchmark:
        "{output}/benchmarks/qc/{sample}_fastqc_merged.txt"
    threads: 1
    wrapper:
        "0.77.0/bio/fastqc"


rule multiqc:
    input:
        expand("output/qc/{sample}-{kind}_fastqc.zip", sample=["readsa", ], kind=["1", "2", "merged"])
    output:
        report="{output}/qc/multiqc.html"
    log:
        "{output}/logs/qc/multiqc.log"
    benchmark:
        "{output}/benchmarks/qc/multiqc.txt"
    wrapper:
        "0.77.0/bio/multiqc"


rule hostremoval:
    input:
        flash_merged="{output}/interleave/{sample}-merged.fq"
    params:
        alignment_id_threshold=70,
        alignment_coverage_threshold=70,
        dbs="mm1,mm2,mm3,mm4,mm5,mm6"
    output:
        host_removal_output="{output}/interleave/{sample}-clean.fq"
    log:
        "{output}/logs/interleave/{sample}-hostremoval.log"
    benchmark:
        "{output}/benchmarks/interleave/{sample}-hostremoval.txt"
    conda:
        "envs/bash.yaml"
    shell:
        "cat {input} > {output}"
        # This rule is off for now due to the lack of databases
        # "perl bin/dqc/deconseq.pl -dbs {params.dbs}" \
        # " -f {input}" \
        # " -i {params.alignment_id_threshold}" \
        # " -c {params.alignment_coverage_threshold} " \
        # " -id {wildcards.sample}"


rule megahit:
    input:
        host_removal_output="{output}/interleave/{sample}-clean.fq"
    output:
        contigs="{output}/megahit/{sample}/{sample}.contigs.fa"
    params:
        out_dir=lambda w, output: str(Path(output.contigs).parent.parent),  # this is equivalent to "{output}/megahit"
        min_contig_len=200,
        k_list="21,29",
        memory=0.5,
        cpus=workflow.cores
    log:
        "{output}/logs/megahit/{sample}-megahit.log"
    benchmark:
        "{output}/benchmarks/megahit/{sample}.txt"
    conda:
        "envs/megahit.yaml"
    
    # Using the '--12' flag yielded slightly better results than the '-r' flag
    shell:
        """
        # MegaHit has no --force flag, so we must remove the created directory prior to running
        rm -rf {params.out_dir}/{wildcards.sample}

        megahit --12 {input} -o {params.out_dir}/{wildcards.sample} \
                --out-prefix {wildcards.sample} \
                --min-contig-len {params.min_contig_len}  \
                -t {params.cpus}  \
                -m {params.memory} \
                --k-list {params.k_list} &> {log}
        """


rule vamb_concatenate:
    input:
        contigs=expand("output/megahit/{sample}/{sample}.contigs.fa", sample=["readsa",]),
        reads=expand("output/interleave/{sample}-clean.fq")
    output:
        catalogue_fna="{output}/vamb/catalogue.fna.gz",
    params:
        threads=workflow.cores,
        minimap2_N=50,
        minimap2_preset="sr",
    log:
        "{output}/logs/vamb/vamb_concatenate.log"
    benchmark:
        "{output}/benchmarks/vamb/vamb_concatenate.log"
    conda:
        "envs/vamb.yaml" 
    shell: 
        """
        concatenate.py {output.catalogue} {input.contigs}
        """


rule vamb_create_map:
    input: 
        catalogue_fna="{output}/vamb/catalogue.fna.gz"
    output: 
        catalogue_idx="{output}/vamb/catalogue.mmi"
    log:
        "{output}/logs/vamb/vamb_create_map.log"
    benchmark:
        "{output}/benchmarks/vamb/vamb_create_map.log"
    conda:
        "envs/vamb.yaml"
    shell:
        """
        minimap2 -d {output} {input}
        """


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
        "envs/prodigal.yaml"
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
