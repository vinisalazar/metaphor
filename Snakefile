# MetaSnakePipe
# Original MetaGenePipe workflow by Bobbie Shaban
# Snakemake port by Vini Salazar (with modifications)
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
        {{ bash scripts/interleave_fastq.sh {input.flash_notcombined1} {input.flash_notcombined2} > {output.interleaved} ; }} &> {log}

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


rule concatenate_contigs:
    input:
        contigs=expand("output/megahit/{sample}/{sample}.contigs.fa", sample=["readsa",])
    output:
        catalogue="{output}/vamb/catalogue.fna.gz"
    params:
        sequence_length_cutoff=2000  # vamb's default
    log:
        "{output}/logs/vamb/concatenate_contigs.log"
    benchmark:
        "{output}/benchmarks/vamb/concatenate_contigs.txt"
    conda:
        "envs/vamb.yaml" 
    shell: 
        """
        concatenate.py -m {params.sequence_length_cutoff} {output} {input} &> {log}
        """


rule create_mapping:
    input: 
        catalogue_fna="{output}/vamb/catalogue.fna.gz"
    output: 
        catalogue_idx="{output}/vamb/catalogue.mmi"
    log:
        "{output}/logs/vamb/create_mapping.log"
    benchmark:
        "{output}/benchmarks/vamb/create_mapping.txt"
    conda:
        "envs/bwa.yaml"
    shell:
        """
        minimap2 -d {output} {input} &> {log}
        """


rule map_reads:
    input: 
        catalogue_idx="{output}/vamb/catalogue.mmi",
        reads="{output}/interleave/{sample}-clean.fq"
    output: 
        bam="{output}/vamb/bam/{sample}.bam"
    params:
        threads=workflow.cores,
        N=50,
        preset="sr",
        flags=3584
    log:
        "{output}/logs/vamb/{sample}_map_reads.log"
    benchmark:
        "{output}/benchmarks/vamb/{sample}_map_reads.txt"
    conda:
        "envs/bwa.yaml"
    shell:
        """
        {{ minimap2 -t {params.threads} -N {params.N} -a -x {params.preset} \
                 {input.catalogue_idx} {input.reads} | samtools view \
                 -F {params.flags} -b --threads {params.threads} > {output} ; }} &> {log}
        """


rule vamb:
    input: 
        bamfiles=expand("output/vamb/bam/{sample}", sample=["readsa",]),
        catalogue="{output}/vamb/catalogue.fna.gz"
    output: 
        outdir=directory("output/vamb/vamb")
    params:  # defaults in vamb's README
        binsplit_sep="C",
        minfasta=200000
    log:
        "output/log/vamb/vamb.log"
    benchmark:
        "output/benchmarks/vamb/vamb.txt"
    conda:
        "envs/vamb.yaml"
    shell: 
        """
        rm -rf {output}

        vamb --outdir {output} \
             --fasta {input.catalogue} \
             --bamfiles {input.bamfiles} \
             -o {params.binsplit_sep} \
             --minfasta {params.minfasta} &> {log}
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
        "envs/bash.yaml"
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
        "envs/bash.yaml"
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
