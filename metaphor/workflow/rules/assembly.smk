"""
assembly.smk

    Assemble short reads with MEGAHIT. Additionally, evalute reads with MetaQuast if a reference is provided.

Assembly rules:
    - concatenate_merged_reads: prepare MEGAHIT input if there coassembly is True in config file
    - megahit: assemble preprocessed reads with MEGAHIT
    - metaquast: evaluate assembly results with MetaQuast
"""


rule concatenate_merged_reads:
    input:
        R1=expand("output/qc/merged/{sample}_R1.fq.gz", sample=sample_IDs),
        R2=expand("output/qc/merged/{sample}_R2.fq.gz", sample=sample_IDs),
    output:
        R1_concat="output/qc/merged/all_samples_R1.fq.gz",
        R2_concat="output/qc/merged/all_samples_R2.fq.gz",
    log:
        "output/logs/assembly/concatenate_merged_reads.log",
    benchmark:
        "output/benchmarks/assembly/concatenate_merged_reads.txt"
    conda:
        "../envs/utils.yaml"
    shell:
        """
        {{ cat {input.R1} > {output.R1_concat} ; }} > {log}
        {{ cat {input.R2} > {output.R2_concat} ; }} >> {log}
        """


rule megahit:
    input:
        fastq1=lambda w: get_fastq_groups(w, "R1"),
        fastq2=lambda w: get_fastq_groups(w, "R2"),
    output:
        contigs=get_contigs_input(),
    params:
        # Turn 'remove_intermediates' on/off in config['megahit']
        fastq1=lambda w, input: ",".join(input.fastq1),
        fastq2=lambda w, input: ",".join(input.fastq2),
        out_dir=lambda w, output: get_parent(get_parent(output.contigs)),  # this is equivalent to "{output}/megahit"
        min_contig_len=config["megahit"]["min_contig_len"],
        k_list="21,29,39,59,79,99,119,141",
        preset=config["megahit"]["preset"],
        cleanup=lambda w, output: cleanup_rule(
            "megahit", Path(output.contigs).parent.joinpath("intermediate_contigs")
        ),
        sample=lambda w: w.group,
    threads: get_threads_per_task_size("big")
    wildcard_constraints:
        group="|".join(group_names),
    resources:
        mem_mb=get_mb_per_cores,
    log:
        get_group_benchmark_or_log("log", "assembly", "megahit"),
    benchmark:
        get_group_benchmark_or_log("benchmark", "assembly", "megahit")
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        # MegaHit has no --force flag, so we must remove the created directory prior to running
        rm -rf {params.out_dir}/{params.sample}

        megahit -1 {params.fastq1} -2 {params.fastq2}         \
                -o {params.out_dir}/{params.sample}         \
                --presets {params.preset}                   \
                --out-prefix {params.sample}                \
                --min-contig-len {params.min_contig_len}    \
                -t {threads}                                \
                --k-list {params.k_list} &> {log}

        {params.cleanup}
        """


rule assembly_report:
    """
    Get metrics for each assembly.
    """
    input:
        contigs=get_contigs_input(expand_=True),
    params:
        fastas=lambda w, input: " ".join(
            input.contigs
            if not isinstance(input.contigs, str)
            else [
                input.contigs,
            ]
        ),
    output:
        report(get_assembly_report("avg_length"), category="Assembly"),
        report(get_assembly_report("max_length"), category="Assembly"),
        report(get_assembly_report("median_length"), category="Assembly"),
        report(get_assembly_report("n_bp"), category="Assembly"),
        report(get_assembly_report("n_contigs"), category="Assembly"),
        report(get_assembly_report("n50"), category="Assembly"),
        assembly_report=get_assembly_report(),
    log:
        "output/logs/assembly/assembly_report.log",
    benchmark:
        "output/benchmarks/assembly/assembly_report.txt"
    conda:
        "../envs/utils.yaml"
    script:
        "../scripts/assembly_report.py"


rule metaquast:
    input:
        contigs=get_contigs_input(),
        reference=get_metaquast_reference,
    output:
        outfile=get_group_or_sample_file(
            "assembly", "metaquast", suffix="report.html", add_sample_to_suffix=False
        ),
    params:
        mincontig=500,
        outdir=lambda w, output: str(Path(output.outfile).parent),
        extra_params="--fragmented --unique-mapping --no-icarus --no-plots --no-gc --no-sv",
    threads: get_threads_per_task_size("big")
    resources:
        mem_mb=get_mb_per_cores,
    log:
        get_group_benchmark_or_log("log", "assembly", "metaquast"),
    benchmark:
        get_group_benchmark_or_log("benchmark", "assembly", "metaquast")
    conda:
        "../envs/quast.yaml"
    shell:
        """
        metaquast.py -t {threads}               \
                     -o {params.outdir}         \
                     -m {params.mincontig}      \
                     -r {input.reference}       \
                     {params.extra_params}      \
                     {input.contigs} &> {log}
        """
