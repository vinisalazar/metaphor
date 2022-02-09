"""
Assembly rules:
    - concatenate_merged_reads: prepare MegaHIT input if there coassembly is True in config file
    - megahit: assemble preprocessed reads with Megahit
    - megahit_coassembly: perform coassembly of pool of reads of all samples with Megahit
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
        "../envs/bash.yaml"
    shell:
        """
        {{ cat {input.R1} > {output.R1_concat} ; }} > {log}
        {{ cat {input.R2} > {output.R2_concat} ; }} >> {log}
        """


rule megahit:
    input:
        fastq1="output/qc/merged/{sample}_R1.fq.gz"
        if not config["coassembly"]
        else "output/qc/merged/all_samples_R1.fq.gz",
        fastq2="output/qc/merged/{sample}_R2.fq.gz"
        if not config["coassembly"]
        else "output/qc/merged/all_samples_R2.fq.gz",
    output:
        contigs=get_contigs_input(),
    params:
        # Turn 'remove_intermediates' on/off in config['megahit']
        out_dir=lambda w, output: get_parent(get_parent(output.contigs)),  # this is equivalent to "{output}/megahit"
        min_contig_len=200,
        k_list="21,29,39,59,79,99,119,141",
        preset=config["megahit"]["preset"],
        remove_intermediate=lambda w, output: cleanup_megahit(output.contigs),
        sample=lambda w: w.sample if getattr(w, 'sample', None) else 'coassembly'
    threads: round(workflow.cores * 0.75)
    log:
        get_coassembly_benchmark_or_log("log", "assembly", "megahit"),
    benchmark:
        get_coassembly_benchmark_or_log("benchmark", "assembly", "megahit")
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        # MegaHit has no --force flag, so we must remove the created directory prior to running
        rm -rf {params.out_dir}/{params.sample}

        megahit -1 {input.fastq1} -2 {input.fastq2}         \
                -o {params.out_dir}/{params.sample}         \
                --presets {params.preset}                   \
                --out-prefix {params.sample}                \
                --min-contig-len {params.min_contig_len}    \
                -t {threads}                                \
                --k-list {params.k_list} &> {log}

        {params.remove_intermediate}
        """


rule assembly_report:
    input:
        contigs=get_contigs_input(expand_=True),
    params:
        fastas=lambda w, input: " ".join(input.contigs if not isinstance(input.contigs, str) else [input.contigs, ]),
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
        "../envs/bash.yaml"
    script:
        "../scripts/assembly_report.py"


rule metaquast:
    input:
        contigs=get_contigs_input(),
        reference=get_metaquast_reference,
    output:
        outfile=get_coassembly_or_sample_file(
            "assembly", "metaquast", suffix="report.html", add_sample_to_suffix=False
        ),
    params:
        mincontig=500,
        outdir=lambda w, output: str(Path(output.outfile).parent),
        extra_params="--fragmented --no-icarus --no-plots --no-gc --no-sv",
    threads: round(workflow.cores * 0.75)
    resources:
        mem_mb=get_mem_mb,
    log:
        get_coassembly_benchmark_or_log("log", "assembly", "metaquast"),
    benchmark:
        get_coassembly_benchmark_or_log("benchmark", "assembly", "metaquast")
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
