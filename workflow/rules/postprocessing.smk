"""
Post-processing module.

After all third-party programs are done running, this module performs additional steps such
as visualizations, removal of unnecessary intermediary files, and so on.

Rules in this module are submitted as a single 'group', see
https://snakemake.readthedocs.io/en/stable/executing/grouping.html;

this way they are submitted to a single node and can be run in succession.
"""


rule concatenate_benchmarks:
    input:
        get_final_output(),
    output:
        outfile="output/benchmarks/processing_benchmarks.csv",
    params:
        benchmarks_dir=lambda w, output: str(Path(output.outfile).parent),
    group:
        "postprocessing"
    log:
        "output/logs/postprocessing/concatenate_benchmarks.log",
    benchmark:
        "output/benchmarks/postprocessing/concatenate_benchmarks.txt"
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/concatenate_benchmarks.py"


rule plot_benchmarks:
    input:
        benchmarks_df="output/benchmarks/processing_benchmarks.csv",
    output:
        runtime_barplot_sum="output/postprocessing/runtime_barplot_sum.png",
        runtime_barplot_errorbar="output/postprocessing/runtime_barplot_errorbar.png",
        memory_barplot_sum="output/postprocessing/memory_barplot_sum.png",
        memory_barplot_errorbar="output/postprocessing/memory_barplot_errorbar.png",
    params:
        time_unit=config["postprocessing"]["runtime_unit"],
        time_cutoff=config["postprocessing"]["runtime_cutoff"],
        memory_unit=config["postprocessing"]["memory_unit"],
        memory_cutoff=config["postprocessing"]["memory_cutoff"],
        gb="--gb" if config["postprocessing"]["memory_gb"] else "",
    group:
        "postprocessing"
    log:
        "output/logs/postprocessing/plot_benchmarks.log",
    benchmark:
        "output/benchmarks/postprocessing/plot_benchmarks.txt"
    conda:
        "../envs/bash.yaml"
    script:
        "../scripts/plot_benchmarks.py"


# TODO
# rule remove_intermediates:
#     input:
#     output:
#     params:
#     log:
#     benchmark:
#     conda:
#     shell:
# rule compress_outputs:
#     input:
#     output:
#     params:
#     log:
#     benchmark:
#     conda:
#     shell:
