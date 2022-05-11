# Output

Metaphor's output is organised in accordance with each module:

* Quality control (QC)
* Assembly
* Annotation
* Mapping
* Binning
* Postprocessing

Additionally, there are the `logs` and `benchmarks` directory which contain the log and benchmark files of each rule.
Log files contain the rules' standard error (stderr), so they are useful for debugging in case your workflow errors out.
Benchmark files, on the other hand, contain the output of the [`psutil`](https://psutil.readthedocs.io/en/latest/)
command, which is the memory consumption and runtime of each rule.

The output directory of each module contains subdirectories named after each rule. These rule directories may contain
yet another level of subdirectories named after each sample. So, the structure of Metaphor's output in general follows
this pattern:

```{code-block} console
(metaphor)$ tree output/
output/
└── module_1/
    └── rule_a/
        └── file_aa
        └── file_ab
└── module_2/
    └── rule_b/
        └── sample_1/
            └── file_ba
            └── file_bb
        └── sample_2/
```

And so on.

## QC

The QC module contains four directories and the `multiqc.html` report:

* `cutadapt`: output of the cutadapt rule, which removes the adapters and performs quality trimming on the raw
reads.
* `fastqc`: output of the FastQC rule. Has subdirectories for raw, trimmed, and merged reads.
* `merged`: merged reads after QC.
* `multiqc_data`: data for the `multiqc.html` report.

## Assembly

The assembly module contains two directories:

* `assembly_report`: report with assembly metrics for all samples, along with plots.
* `megahit`: results of the MegaHIT assembly. Contains subdirectories for each sample.

## Annotation

This module contains three directories:

* `cog`: final output of the annotation with the NCBI COG database. Generated from the DIAMOND tables along with the
COG reference files. This contains subdirectories for each sample with the taxonomy and annotation tables for each
sample, and two extra subdirectories:
    * `tables`: the concatenated (with all samples) tables of taxonomy and functional annotation.
    * `plots`: the plots generated from the concatenated tables.
* `diamond`: output tables of the DIAMOND annotation. These tables are similar to
[BLAST tabular](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6) format.
* `prodigal`: output of gene prediction with Prodigal. Contains subdirectories for each sample. Each subdirectory
contains 3-4 files:
    * `{sample}_genbank.gbk`: the [GenBank Flat File format](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html).
    * `{sample}_genes.fna`: predicted coding sequences as nucleotides.
    * `{sample}_proteins.faa`: predicted coding sequences as amino acids
    * `{sample}_scores.cds` (optional):  is the tabular format with the scores for all possible genes.

## Mapping

This module contains the files that map the original reads back to the contigs. These files are required for the
binning process:
* `bam`: directory with BAM (four different kinds) files for all samples.
* `catalogue.fna.gz`: concatenated contigs for all samples.
* `catalogue.mmi`: index file of the concatenated contigs.
* `bam_contig_depths.txt`: coverage of each contig calculated from BAM files.

## Binning

This module contains one directory for each of the binners:
* `metabat2`: bins are inside this directory as `.fa` files.
* `vamb`: contains a `bins` directory.
* `concoct`: contains a `fasta_bins` directory.
* `DAS_tool`: contains the refined bins inside the `DASTool_bins` directory, and the table with quality score for each
bin is named as `DASTool_summary.tsv`

## Postprocessing

The postprocessing module contains four different plots:
* `memory_barplot_errorbar.png`
* `memory_barplot_sum.png`
* `runtime_barplot_errorbar.png`
* `runtime_barplot_sum.png`

These plots show the total (`_sum`) and average (`_errorbar`) runtime and memory consumption for each rule. You can use
them to identify computational bottlenecks in your analyses, that is, if any rule in particular is taking more time than
it should.

