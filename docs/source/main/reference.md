
# Reference

This is an automatically generated list of all supported rules, their docstrings, and command. At the start of each 
workflow run a list is printed of which rules will be run. And while the workflow is running it prints which rules are
being started and finished. This page is here to give an explanation to the user about what each rule does, and for
developers to find what is, and isn't yet supported. Not all Metaphor rules are listed here, only the ones with a
`shell` directive. Rules with `script` or `wrapper` directives are not included. To see all rules in Metaphor, 
please refer to the [workflow source code](https://github.com/vinisalazar/metaphor/tree/main/workflow).
## qc.smk
**cutadapt_pipe**

Pipe reads into cutadapt.

```
cat {input} > {output} 2> {log}
```

**cutadapt_pe**

Trim paired end reads with cutadapt.

**merge_fastqs**

Concatenate paired-end reads from different units.

```
cat {input} > {output} 2> {log}
```

**fastqc_raw**

**host_removal**

```
{{ minimap2 -t {threads}           \
         -ax {params.preset}       \
         {input.reference}         \
         {input.fastqs} ; }}        2>> {log}   | 
{{ samtools view -buSh -f 4 ; }}    2>> {log}   |
{{ samtools fastq - ; }}            2>> {log}   |
{{ gzip -c - > {output} ; }}        2>> {log}
```

## assembly.smk
**concatenate_merged_reads**

```
{{ cat {input.R1} > {output.R1_concat} ; }} > {log}
{{ cat {input.R2} > {output.R2_concat} ; }} >> {log}
```

**megahit**

```
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
```

**assembly_report**

Get metrics for each assembly.

**metaquast**

```
metaquast.py -t {threads}               \
             -o {params.outdir}         \
             -m {params.mincontig}      \
             -r {input.reference}       \
             {params.extra_params}      \
             {input.contigs} &> {log}
```

## annotation.smk
**prodigal**

```
prodigal {params.quiet}         \
         -p {params.mode}       \
         -i {input}             \
         -a {output.proteins}   \
         -o {output.genbank}    \
         {params.genes}         \
         {params.scores} &> {log}
```

**prokka**

```
# Get kingdom from bin eval file
kingdom=$({params.kingdom_cmd} | cut -f 5)
kingdom=$(echo $kingdom | head -c 1 | tr '[a-z]' '[A-Z]'; echo $kingdom | tail -c +2)

prokka --outdir {params.outdir}     \
       --kingdom $kingdom           \
       --cpus {threads}             \
       --prefix {wildcards.bin}     \
       {params.args}                \
       {input.genome_bin} &> {log}
```

**download_COG_database**

```
for file in {output}; do
    wget https://ftp.ncbi.nih.gov/pub/COG/COG2020/data/$(basename $file) -O $file 2>> {log};
done
```

**generate_COG_taxonmap**

**diamond**

```
echo {params.output_format} | sed -e 's/ /\t/g' > {output.fname}
{{ diamond blastp -q {input.fname_fasta}                \
           -p {threads}                                 \
           -d {input.fname_db}                          \
           -f {params.output_type}                      \
           {params.output_format}                       \
           {params.extra}                               \
           >> {output.fname} ; }} &> {log}
```

## mapping.smk
**concatenate_contigs**

```
concatenate.py -m {params.sequence_length_cutoff} {output} {input} &> {log}
```

**decompress_catalogue**

```
pigz -d -f -p {threads} -k {input.catalogue_gz} &> {log}
```

**concatenate_proteins**

Used by DAS_Tool (skips the Prodigal run).

```
cat {input} > {output}
```

**create_contigs_index**

```
minimap2 -d {output} {input} &> {log}
```

**create_genes_index**

```
minimap2 -d {output} {input} &> {log}
```

**map_reads**

```
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
```

## binning.smk
**vamb**

```
rm -rf {params.outdir}

vamb --outdir {params.outdir}           \
     --fasta {input.catalogue}          \
     --jgi {input.bam_contig_depths}    \
     -p {threads}                       \
     -o {params.binsplit_sep}           \
     -t {params.batchsize}              \
     --minfasta {params.minfasta} &> {log}

{{ awk -v OFS='\t' '{{ print $2, $1 }}' {output.clusters} |  \
sed "s/$(echo '\t')/$(echo '\t')vamb./g" >          \
{output.scaffolds2bin} ; }} >> {log} 2>&1
```

**metabat2**

```
rm -rf {output.outdir} && mkdir {output.outdir}

metabat2 -i {input.contigs}             \
         -a {input.depths}              \
         -m {params.minContig}          \
         -t {threads}                   \
         --seed {params.seed}           \
         --saveCls                      \
         -o {params.outfile} &> {log}

sed "s/$(echo '\t')/$(echo '\t')metabat2./g" {params.outfile} > {output.scaffolds2bin}
```

**concoct**

```
 
rm -rf {output.outdir}
mkdir {output.outdir} 

{{ cut_up_fasta.py {input.catalogue}                        \
                   -c {params.contig_size}                  \
                   -o 0                                     \
                   -b {params.bed}                          \
                   --merge_last                             \
                   > {params.contigs}  ; }} 2>> {log}

{{ concoct_coverage_table.py {params.bed}                   \
                             {input.bams}                   \
                             > {params.coverage_table} ; }} 2>> {log}

{{ concoct --composition_file {params.contigs}              \
           --coverage_file {params.coverage_table}          \
           -b {output.outdir}                               \
           -t {threads}  ; }} 2>> {log}

{{ merge_cutup_clustering.py {params.clustering_gt}         \
                             > {params.clustering_merged}  ; }} 2>> {log}

mkdir {params.fasta_bins}

{{ extract_fasta_bins.py {input.catalogue}                  \
                         {params.clustering_merged}         \
                         --output_path {params.fasta_bins} ; }} 2>> {log}

sed "s/,/$(echo '\t')concoct./g" {params.clustering_merged} | tail -n +2 > {output.scaffolds2bin}
```

**DAS_tool**

Refine bins assembled with one or more binners.

```
DAS_Tool -i {params.fmt_scaffolds2bin}                      \
         -l {params.binners}                                \
         -c {input.contigs}                                 \
         -o {params.outpreffix}                             \
         -p {input.proteins}                                \
         --score_threshold {params.score_threshold}         \
         --search_engine diamond                            \
         --write_bins                                       \
         --write_bins_eval                                  \
         {params.extra}                                     \
         --threads {threads} &> {log}
```

**bins_report**

Plots binning metrics generated by DAS Tool.


---

**Disclaimer**

This page was generated with a script adapted from the 
[seq2science repository](https://github.com/vanheeringen-lab/seq2science).

MIT License

Copyright (c) 2019 Maarten-vd-Sande (vanheeringen-lab)

For the full license, please see the
[script source code](https://github.com/vinisalazar/metaphor/blob/main/docs/scripts/rule_description.py).
