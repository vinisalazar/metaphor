
# Configuration

This page explains each value of Metaphor's config settings, that is, the values defined in the config YAML file.


**TOP LEVEL**

These settings are valid for all steps in the workflow.

**`samples:`** `samples.csv`    


**QC**

**`fastp:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`length_required:`** `50`    
&nbsp;&nbsp;&nbsp;**`cut_mean_quality:`** `30`    
&nbsp;&nbsp;&nbsp;**`extra:`** `"--detect_adapter_for_pe"`    


**`merge_reads:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    


**`host_removal:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `False`    
&nbsp;&nbsp;&nbsp;**`reference:`** `""`    


**`fastqc:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    


**`multiqc:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    


**ASSEMBLY**

**`coassembly:`** `False`   Whether to perform coassembly (also known as pooled assembly). If this is true, all samples are merged together and assembled into a single file of contigs.  


**`megahit:`**    
&nbsp;&nbsp;&nbsp;**`preset:`** `"meta-large"`    
&nbsp;&nbsp;&nbsp;**`cleanup:`** `True`    
&nbsp;&nbsp;&nbsp;**`min_contig_len:`** `1000`    


**`rename_contigs:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`   Whether to rename contigs so contigs and mapping files (.bam) can be imported into Anvi'o. We suggest you keep this on.  
&nbsp;&nbsp;&nbsp;**`awk_command:`** `awk '/^>/{{gsub(" |\\\\.|=", "_", $0); print $0; next}}{{print}}' {input} > {output}`   This is to prevent errors with the Snakemake --lint command. Don't change it unless you know what you're doing.  


**`metaquast:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `False`    
&nbsp;&nbsp;&nbsp;**`coassembly_reference:`** `""`   Reference FASTA file for Metaquast to use as reference. Only required if `coassembly` is True.  


**ANNOTATION**

**`prodigal:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`mode:`** `"meta"`    
&nbsp;&nbsp;&nbsp;**`quiet:`** `True`    
&nbsp;&nbsp;&nbsp;**`genes:`** `False`    
&nbsp;&nbsp;&nbsp;**`scores:`** `False`    


**`prokka:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `False`    
&nbsp;&nbsp;&nbsp;**`args:`** `"--quiet --force"`    


**`diamond:`**    
&nbsp;&nbsp;&nbsp;**`db:`** `"COG2020/cog-20.dmnd"`   Will try to create from db_source if it doesn't exist.  
&nbsp;&nbsp;&nbsp;**`db_source:`** `"COG2020/cog-20.fa.gz"`    
&nbsp;&nbsp;&nbsp;**`output_type:`** `6`    
&nbsp;&nbsp;&nbsp;**`output_format:`** `"qseqid sseqid stitle evalue bitscore staxids sscinames"`    


**`cog_functional_parser:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`db:`** `"COG2020"`    


**`lineage_parser:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`taxonmap:`** `"COG2020/cog-20.taxonmap.tsv"`    
&nbsp;&nbsp;&nbsp;**`rankedlineage:`** `"taxonomy/rankedlineage.dmp"`    
&nbsp;&nbsp;&nbsp;**`names:`** `"taxonomy/names.dmp"`   Path of names file of NCBI Taxonomy  
&nbsp;&nbsp;&nbsp;**`nodes:`** `"taxonomy/nodes.dmp"`   Path of nodes file of NCBI Taxonomy  
&nbsp;&nbsp;&nbsp;**`download_url:`** `"https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"`   URL to download NCBI Taxonomy database  


**`plot_cog_functional:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`filter_categories:`** `True`    
&nbsp;&nbsp;&nbsp;**`categories_cutoff:`** `0.01`   Remove categories with mean abundance across samples smaller than this value  


**`plot_taxonomies:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`tax_cutoff:`** `20`   Only show the N most abundant taxa for any rank. Leave as 0 for no filtering. Low abundance taxa will be grouped as 'Low abundance'.  
&nbsp;&nbsp;&nbsp;**`colormap:`** `"tab20c"`   Which matplotlib colormap to use  


**BINNING**

**`cobinning:`** `True`   Whether to perform cobinning. When this is true, only one binning group will be used. If False, samples will be binned according to their 'group' column.  


**`vamb:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`minfasta:`** `100000`    
&nbsp;&nbsp;&nbsp;**`batchsize:`** `256`    


**`metabat2:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`seed:`** `0`    
&nbsp;&nbsp;&nbsp;**`preffix:`** `"bin"`   Preffix of each bin, e.g. bin.1.fa, bin.2.fa, etc.  


**`concoct:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    


**`das_tool:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`score_threshold:`** `0.5`    
&nbsp;&nbsp;&nbsp;**`bins_report:`** `True`    


**POSTPROCESSING**

**`postprocessing:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`runtime_unit:`** `"m"`    
&nbsp;&nbsp;&nbsp;**`runtime_cutoff:`** `5`    
&nbsp;&nbsp;&nbsp;**`memory_unit:`** `"max_vms"`    
&nbsp;&nbsp;&nbsp;**`memory_cutoff:`** `1`    
&nbsp;&nbsp;&nbsp;**`memory_gb:`** `True`    


