
# Configuration

This page explains each value of Metaphor's config settings, that is, the values defined in the config YAML file.


**TOP LEVEL**

These settings are valid for all steps in the workflow.

**`samples:`** `samples.csv`    


**`mb_per_core:`** `2048`   How many MBs of RAM to use per processor. For more information about this, read the "Advanced" page of the docs.  


**QC**

**`cutadapt:`**   Settings for the Cutadapt cutadapt tool.  
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`phred:`** `33`    
&nbsp;&nbsp;&nbsp;**`minimum_length:`** `50`    
&nbsp;&nbsp;&nbsp;**`quality_cutoff:`** `30`    
&nbsp;&nbsp;&nbsp;**`clip_r5:`** `10`    
&nbsp;&nbsp;&nbsp;**`clip_r3:`** `5`    


**`merge_reads:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    


**`fastqc:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    


**`multiqc:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    


**ASSEMBLY**

**`coassembly:`** `False`   Whether to perform coassembly (also known as pooled assembly). If this is true, all samples are merged together and assembled into a single file of contigs.  


**`megahit:`**    
&nbsp;&nbsp;&nbsp;**`preset:`** `"meta-large"`    
&nbsp;&nbsp;&nbsp;**`cleanup:`** `True`    


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
&nbsp;&nbsp;&nbsp;**`args:`** `"--metagenome --quiet --force"`    
&nbsp;&nbsp;&nbsp;**`kingdom:`** `"Bacteria"`    


**`diamond:`**    
&nbsp;&nbsp;&nbsp;**`db:`** `"data/COG2020/cog-20.dmnd"`   Will try to create from db_source if it doesn't exist.  
&nbsp;&nbsp;&nbsp;**`db_source:`** `"data/COG2020/cog-20.fa.gz"`    
&nbsp;&nbsp;&nbsp;**`output_type:`** `6`    
&nbsp;&nbsp;&nbsp;**`output_format:`** `"qseqid sseqid stitle evalue bitscore staxids sscinames"`    


**`cog_functional_parser:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`db:`** `"data/COG2020"`    


**`lineage_parser:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`taxonmap:`** `"data/COG2020/cog-20.taxonmap.tsv"`    
&nbsp;&nbsp;&nbsp;**`rankedlineage:`** `"data/taxonomy/rankedlineage.dmp"`    
&nbsp;&nbsp;&nbsp;**`names:`** `"data/taxonomy/names.dmp"`   Path of names file of NCBI Taxonomy  
&nbsp;&nbsp;&nbsp;**`nodes:`** `"data/taxonomy/nodes.dmp"`   Path of nodes file of NCBI Taxonomy  
&nbsp;&nbsp;&nbsp;**`download_url:`** `"https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"`   URL to download NCBI Taxonomy database  


**`plot_cog_functional:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`filter_categories:`** `True`    
&nbsp;&nbsp;&nbsp;**`categories_cutoff:`** `0.01`    


**`plot_taxonomies:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`tax_cutoff:`** `0.000001`   1e-6  


**BINNING**

**`vamb:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`minfasta:`** `10000`    


**`metabat2:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`seed:`** `0`    
&nbsp;&nbsp;&nbsp;**`preffix:`** `"bin"`   Preffix of each bin, e.g. bin.1.fa, bin.2.fa, etc.  


**`concoct:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    


**`das_tool:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`score_threshold:`** `0.5`    


**POSTPROCESSING**

**`postprocessing:`**    
&nbsp;&nbsp;&nbsp;**`activate:`** `True`    
&nbsp;&nbsp;&nbsp;**`runtime_unit:`** `"m"`    
&nbsp;&nbsp;&nbsp;**`runtime_cutoff:`** `5`    
&nbsp;&nbsp;&nbsp;**`memory_unit:`** `"max_vms"`    
&nbsp;&nbsp;&nbsp;**`memory_cutoff:`** `1`    
&nbsp;&nbsp;&nbsp;**`memory_gb:`** `True`    


