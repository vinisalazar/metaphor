# Metaphor
## Metagenomics Pipeline for Short Reads

Metaphor is a Snakemake-based workflow for analysis of metagenomics short reads data. It includes the following steps:
- Quality control (with [FastQC](https://github.com/s-andrews/FastQC/), [Cutadapt](https://github.com/marcelm/cutadapt))
- Assembly (with [Megahit](https://github.com/voutcn/megahit), and evaluation is done with [MetaQUAST](https://github.com/ablab/quast))
- Read mapping (with [Minimap2](https://github.com/lh3/minimap2), [Samtools](https://github.com/samtools/samtools))
- Binning (binners include [Vamb](https://github.com/RasmussenLab/vamb/), [MetaBAT](https://bitbucket.org/berkeleylab/metabat), [CONCOCT](https://github.com/BinPro/CONCOCT), [GraphBin](https://github.com/Vini2/GraphBin); refinement is done with [DAS Tool](https://github.com/cmks/DAS_Tool))
- Annotation (with [Prodigal](https://github.com/hyattpd/Prodigal), [Diamond](https://github.com/bbuchfink/diamond), and the [NCBI COG database](https://www.ncbi.nlm.nih.gov/research/cog-project/))
- Postprocessing (with custom scripts)

Metaphor aims to be concise, portable, and sustainable. It only includes third-party software that is properly packaged and easily installable. Metaphor originally started as Snakemake port of [MetaGenePipe](https://gitlab.unimelb.edu.au/bshaban/metaGenePipe/), developed by the [Melbourne Data Analytics Platform](https://mdap.unimelb.edu.au/), specifically by Bobbie Shaban *et al*, but  has since seen significant modifications (including the name change).

If you have any questions regarding Metaphor, please don't hesitate to open an issue.

<!-- 
TO-DO:
### Installation

### Usage

### Documentation 
-->
