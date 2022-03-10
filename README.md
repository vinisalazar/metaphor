# Metaphor
## Metagenomics Pipeline for Short Reads

Metaphor is a Snakemake-based workflow for analysis of metagenomics short reads data. It includes the following steps:
- Quality control (with [FastQC](https://github.com/s-andrews/FastQC/), [Cutadapt](https://github.com/marcelm/cutadapt))
- Assembly (with [Megahit](https://github.com/voutcn/megahit), and evaluation is done with [MetaQUAST](https://github.com/ablab/quast))
- Read mapping (with [Minimap2](https://github.com/lh3/minimap2), [Samtools](https://github.com/samtools/samtools))
- Binning (binners include [Vamb](https://github.com/RasmussenLab/vamb/), [MetaBAT](https://bitbucket.org/berkeleylab/metabat), [CONCOCT](https://github.com/BinPro/CONCOCT)<!--, [GraphBin](https://github.com/Vini2/GraphBin)-->; refinement is done with [DAS Tool](https://github.com/cmks/DAS_Tool))
- Annotation (with [Prodigal](https://github.com/hyattpd/Prodigal), [Diamond](https://github.com/bbuchfink/diamond), and the [NCBI COG database](https://www.ncbi.nlm.nih.gov/research/cog-project/))
- Postprocessing (with custom scripts)

Metaphor aims to be concise, portable, and sustainable. It only includes third-party software that is properly packaged and easily installable. Metaphor originally started as Snakemake port of [MetaGenePipe](https://gitlab.unimelb.edu.au/bshaban/metaGenePipe/), developed by by Bobbie Shaban *et al* at the [Melbourne Data Analytics Platform](https://mdap.unimelb.edu.au/), but has since seen significant modifications (including the name change).

If you have any questions regarding Metaphor, please don't hesitate to open an issue.

### Installation
The first thing that you need to install Metaphor is [conda](https://docs.conda.io/). To install it, please follow their [user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

Once you have conda, we highly recommend that you use [mamba](https://mamba.readthedocs.io/en/latest/installation.html) for installing Metaphor. If mamba is not available, replace all `mamba` commands for `conda`.

Run the following commands to install Metaphor from source.
```bash

# Install mamba (if you haven't already)
$ conda install mamba -n base -c conda-forge

# Copy the code to your machine
$ git clone https://github.com/vinisalazar/metaphor && cd metaphor

# Create the Metaphor environment
$ mamba env create -n metaphor -f environment.yaml && conda activate metaphor

# Install Metaphor with pip
pip install .

# Check that the installation works
metaphor -h

# Test Metaphor (follow the screen prompts)
metaphor test
```

Testing may take a long time (a couple of hours), so please be patient.

### Usage
To run Metaphor on your data, we recommend that you create a configuration profile specific to your needs, and then run Metaphor on your directory of FASTQ files:

```bash
# Run this and follow the screen prompts
$ metaphor config settings

# To execute with your config simply type 
$ metaphor execute -i path/to/directory/of/fastq
```

If you receive any errors, feel free to open an issue describing your problem.

<!-- 
TO-DO:
### Documentation 
-->
