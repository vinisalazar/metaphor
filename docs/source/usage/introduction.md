# Introduction
<!-- 
- Install
- Structure
- Input
- Configuring
- Output
 -->
Metaphor - Metagenomic Pipeline for Short Reads - is a [Snakemake](https://snakemake.readthedocs.io/)-based workflow for assembly and binning of metagenomes. It also performs quality control, and basic functional and taxonomic annotation of the assembled contigs, using the [NCBI COG](https://www.ncbi.nlm.nih.gov/research/cog/) database.

Metaphor is designed to be lightweight, flexible, and sustainable. That means that it strives to produce the desired output with the minimal amount of dependencies and overhead, it is suitable for a wide range of use cases, and it is easy to maintain and modify. Metagenomic analyses are usually quite complex, with numerous steps and dependencies. Metaphor's goal is to simplify that, and to provide users with a final output that enables exploratory data analysis and that can be "plugged" into other downstream pipelines, such as phylogenomics or advanced genome annotation pipelines.

## Installation

To install Metaphor, the only thing you'll really need is [conda](https://docs.conda.io/). To install it, please follow their [user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

Once you have conda, we highly recommend that you use [mamba](https://mamba.readthedocs.io/en/latest/installation.html) for installing Metaphor. If mamba is not available, replace all the following `mamba` commands for `conda`.

The following command creates a new environment containing Metaphor:

```{code-block} console
$ mamba create -n metaphor -c conda-forge -c bioconda metaphor
```

The `-n` flag indicates the environment name, which can be replaced for whatever you want. The `-c` flags indicate that we will install the dependencies for Metaphor from the [conda-forge](https://conda-forge.org/) and [Bioconda](https://bioconda.github.io/) repositories. These repositories contain the most up-to-date software necessary to run Metaphor. Lastly, the `metaphor` argument at the end indicates that we wish to install Metaphor in our environment.

Once you have it installed, activate the newly-created environment environment and check if Metaphor is found:

```{code-block} console
$ conda activate metaphor

$ metaphor -h
```

If you see a help prompt, that means Metaphor is correctly installed.

Let's move on to testing our installation.



% Metaphor is packaged as a Python application that contains  
