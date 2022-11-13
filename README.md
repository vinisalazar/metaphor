# Metaphor
## Metagenomics Pipeline for sHOrt Reads

![https://metaphor-workflow.readthedocs.io/](https://img.shields.io/readthedocs/metaphor-workflow?color=g) ![https://github.com/vinisalazar/metaphor/tags](https://img.shields.io/github/v/tag/vinisalazar/metaphor?color=g&label=release)

[![Version](https://anaconda.org/bioconda/metaphor/badges/version.svg)](https://anaconda.org/bioconda/metaphor)
<!-- [![Bioconda](https://img.shields.io/conda/dn/bioconda/metaphor.svg?label=Bioconda )](https://anaconda.org/bioconda/metaphor) -->

Metaphor is a Snakemake-based workflow for analysis of metagenomics short reads data. It includes the following steps:
- Quality control (with [FastQC](https://github.com/s-andrews/FastQC/), [Cutadapt](https://github.com/marcelm/cutadapt))
- Assembly (with [Megahit](https://github.com/voutcn/megahit), and evaluation is done with [MetaQUAST](https://github.com/ablab/quast))
- Read mapping (with [Minimap2](https://github.com/lh3/minimap2), [Samtools](https://github.com/samtools/samtools))
- Binning (binners include [Vamb](https://github.com/RasmussenLab/vamb/), [MetaBAT](https://bitbucket.org/berkeleylab/metabat), [CONCOCT](https://github.com/BinPro/CONCOCT)<!--, [GraphBin](https://github.com/Vini2/GraphBin)-->; refinement is done with [DAS Tool](https://github.com/cmks/DAS_Tool))
- Annotation (with [Prodigal](https://github.com/hyattpd/Prodigal), [Diamond](https://github.com/bbuchfink/diamond), and the [NCBI COG database](https://www.ncbi.nlm.nih.gov/research/cog-project/))
- Postprocessing (with custom scripts)

Please cite these software if you use Metaphor. The bib files are located [here](./metaphor/workflow/bibs/) for your convenience.
Metaphor will support automatic citation in the future.

Metaphor aims to be concise, portable, and sustainable. It only includes third-party software that is properly packaged and easily installable.
If you have any questions regarding Metaphor, please don't hesitate to open an issue.

Check out our [documentation!](https://metaphor-workflow.readthedocs.io)

### Installation
The first thing that you need to install Metaphor is [conda](https://docs.conda.io/). To install it, please follow their [user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

Once you have conda, we highly recommend that you use [mamba](https://mamba.readthedocs.io/en/latest/installation.html) for installing Metaphor. If mamba is not available, replace all `mamba` commands for `conda`.

To install, either create a new environment or install it in your preferred environment:
```bash
$ mamba create -n metaphor metaphor -c conda-forge -c bioconda
$ conda activate metaphor
```

You should see the `(metaphor)` indicator next to your prompt.

### Testing
After installing, check if your Metaphor installation works:

```bash
# Check that the `metaphor` command works
$ metaphor -h

# You can see available options for testing Metaphor with:
$ metaphor test -h

# To test Metaphor (follow the screen prompts)
$ metaphor test
```

Testing may take a long time (a couple of hours), so please be patient. After testing, you can run the workflow on the same directory, this way you won't have to install the conda dependencies all over again.

### Usage
To run Metaphor on your data, we recommend that you create a configuration profile specific to your needs, and then run Metaphor on your directory of FASTQ files:

```bash
# Run this and follow the screen prompts
$ metaphor config settings

# To execute with your config simply type 
$ metaphor execute -i path/to/directory/of/fastq
```

If you receive any errors, feel free to open an issue describing your problem.


##### DISCLAIMER
Metaphor is a derivative work of [MetaGenePipe](https://gitlab.unimelb.edu.au/bshaban/metaGenePipe/), originally released under the
Apache 2.0 license, developed by [Bobbie Shaban](https://gitlab.unimelb.edu.au/bshaban), Mar Quiroga, Robert Turnbull
and Edoardo Tescari at Melbourne Data Analytics Platform ([MDAP](https://mdap.unimelb.edu.au/)) at the
University of Melbourne.
[Link to MetaGenePipe](https://gitlab.unimelb.edu.au/bshaban/metaGenePipe/). For more information, please see the [license file](./LICENSE.md).

