# Tutorial
<!-- 
- Install
- Structure
- Input
- Configuring
- Output
 -->

Welcome to Metaphor's introductory tutorial. Here you will find instructions for installing, testing, configuring,
and running the workflow.

```{note}
Before you install the workflow, it is important to note that **Metaphor is only fully-supported for Linux.**
There is limited support for OS X, and there are known issues with Mac M1 systems. Metaphor is not supported for Windows,
although we have managed to run it using Windows Subsystem for Linux.

The reason for this is that many of Metaphor's dependencies only have Linux builds, so they may not work on different systems.
```

## Installation

To install Metaphor, the only thing you'll really need is [conda](https://docs.conda.io/). To install it, please follow
their [user guide](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

Once you have conda, we highly recommend that you use [mamba](https://mamba.readthedocs.io/en/latest/installation.html)
for installing Metaphor. If mamba is not available, replace all the following `mamba` commands for `conda`.


The following command creates a new environment containing Metaphor:

```{code-block} console
(base)$ mamba create -n metaphor -c conda-forge -c bioconda metaphor
```

The `-n` flag indicates the environment name, which can be replaced for whatever you want. The `-c` flags indicate that
we will install the dependencies for Metaphor from the [conda-forge](https://conda-forge.org/) and
[Bioconda](https://bioconda.github.io/) repositories. These repositories contain the most up-to-date software necessary
to run Metaphor. Lastly, the `metaphor` argument at the end indicates that we wish to install Metaphor in our
environment.

Once you have it installed, activate the newly-created environment environment and check if Metaphor is found:

```{code-block} console
(base)$ conda activate metaphor

(metaphor)$ metaphor -h
```

If you see a help prompt, that means Metaphor is correctly installed.

Let's move on to testing our installation. For that, let's create a separate directory that will be used as the test
directory.

From within that directory, run the test command.

```{code-block} console
(metaphor)$ mkdir tutorial-metaphor

(metaphor)$ cd tutorial-metaphor

(metaphor)$ metaphor test
```

Metaphor will start the download of test data (about 100MB). Once that's done, it will prompt you to start the workflow.
Enter `y` and press `Enter⏎` to confirm the prompt. Metaphor will start to download and install the workflow
dependencies.
Tests can take a while, but usually finish within two hours.

Once tests are finished, most of Metaphor's dependencies will be installed and can be reused for your analysis.

Now that Metaphor is installed and we know it's working as intended, we can try running it on our own data!

## Running

To run Metaphor on your own data, you will need a CSV or TSV file containing three columns: `sample_name`, `R1` and `R2`.
The first column, `sample_name`, must contain a short, unique identifier for each of your samples. The `R1` and `R2`
columns will contain, respectively, the **full path** to your forward and reverse paired-end reads. If you only have
single-end reads, you don't need to provide the `R1` column. This is what your file could look like:

| `sample_name` | `R1`                         | `R2`                         |
|---------------|------------------------------|------------------------------|
| `sample_1`    | `/path/to/sample_1_R1.fq.gz` | `/path/to/sample_1_R2.fq.gz` |
| `sample_2`    | `/path/to/sample_2_R1.fq.gz` | `/path/to/sample_2_R2.fq.gz` |
| `sample_3`    | `/path/to/sample_3_R1.fq.gz` | `/path/to/sample_3_R2.fq.gz` |

It doesn't matter if you use CSV or TSV, as long as you have a `sample_name` and an `R1` column.

You can generate that file manually or using the `metaphor config input` command. For example, if you have a `fastq/`
directory with the following files:

```{code-block} console
(metaphor)$ ls fastq/
fastq
└── sample_1_R1.fq.gz
└── sample_1_R2.fq.gz
└── sample_2_R1.fq.gz
└── sample_2_R2.fq.gz
└── sample_3_R1.fq.gz
└── sample_3_R2.fq.gz

(metaphor)$ metaphor config input -i fastq/
6 files detected:
        sample_1_R1.fq.gz
        sample_1_R2.fq.gz
        sample_2_R1.fq.gz
        sample_2_R2.fq.gz
        sample_3_R1.fq.gz
        sample_3_R2.fq.gz

3 samples detected, 2 files per sample.
sample_1:
                sample_1_R1.fq.gz
                sample_1_R2.fq.gz
sample_2:
                sample_2_R1.fq.gz
                sample_2_R2.fq.gz
sample_3:
                sample_3_R1.fq.gz
                sample_3_R2.fq.gz

Generated input table 'samples.csv'.
```

The generated `samples.csv` file will be the entry point for your Metaphor run. Always inspect the samples file after
creating it, to make sure that it's properly formatted.

We need one additional file to run Metaphor: a configuration file. While the `samples.csv` contains the path to our
input files, this second file will contain the settings used to run Metaphor. You can edit this file to your own
preferences to customise the parameters of your Metaphor run.

To create the configuration file, use the `config settings` command and follow the on-screen prompt. Press `Enter⏎` for
to use the default values. Input your values carefully! Specifically, the number of cores and amount of memory will
affect how long your run takes. You can read about each config value in the [Configuration](configuration.md) page.

```{code-block} console
(metaphor)$ metaphor config settings

[...]

Wrote YAML to 'metaphor_settings.yaml'.
```

Now you have everything to run Metaphor! You can use the `execute` command for that. If you run it in the test
directory, there's no need to install dependencies again.

```{code-block} console
(metaphor)$ metaphor execute
```

It's simple as that. When running the `execute` command without any flags, Metaphor will automatically pick up on the
`samples.csv` and `metaphor_settings.yaml` generated by the `config` commands.

## Next steps

Your Metaphor run can take a few hours up to several days to complete, depending on the size of your dataset and the
resources you've provided. Move on to the [Output](output.md) section to understand the output of your Metaphor run,
or to the [Advanced](advanced.md) section for an in-depth explanation of how Metaphor works.
