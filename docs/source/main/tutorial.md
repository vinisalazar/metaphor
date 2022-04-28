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

Once tests are finished, you can make use this directory to run your own analyses, since the workflow dependencies will
already be installed here. Your directory should contain the following:

```{code-block} console
(metaphor)$ ls -a

.snakemake/     data/       output/     test_data_metaphor/     
```

The `.snakemake/` and `data/` directory contain respectively the dependencies and the reference databases required to
run Metaphor. We can get rid of the `output/` and `test_data_metaphor/` directories:

```{code-block} console
(metaphor)$ rm -rf output/ test_data_metaphor/
```

Now that Metaphor is installed and we know it's working as intended, we can try running it on our own data!

## Running

To run Metaphor on your own data, you will need a CSV file containing three columns:


% ```{code-block} console
% ```

% Metaphor is packaged as a Python application that contains  
