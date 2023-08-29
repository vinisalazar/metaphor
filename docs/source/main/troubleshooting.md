# Troubleshooting

The first thing to consider when troubleshooting your analysis is that **bioinformatics is hard**. This is not to desencourage the reader, but rather to set expectations that, as much as workflow developers strive to make their pipelines as seamless as possible, there are many factors that come into play for a successful run. Dependencies may break, the data may not be in the specific format that the program wants it, computational requirements may exceed capacity, the list goes on.

The first thing to do when your analysis does not work is to **find the cause of the problem**. It is not unusual to spend most of our time trying to figure out what actually went wrong. Once we correctly identify our problem, it becomes much easier to fix.

This page is intended to document some common problems that may occur when running Metaphor, and we hope it may even help to debug other Snakemake pipelines. Hopefully it will help user to identify some common problems and quickly fix them. If you cannot find your problem on this page, we encourage you to please [**open an issue**](https://github.com/vinisalazar/metaphor/issues/new/choose) describing it.

## Formatting input data

Snakemake calculates the workflow's execution graph before starting to run each rule. So, if your analysis does not start, it is likely that there is a problem with your **input data**.

Metaphor requires two files to run: the **settings YAML file** (default name `metaphor_settings.yaml`) and the **sample file paths tabular file** (default name `samples.csv`). When you don't specify any of these for the `metaphor execute` command, Metaphor is going to look for `metaphor_settings.yaml` to start your analysis. Your tabular file should be listed in that settings file, in the `samples` field.

It is important that your samples file has at least three columns: `sample_name` with the unique identifier of each sample, `R1` with the file path to the **forward reads** and `R2` with the file path to the **reverse reads**. You may also want to include a `group` column to perform different strategies for assembly and binning. Please see the [Assembly and Binning section](./advanced.md#assembly-and-binning) for a description of how to choose your strategy.

### Creating your samples table
Metaphor offers the `metaphor config input` command to generate the sample file paths table, but it may not work depending on how the input files are structured or named. The best way for this command to work is to have all files in a directory containing no other files, and files from the same sample have the same prefix and an `_R1.fq.gz` or `_R2.fq.gz` suffix.
If the command does not work as intended, it is preferable to build your samples file manually.

Metaphor also supports having multiple files for a single sample. These are referred to as **units** and may be added to the samples file as the `unit_name` column. For example, I may have `SAMPLE001_L001_R1.fq.gz`, `SAMPLE001_L001_R2.fq.gz` and
`SAMPLE001_L002_R1.fq.gz`, `SAMPLE001_L002_R2.fq.gz`. These come from the same sample, but different lanes on the sequencer. The sample table will look like this:

| `sample_name` | `unit_name` | `R1`                               | `R2`                               |
|---------------|-------------|------------------------------------|------------------------------------|
| `SAMPLE001`   | `L001`      | `/path/to/SAMPLE001_L001_R1.fq.gz` | `/path/to/SAMPLE001_L001_R2.fq.gz` |
| `SAMPLE001`   | `L002`      | `/path/to/SAMPLE001_L002_R1.fq.gz` | `/path/to/SAMPLE001_L002_R2.fq.gz` |

Metaphor will perform QC on each unit on its own, and then (by default), it will merge all units from each sample. This can be deactivated with the `merge_reads` parameter in the YAML settings file.

## Finding logs

Whenever we encounter a problem, the first reaction is usually to read the error message in the screen. However, we may have "lost"
the error message, due to the terminal being cleared, our program being interrupted, or any other reason. We must then turn to log files.
Thankfully for us, Snakemake has great logging features. There are two different types of log files:

1. The 'workflow-level' log: this is what is printed on the screen when we run Snakemake. It shows the rules to be executed, when they are completed, the shell commands which are being run, and if there are any errors. **These logs live in the `.snakemake/logs` folder inside the work directory** where the analysis is being run. Every time you run the `metaphor execute` command (or the `snakemake` command for that matter) a new log is created with the same content that is being printed on the screen. These files are named with a timestamp
indicating when the analysis started. You can easily access the last created file using this bash command from the analysis work directory:

```{code-block} console
(metaphor)$ less $(ls -Art .snakemake/log | tail -n 1)
```

You can replace `less` with other commands, such as `more`, `nano`, or `vi`. 

## Conda environment problems

A big part of running bioinformatics workflows is setting up the programs that will be used in the analysis. Each program may have a different installation method, several dependencies, and OS-system specific requirements. 

Thankfully, [conda](https://docs.conda.io/en/latest/) solves most of our installation problems. Snakemake has [native conda support](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#integrated-package-management), and Metaphor has [individual environment files](https://github.com/vinisalazar/metaphor/tree/main/metaphor/workflow/envs) for each program that it uses. This allows Snakemake to create individual conda environments for these YAML files, and run each rule from its determined environment.

However, it may still be the case that we run into problems, and have to debug an environment installation. This is perfectly feasible if we know the right steps.

1. **Identify the environment/rule causing the problem.**

    Like mentioned at the start of this page, we must first identify what is causing the error. This can be done by looking at the error message which crashed the analysis. These are available in the directory you are running the workflow, under the `.snakemake/log` directory. For example, take this error message:

    ```{code-block} console
    [Wed Nov 30 20:16:17 2022]
    Error in rule DAS_tool:
        jobid: 1044
        output: output/binning/DAS_tool/myGroup/myGroup_DASTool_summary.tsv, output/binning/DAS_tool/myGroup/myGroup_allBins.eval
        log: output/logs/binning/DAS_tool/myGroup.log (check log file(s) for error message)
        conda-env: /data/scratch/projects/myProject/metaphor/metaphor/config/conda/7667849902b02fd8be597a1830a293cd
        shell:

            DAS_Tool -i output/binning/concoct/myGroup/concoct_scaffolds2bin.tsv,output/binning/metabat2/myGroup/metabat2_scaffolds2bin.tsv,output/binning/vamb/myGroup/vamb_scaffolds2bin.tsv \
                -c output/mapping/myGroup/myGroup_contigs_catalogue.fna \
                -o output/binning/DAS_tool/myGroup/myGroup \
                -l concoct,metabat2,vamb \
                -p output/annotation/prodigal/myGroup/myGroup_proteins.faa \
                --score_threshold 0.25 \
                --search_engine \
                --write_bins \
                --write_bin_evals \
                --threads 64 &> output/logs/binning/DAS_tool/myGroup.log

            (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)
        cluster_jobid: 42056233

    Error executing rule DAS_tool on cluster (jobid: 1044, external: 42056233, jobscript: /data/scratch/projects/myProject/results/.snakemake/tmp.vd9zix7u/snakejob.DAS_tool.1044.sh). For error details see the cluster log and the log files of the involved rule(s).
    ```

    This error message shows that the rule causing the problem is `DAS_tool`. But if we think that our problem is with the DAS Tool installation environment (by looking at the rule log file, for example), there's another piece of information that can help us:

    ```{code-block} console
    conda-env: /data/scratch/projects/myProject/metaphor/metaphor/config/conda/7667849902b02fd8be597a1830a293cd
    ```

2. **Activate the broken environment.**

    That line shows the exact environment where DAS Tool is installed. So, if you run:

    ```{code-block} console
    (metaphor)$ conda activate /data/scratch/projects/myProject/metaphor/metaphor/config/conda/7667849902b02fd8be597a1830a293cd/  # Don't forget the / at the end.
    ```

    You can activate that environment! From there, it's easier to identify the problem, to install any further dependencies (just use `conda install` since you will have the environment active already), or do any other adjustments that may be required to run that rule correctly.

    Another way to identify your conda environment is to look into the `conda/` folder. By default, Metaphor installs the conda environments inside the `config` folder where Metaphor is installed. The path to this file can be found by running the `metaphor config show` command:

    ```{code-block} console
    (metaphor)$ metaphor config show --conda
    /data/scratch/projects/myProject/metaphor/metaphor/config/conda
    ```

    If you list the contents of that directory, you will find the different environment YAML files and corresponding directories:

    ```{code-block} console
    (metaphor)$ ls $(metaphor config show --conda)
    011169a94c97a65ec5a1d60b66873dcf       4915c26470b96bb288623853b78cdc8e       7667849902b02fd8be597a1830a293cd       ba9a1a88b1b205f71adc49915258c2c6
    011169a94c97a65ec5a1d60b66873dcf.yaml  4915c26470b96bb288623853b78cdc8e.yaml  7667849902b02fd8be597a1830a293cd.yaml  ba9a1a88b1b205f71adc49915258c2c6.yaml
    ...
    ```

    We can find the DAS Tool environment by running `grep` on the YAML files:

    ```{code-block} console
    (metaphor)$ grep "das_tool" $(metaphor config show --conda)/*.yaml
    /data/scratch/projects/myProject/metaphor/metaphor/config/conda/7667849902b02fd8be597a1830a293cd.yaml:  - das_tool=1.1.4

    # Activate the environment
    (metaphor)$ conda activate /data/scratch/projects/myProject/metaphor/metaphor/config/conda/7667849902b02fd8be597a1830a293cd/

    (/data/scratch/projects/myProject/metaphor/metaphor/config/conda/7667849902b02fd8be597a1830a293cd/)$ 
    ```

3. **Make required fixes**

    Once the environment is active, you can run commands like `conda install` (or `mamba`), `conda list`, `DAS_Tool` and other programs will be on the `$PATH`, etc.

### Environment errors when starting Metaphor

Sometimes, depending on your conda configuration, you may get errors like this:

```
Downloading and installing remote packages.
CreateCondaEnvironmentException:
Could not create conda environment from /path/to/metaphor/metaphor/workflow/rules/../envs/concoct.yaml:
Command:
mamba env create --quiet --file "/path/to/metaphor/metaphor/config/conda/113f1c07c1e471e26f0161cfd3df2c4e_.yaml" --prefix "/path/to/metaphor/metaphor/config/conda/113f1c07c1e471e26f0161cfd3df2c4e_"
Output:
Could not solve for environment specs
Encountered problems while solving:
  - package concoct-1.1.0-py27h88e4a8a_0 requires python >=2.7,<2.8.0a0, but none of the providers can be installed

The environment can't be solved, aborting the operation
```

This is likely due to strict channel priorities in conda. Please edit your conda configuration file (usually `~/.condarc`) and set channel priorities to flexible while running Metaphor.

## Adding additional flags to Snakemake

Snakemake has [many different command line options.](https://snakemake.readthedocs.io/en/stable/executing/cli.html) It wouldn't be practical to add every single option to Metaphor, so, it's possible to pass any additional settings to Snakemake using the `-e` or `--extras` flag in the `metaphor execute` or `metaphor test` commands. For example, if I would like to pass the Snakemake `--rerun-incomplete` setting, I should run:

```{code-block} console
(metaphor)$ metaphor execute -y -c 32 -e " --rerun-incomplete"
```

**It is necessary to add a space** before the `--` characters in the string passed to the `-e` flag, as [argparse](https://docs.python.org/3/library/argparse.html) does not accept arguments starting with dash. This is a [documented issue](https://bugs.python.org/issue9334) with the module.
**Edit:** this has been fixed as of Metaphor v1.7.8.

## Unlocking directories

If you get the following error when trying to resume your analysis:

```{code-block} console
Config file /data/scratch/projects/myProject/metaphor/metaphor/config/default-config.yaml is extended by additional config specified via the command line.
Building DAG of jobs...
Error: Directory cannot be locked. Please make sure that no other Snakemake process is trying to create the same files in the following directory:
/data/scratch/projects/myProject/
If you are sure that no other instances of snakemake are running on this directory, the remaining lock was likely caused by a kill signal or a power loss. It can be removed with the --unlock argument.
```

Snakemake places a "lock" on work directories to ensure that two executions don't run at the simultaneously, as that may cause problems (*e.g.* two processes writing to the same file). If the analysis is interrupted suddenly, that lock may remain in place. If we read the error message, we see that Snakemake asks the lock to be removed with the `--unlock` argument. We can pass that directly to the `metaphor execute` command:

```{code-block} console
(metaphor)$ metaphor execute -y -c 32 --unlock

Unlocking working directory.
Metaphor finished successfully.
```

## Misc problems
If you run into any other problems while using Metaphor, we highly encourage you to please [**open an issue**](https://github.com/vinisalazar/metaphor/issues/new/choose) describing it. This helps us to improve it for all users and is a great way to contribute to Metaphor as an open-source package.

<!-- 
- Running workflows is not easy
- Identifying errors
- Dealing with conda environments
- Hanging processes
- Using `snakemake --unlock`
    - argparse bug https://bugs.python.org/issue47002
 -->
