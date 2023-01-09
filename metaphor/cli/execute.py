#!/usr/bin/env python

__doc__ = """
    Executes Metaphor on real data.
    """


import sys
from argparse import Namespace
from pathlib import Path
from shutil import copyfile
from datetime import datetime

from metaphor import wrapper_prefix
from metaphor.workflow import snakefile
from metaphor.config import default_config, conda_prefix
from metaphor.utils import (
    confirm_message,
    get_successful_completion,
    load_yaml,
    run_cmd,
)

from .create_input_table import main as create_input_table


def main(args):

    # Initial sanity check
    config_file = args.configfile
    input_dir = args.input_dir
    join_units = args.join_units
    cores = args.cores
    max_mb = args.max_mb
    confirm = args.confirm
    conda_prefix_arg = args.conda_prefix
    extras = args.extras
    profile = args.profile
    skip_report = args.skip_report

    if not Path(config_file).exists():
        if not confirm:
            yn = input(
                f"Config file '{config_file}' does not exist yet. Metaphor will create one based on the default config contained on: '{default_config}'.\n"
                "This does not include all of Metaphor's features. We recommend you create your own config file with the 'metaphor config settings' command.\n"
                "Ok to continue? [y/N]\n"
            )
        if yn.lower() == "y":
            print(f"Copying default config {config_file} to current directory.")
            copyfile(default_config, config_file)
        else:
            print("Metaphor execution cancelled.")
            sys.exit()

    # Load config file to check for input_dir
    config = load_yaml(config_file)

    samples = config["samples"]
    if max_mb is None:
        max_mb = config["max_mb"]
    if not input_dir:
        assert Path(
            samples
        ).exists(), f"Samples file {samples} doesn't exist. Please provide an input directory with FASTQ or a valid samples file."
    else:
        create_input_table_args = Namespace()
        create_input_table_args.input_dir = input_dir
        create_input_table_args.output_file = samples = "samples.csv"
        create_input_table_args.join_units = join_units
        print("\nCreating input table for test files.\n")
        create_input_table(create_input_table_args)

    print(f"Your config file is '{config_file}'.\n")

    print("Starting Snakemake.")
    print(
        "This may require the installation of conda environments which should take a while.\n"
    )
    if not confirm:
        confirm_message(cores, max_mb)

    cmd = f"""
    snakemake --snakefile {snakefile}           \
              --configfile {config_file}        \
              --cores {cores}                   \
              -p -r                             \
              --use-conda                       \
              --wrapper-prefix {wrapper_prefix}
    """

    if conda_prefix_arg:
        cmd += f"  --conda-prefix {conda_prefix_arg}"
    else:
        cmd += f"  --conda-prefix {conda_prefix}"

    for arg in ("samples", "max_mb"):
        if value := eval(arg):
            cmd += f" --config {arg}={value} "

    if profile:
        cmd += f" --profile {profile} "

    cmd += f" {extras} "

    retcode = run_cmd(cmd)
    msg = "Metaphor finished successfully."
    if config["prokka"]["activate"]:
        msg += " Please rerun it if you would like to annotate the genome bins."

    get_successful_completion(retcode, msg)

    # Don't run report if running unlock, lint or cleanup metadata option
    options = ["unlock", "lint", "cleanup-metadata", "dry-run", "dryrun"]
    skip_report = True if any(f"--{o}" in extras for o in options) else skip_report

    if not skip_report:
        timestamp = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        fileout = f"metaphor_report_{timestamp}.html"
        cmd += f" --report {fileout}"
        format_cmd = cmd.split()
        format_cmd = list(zip(format_cmd[1::2], format_cmd[2::2]))
        format_cmd = "snakemake\t\\\n\t" + "\t\t\\\n\t".join([" ".join(t) for t in format_cmd])
        print("Your command is:")
        print(format_cmd)
        retcode = run_cmd(cmd)
        get_successful_completion(
            retcode, f"Metaphor finished successfully and generated the report above."
        )
