#!/usr/bin/env python

__doc__ = """
    Tests Metaphor.

    1. Checks for test data.
        1a. Downloads if doesn't match the checksum.
    2. Create input file.
    3. Run workflow with default config.
    4. Generate report.
    5. Delete output.
    """


import os
import sys
from argparse import Namespace
from pathlib import Path

import yaml
from snakemake import snakemake

from metaphor import get_successful_completion, snakefile, default_config, ascii_art
from .create_input_table import main as create_input_table


def load_yaml(file_path):
    with open(file_path) as f:
        try:
            yaml_dict = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(f"Something wrong with your config file: '{file_path}.'")
            print("Please check if it's valid YAML.")
            raise
    return yaml_dict


def main(args):

    # Initial sanity check
    config_file = args.configfile
    input_dir = args.input_dir
    join_units = args.join_units
    cores = int(args.cores)
    coassembly = args.coassembly
    confirm = args.confirm
    assert Path(config_file).exists(), f"Could not find config file: {config_file}"

    profile_args = (
        "jobscript",
        "cluster",
        "cluster_config",
        "cluster_sync",
        "cluster_status",
        "cluster_cancel",
        "cluster_sidecar",
        "report_stylesheet",
    )

    if args.profile:
        # Code taken from snakemake.__init__.py main function
        from snakemake import get_profile_file

        # reparse args while inferring config file from profile
        profile_file = get_profile_file(args.profile, "config.yaml")
        profile = load_yaml(profile_file)

        def adjust_path(f):
            if os.path.exists(f) or os.path.isabs(f):
                return f
            else:
                return get_profile_file(args.profile, f, return_default=True)

        # update file paths to be relative to the profile
        # (if they do not exist relative to CWD)
        for key, value in profile.items():
            if key in profile_args:
                if isinstance(value, list):
                    setattr(args, key, [adjust_path(cfg) for cfg in value])
                else:
                    setattr(args, key, adjust_path(value))
            else:
                setattr(args, key, None)
    else:
        for arg_ in profile_args:
            setattr(args, arg_, None)

    if config_file == default_config:
        if not confirm:
            yn = input(
                f"You haven't specified a custom config file. Metaphor will run with the default config: {config_file}.\n"
                "This does not include all of Metaphor's features. We recommend you create your own config file.\n"
                "Ok to continue? [y/N]"
            )
        if yn.lower() != "y":
            print("Metaphor execution cancelled.")
            sys.exit()

    # Load config file to check for input_dir
    config = load_yaml(config_file)

    samples_file = config["samples"]
    mem_mb = config["resources"]["mb_per_thread"]
    if coassembly is None:
        coassembly = config["coassembly"]
    if not input_dir:
        assert Path(
            samples_file
        ).exists(), f"Samples file {samples_file} doesn't exist. Please provide an input directory with FASTQ or a valid samples file."
    else:
        create_input_table_args = Namespace()
        create_input_table_args.input_dir = input_dir
        create_input_table_args.output_file = samples_file = str(
            Path(input_dir).joinpath("samples.csv")
        )
        create_input_table_args.join_units = join_units
        print("\nCreating input table for test files.\n")
        create_input_table(create_input_table_args)

    print("Starting Snakemake.")
    print(
        "This may require the installation of conda environments which should take a while.\n"
    )
    if not confirm:
        yn = input(
            f"Snakemake will start with {cores} cores and {mem_mb} MB RAM PER THREAD. Ok to continue? [y/N]"
        )
        if yn.lower() != "y":
            print("Metaphor execution cancelled.")
            sys.exit()
    print(ascii_art)
    sm_exit = snakemake(
        snakefile=snakefile,
        configfiles=[
            config_file,
        ],
        config={
            "samples": samples_file,
            "coassembly": coassembly,
        },
        cores=cores,
        use_conda=True,
        printshellcmds=True,
        # profile settings
        report_stylesheet=args.report_stylesheet,
        cluster=args.cluster,
        cluster_config=args.cluster_config,
        cluster_sync=args.cluster_sync,
    )
    get_successful_completion(sm_exit, "Metaphor finished successfully.")
