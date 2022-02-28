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


def main(args):

    # Initial sanity check
    config_file = args.configfile
    input_dir = args.input_dir
    join_units = args.join_units
    cores = int(args.cores)
    coassembly = args.coassembly
    confirm = args.confirm
    assert Path(config_file).exists(), f"Could not find config file: {config_file}"

    if args.profile:
        # Code taken from snakemake.__init__.py main function
        from snakemake import get_argument_parser, get_profile_file

        # reparse args while inferring config file from profile
        parser = get_argument_parser(args.profile)
        args = parser.parse_args()

        def adjust_path(f):
            if os.path.exists(f) or os.path.isabs(f):
                return f
            else:
                return get_profile_file(args.profile, f, return_default=True)

        # update file paths to be relative to the profile
        # (if they do not exist relative to CWD)
        if args.jobscript:
            args.jobscript = adjust_path(args.jobscript)
        if args.cluster:
            args.cluster = adjust_path(args.cluster)
        if args.cluster_config:
            if isinstance(args.cluster_config, list):
                args.cluster_config = [adjust_path(cfg) for cfg in args.cluster_config]
            else:
                args.cluster_config = adjust_path(args.cluster_config)
        if args.cluster_sync:
            args.cluster_sync = adjust_path(args.cluster_sync)
        for key in "cluster_status", "cluster_cancel", "cluster_sidecar":
            if getattr(args, key):
                setattr(args, key, adjust_path(getattr(args, key)))
        if args.report_stylesheet:
            args.report_stylesheet = adjust_path(args.report_stylesheet)
    else:
        for arg_ in (
            "jobscript",
            "cluster",
            "cluster_config",
            "cluster_sync",
            "cluster_status",
            "cluster_cancel",
            "cluster_sidecar",
            "report_stylesheet",
        ):
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
    with open(config_file) as f:
        try:
            config = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(f"Something wrong with your config file: '{config_file}.'")
            print("Please check if it's valid YAML.")
            raise

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
