#!/usr/bin/env python

__doc__ = """
    Executes Metaphor on real data.
    """


import os
import sys
from argparse import Namespace
from pathlib import Path
from shutil import copyfile

from snakemake import snakemake

from metaphor import wrapper_prefix
from metaphor.workflow import snakefile
from metaphor.config import default_config
from metaphor.utils import get_successful_completion, load_yaml

from .create_input_table import main as create_input_table


def main(args):

    # Initial sanity check
    config_file = args.configfile
    input_dir = args.input_dir
    join_units = args.join_units
    cores = int(args.cores)
    coassembly = args.coassembly
    confirm = args.confirm
    until = args.until
    if until and not isinstance(until, list):
        until = until.split()

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

    if not Path(config_file).exists():
        if not confirm:
            yn = input(
                f"Config file '{config_file}' does not exist yet. Metaphor will create one based on the default config contained on: '{default_config}'.\n"
                "This does not include all of Metaphor's features. We recommend you create your own config file with the 'metaphor config' command.\n"
                "Ok to continue? [y/N]"
            )
        if yn.lower() == "y":
            print(f"Copying default config {config_file} to current directory.")
            copyfile(default_config, config_file)
        else:
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
        wrapper_prefix=wrapper_prefix,
        until=[] if not until else until,
        # profile settings
        report_stylesheet=vars(args).get("report_stylesheet", None),
        cluster=vars(args).get("cluster", None),
        cluster_config=vars(args).get("cluster_config", None),
        cluster_sync=vars(args).get("cluster_sync", None),
    )
    get_successful_completion(sm_exit, "Metaphor finished successfully.")
