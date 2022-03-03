#!/usr/bin/env python


import argparse

from metaphor import __version__
from .test import main as metaphor_test
from .execute import main as metaphor_execute
from .create_input_table import main as create_input_table_main
from .create_input_table import __doc__ as create_input_table_doc
from .create_config_yaml import main as create_config_yaml_main
from .create_config_yaml import __doc__ as create_config_yaml_doc

__doc__ = f"""
Metaphor v{__version__}  CLI - wraps commands for easier execution.

Commands: 
    metaphor
        test
        config
        execute
"""


def main():
    parser = argparse.ArgumentParser(prog="metaphor", description=__doc__)
    subparsers = parser.add_subparsers(help="Command to be executed.")

    ###############################################################
    # Execute
    # Execute command subparser
    ###############################################################

    execute = subparsers.add_parser(
        "execute",
        help="Execute Metaphor on real data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    execute.add_argument(
        "-i", "--input-dir", help="Input directory containing FASTQ files."
    )
    execute.add_argument(
        "-c", "--configfile", help="Configuration file to run the workflow."
    )
    execute.add_argument(
        "-j",
        "--join-units",
        action="store_true",
        help="Whether to join units (S001, S002) with the same preffix as the same file.",
    )
    execute.add_argument("-p", "--cores", help="Number of processors to use in tests.")
    execute.add_argument("-l", "--profile", help="Profile to be used to run Metaphor.")
    execute.add_argument(
        "-co",
        "--coassembly",
        action="store_true",
        help="Whether to run tests in coassembly mode, "
        "i.e. all samples are pooled together and assembled.",
    )
    execute.add_argument(
        "-y",
        "--confirm",
        action="store_true",
        help="Don't ask for confirmation when running tests.",
    )
    execute.add_argument(
        "--until",
        help="Only run workflow until the listed files are generated (see Snakemake docs). "
        "Files must be space separated with the complete output path.",
        nargs="+",
    )

    execute.set_defaults(
        func=metaphor_execute,
        input_dir=None,
        configfile="metaphor_settings.yaml",
        join_units=False,
        cores=8,
        profile=None,
        coassembly=None,
    )

    ###############################################################
    # TEST
    # Test command subparser
    ###############################################################

    test = subparsers.add_parser(
        "test",
        help="Test Metaphor with example data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    test.add_argument("-d", "--directory", help="Directory to run tests.")
    test.add_argument("-p", "--cores", help="Number of processors to use in tests.")
    test.add_argument("-m", "--mem_mb", help="Amount of RAM to use in tests.")
    test.add_argument(
        "-co",
        "--coassembly",
        action="store_true",
        help="Whether to run tests in coassembly mode, "
        "i.e. all samples are pooled together and assembled.",
    )
    test.add_argument(
        "-y",
        "--confirm",
        action="store_true",
        help="Don't ask for confirmation when running tests.",
    )
    test.add_argument(
        "--remove-conda",
        help="If this option is selected, conda environments will be created in the test "
        "directory instead of the current directory (and therefore are deleted when tests finish).",
        action="store_true",
    )
    test.add_argument(
        "-dry",
        "--dry-run",
        action="store_true",
        help="Whether to run tests as a dry-run only (used for CI).",
    )
    test.set_defaults(
        func=metaphor_test,
        directory="metaphor_test",
        cores=2,
        mem_mb=4096,
        dry_run=False,
    )

    ###############################################################
    # CONFIG
    # Config command subparser
    ###############################################################

    config = subparsers.add_parser(
        "config", help="Sets up Metaphor input data and settings file."
    )
    config_subparsers = config.add_subparsers(help="Config command to be executed.")

    # Input table
    ###############################################################
    create_input_table = config_subparsers.add_parser(
        "input",
        help=create_input_table_doc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    create_input_table.add_argument(
        "-i",
        "--input-dir",
        help="Input directory containing FASTQ files.",
        required=True,
    )
    create_input_table.add_argument(
        "-j",
        "--join-units",
        help="If this option is on, files with the same preffix but with "
        "S001, S002, S00N distinctions in the filenames will be treated as different units of the same sample, "
        "i.e. they will be joined into a single file.",
        action="store_true",
    )
    create_input_table.add_argument("-o", "--output-file", help="Path to output file.")
    create_input_table.set_defaults(func=create_input_table_main)

    # Config YAML
    ###############################################################
    create_config_yaml = config_subparsers.add_parser(
        "settings",
        help=create_config_yaml_doc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    create_config_yaml.add_argument(
        "-o", "--outfile", help="Output path of config YAML file."
    )
    create_config_yaml.set_defaults(
        func=create_config_yaml_main, outfile="metaphor_settings.yaml"
    )

    ###############################################################
    # PARSE ALL ARGS
    ###############################################################
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
