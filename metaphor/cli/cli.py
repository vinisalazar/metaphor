#!/usr/bin/env python


import argparse

from metaphor import __version__
from .test import main as metaphor_test
from .execute import main as metaphor_execute
from .create_input_table import main as create_input_table
from .create_input_table import __doc__ as create_input_table_doc
from .create_config_yaml import main as create_config_yaml
from .create_config_yaml import __doc__ as create_config_yaml_doc
from .config_show import main as config_show
from .config_show import __doc__ as config_show_doc

__doc__ = f"""
Metaphor v{__version__}  CLI - wraps commands for easier execution.
"""


def main():
    """
    Main function of the Metaphor CLI application.

    This function defines the parsers and subparsers for each of the
    subcommands in the metaphor executable, and calls functions
    from the other scripts defined in this subpackage.
    """
    parser = argparse.ArgumentParser(prog="metaphor", description=__doc__)
    subparsers = parser.add_subparsers(help="Command to be executed.", dest="subparser")

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
        "-f", "--configfile", help="Configuration file to run the workflow."
    )
    execute.add_argument(
        "-j",
        "--join-units",
        action="store_true",
        help="Whether to join units (S001, S002) with the same preffix as the same file.",
    )
    execute.add_argument(
        "-c", "--cores", help="Number of processors to be used.", type=int
    )
    execute.add_argument("-l", "--profile", help="Profile to be used to run Metaphor.")
    execute.add_argument(
        "-m",
        "--mem_mb",
        help="Amount of MB RAM to be used PER CORE. "
        "i.e. if you set 1024 and 2 cores, it will use up to 2048 MB of RAM.",
        type=int,
    )
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
        cores=4,
        mem_mb=None,
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
    test.add_argument(
        "-c", "--cores", help="Number of processors to use in tests.", type=int
    )
    test.add_argument("-m", "--mem_mb", help="Amount of RAM to use in tests.", type=int)
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
        directory="test_data_metaphor",
        cores=3,
        mem_mb=1024,
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
    create_input_table_parser = config_subparsers.add_parser(
        "input",
        help=create_input_table_doc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    create_input_table_parser.add_argument(
        "-i",
        "--input-dir",
        help="Input directory containing FASTQ files.",
        required=True,
    )
    create_input_table_parser.add_argument(
        "-j",
        "--join-units",
        help="If this option is on, files with the same preffix but with "
        "S001, S002, S00N distinctions in the filenames will be treated as different units of the same sample, "
        "i.e. they will be joined into a single file.",
        action="store_true",
    )
    create_input_table_parser.add_argument(
        "-o", "--output-file", help="Path to output file."
    )
    create_input_table_parser.set_defaults(func=create_input_table)

    # Config YAML
    ###############################################################
    create_config_yaml_parser = config_subparsers.add_parser(
        "settings",
        help=create_config_yaml_doc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    create_config_yaml_parser.add_argument(
        "-o", "--outfile", help="Output path of config YAML file."
    )
    create_config_yaml_parser.set_defaults(
        func=create_config_yaml, outfile="metaphor_settings.yaml"
    )

    # Show filepaths
    ###############################################################
    show_config_paths_parser = config_subparsers.add_parser(
        "show",
        help=config_show_doc,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    show_config_paths_exc_group = (
        show_config_paths_parser.add_mutually_exclusive_group()
    )
    show_config_paths_exc_group.add_argument(
        "--metaphor-path",
        help="Show path of the Metaphor package installation.",
        action="store_true",
    )
    show_config_paths_exc_group.add_argument(
        "--snakefile", help="Show path of the Metaphor Snakefile.", action="store_true"
    )
    show_config_paths_exc_group.add_argument(
        "--test_config",
        help="Show path of the Metaphor test configuration file.",
        action="store_true",
    )
    show_config_paths_exc_group.add_argument(
        "--default_config",
        help="Show path of the Metaphor default configuration file.",
        action="store_true",
    )
    show_config_paths_exc_group.set_defaults(func=config_show)

    ###############################################################
    # PARSE ALL ARGS
    ###############################################################
    args = parser.parse_args()
    try:
        args.func(args)
    except AttributeError:
        eval(args.subparser).print_help()  # get subparser from command and print help


if __name__ == "__main__":
    main()
