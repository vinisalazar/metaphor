#!/usr/bin/env python

__doc__ = """
Metaphor CLI - wraps commands for easier execution.

Commands: 
    metaphor
        test
        config
        create
        execute
"""


import argparse

from .metaphor_test import main as metaphor_test


def main():
    parser = argparse.ArgumentParser(prog="metaphor", description=__doc__)
    subparsers = parser.add_subparsers(help="Command to be executed.")

    test = subparsers.add_parser(
        "test",
        help="Test Metaphor with example data.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    test.add_argument("-d", "--directory", help="Directory to run tests.")
    test.add_argument("-p", "--cores", help="Number of processors to use in tests.")
    test.add_argument("-m", "--mem_mb", help="Amount of RAM to use in tests.")
    test.add_argument(
        "-c",
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
    test.set_defaults(
        func=metaphor_test, directory="metaphor_test", cores=4, mem_mb=4096
    )

    config = subparsers.add_parser(
        "config", help="Shows the location of Metaphor installation and config values."
    )

    create = subparsers.add_parser(
        "create", help="Create input sample table or input config."
    )

    execute = subparsers.add_parser("execute", help="Runs the workflow on your data!")
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
