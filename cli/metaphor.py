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

from metaphor_test import metaphor_test


def main():
    parser = argparse.ArgumentParser(prog="metaphor", description=__doc__)
    subparsers = parser.add_subparsers(help="Command to be executed.")

    test = subparsers.add_parser("test", help="Test Metaphor with example data.")
    test.add_argument("-d", "--directory", help="Directory to run tests.")
    test.add_argument(
        "--preserve-conda",
        help="If this option is selected, conda environments will be created in the current "
             "directory instead of the test directory (which is deleted when tests finish).",
        action="store_true",
    )
    test.set_defaults(func=metaphor_test, directory="metaphor_test")

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
