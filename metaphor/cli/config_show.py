#!/usr/bin/env python

from re import L
from metaphor import __path__ as metaphor_path
from metaphor.config import default_config, test_config, example_input
from metaphor.workflow import snakefile

__doc__ = "Show path of Metaphor installation, Snakefile, config files, and a sample input file."


metaphor_path = metaphor_path[0]


def main(*args):
    choices = [
        "metaphor_path",
        "snakefile",
        "test_config",
        "default_config",
        "example_input",
    ]
    choices = {var: (var.replace("_", " ").capitalize(), eval(var)) for var in choices}
    print_all = not any(getattr(args[0], key, False) for key in choices.keys())

    for key, value in choices.items():
        if not print_all:
            if attr := getattr(args[0], key, False):
                print(value[1])
        else:
            print(f"{value[0]}:\t{value[1]}")


if __name__ == "__main__":
    main()
