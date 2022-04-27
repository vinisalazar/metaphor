#!/usr/bin/env python

from metaphor import __path__ as metaphor_path
from metaphor.config import default_config, test_config, example_input
from metaphor.workflow import snakefile

__doc__ = "Show path of Metaphor installation, Snakefile, config files, and a sample input file."


metaphor_path = metaphor_path[0]


def main(*args):
    paths = [
        "metaphor_path",
        "snakefile",
        "test_config",
        "default_config",
        "example_input",
    ]
    for path in paths:
        fmt_name = path.replace("_", " ").capitalize()
        value = eval(path)
        print(f"{fmt_name}:\t{value}")


if __name__ == "__main__":
    main()
