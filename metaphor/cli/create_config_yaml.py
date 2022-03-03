#!/usr/bin/env python

__doc__ = """
Generate Metaphor settings YAML file.
"""
import argparse


from metaphor import ascii_art
from metaphor.config import default_config
from metaphor.utils import load_yaml, write_yaml

config = load_yaml(default_config)


def setting_prompt(message, key, subkey, transform=str, default=None):
    # Get default from config
    if default is None:
        if subkey:
            default = config[key][subkey]
        else:
            default = config[key]

    # Get suggested answer
    def get_suggestion_type(transform):
        suggestion_dict = {bool: "y/n", int: "integer", float: "float", str: "string"}
        return suggestion_dict[transform]

    suggestion = get_suggestion_type(transform)

    # Present prompt
    answer = input(
        "\n".join(
            (
                message,
                f"Expected answer type: {suggestion}",
                f"Default value: {default}",
            )
        )
    ).lower()

    # Validate answer
    suggestion = suggestion.lower().split("/")
    if not answer:
        print(f"Setting default option: {default}\n")
        answer = default
    elif answer not in suggestion:
        breakpoint()
        print(f"Answer not valid. Please answer one of: {suggestion}.")
        answer = ""

    # Transform answer type
    def transform_func(answer):
        if transform == bool:
            return True if answer in (True, "y") else False
        else:
            try:
                return transform(answer)
            except ValueError:
                print(
                    f"Could not transform {answer} to {transform}. Please check your input value."
                )
                print(f"Setting to default value: {default}\n")
                return default

    answer = transform_func(answer)

    # Set new answer
    if subkey:
        config[key][subkey] = answer
    else:
        config[key] = answer

    if answer != default:
        if subkey:
            print(f"Set config value config['{key}']['{subkey}'] to {answer}\n")
        else:
            print(f"Set config value config['{key}'] to {answer}\n")
        return answer


def get_general_settings():
    print("General settings\n")
    setting_prompt(
        "How many MB RAM per thread would you like to use?",
        "resources",
        "mb_per_thread",
        int,
    )


def get_assembly_settings():
    # Assembly
    print("Assembly settings\n")
    setting_prompt(
        "Would you like to perform co- (pooled) assembly instead of individual assembly?",
        "coassembly",
        None,
        bool,
    )
    setting_prompt(
        "Would you like to perform Metaquast evaluation of reads?",
        "metaquast",
        "activate",
        bool,
    )

    if config["metaquast"]["activate"]:
        if config["coassembly"]:
            setting_prompt(
                "Please specify the file to be used as reference for Metaquast",
                "metaquast",
                "coassembly_reference",
            )
        else:
            print(
                "Don't forget to add a 'metaquast_reference' column to the input table to specify references for each sample."
            )


def get_annotation_settings():
    # Annotation
    print("Annotation settings\n")

    custom_diamond_db = setting_prompt(
        "Would you like to specify a custom Diamond database? Press enter to use the default (NCBI COG based).",
        "diamond",
        "db",
    )

    if custom_diamond_db:
        setting_prompt(
            "Would you like to specify the path of a .fa file to generate said custom Diamond database? "
            "Press enter to use the default (NCBI COG based).",
            "diamond",
            "db_source",
        )
        setting_prompt(
            "Would you like to specify a taxonmap to your custom Diamond database? "
            "Press enter to use the default (NCBI COG based).",
            "lineage_parser",
            "taxonmap",
        )

    setting_prompt(
        "Would you like to perform functional annotation with COG? "
        "Only possible with the default COG database.",
        "cog_functional_parser",
        "activate",
        bool,
    )

    setting_prompt(
        "Would you like to plot the functional annotation results? "
        "Only possible with the default COG database.",
        "plot_cog_functional",
        "activate",
        bool,
    )

    setting_prompt(
        "Would you like to plot the taxonomy annotation results? "
        "Requires a Diamond database with taxonomy information.",
        "plot_taxonomies",
        "activate",
        bool,
    )


def get_binning_settings():
    # Binning
    print("Binning settings\n")

    vamb = setting_prompt(
        "Would you like to perform binning with vamb?", "vamb", "activate", bool
    )
    if vamb:
        setting_prompt(
            "Please set the minimum contig length for vamb bins.",
            "vamb",
            "minfasta",
            int,
        )

    setting_prompt(
        "Would you like to perform binning with metabat2?", "metabat2", "activate", bool
    )

    setting_prompt(
        "Would you like to perform binning with CONCOCT?", "concoct", "activate", bool
    )

    setting_prompt(
        "Would you like to refine bins with DAS_Tool?", "das_tool", "activate", bool
    )


def main(args):
    print(ascii_art)
    print("Please define your desired Metaphor settings.")
    print("Press ENTER to input the default setting.")
    print(f"Settings will be written to '{args.outfile}'.\n")
    get_general_settings()
    get_assembly_settings()
    get_annotation_settings()
    get_binning_settings()
    write_yaml(config, args.outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-o", "--outfile")
    args = parser.parse_args()
    main(args)