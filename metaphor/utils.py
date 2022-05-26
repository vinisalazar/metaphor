import sys
import yaml


def get_successful_completion(status, msg):
    """
    Prints a message if the workflow completed successfully.

    args:
        status: int
            0 if the workflow completed successfully, else retcode.
        msg: str
            Message to be printed if bool is True.
    """
    if status:
        raise Exception(
            "The workflow did not complete successfully, "
            "or there is nothing to be done (all expected outputs have been produced)."
        )
    else:
        print(msg)


def load_yaml(file_path):
    """
    Loads a YAML file as a dictionary.

    args:
        file_path: str
            String representing file path to be loaded

    returns:
        yaml_dict: dict
    """
    with open(file_path) as f:
        try:
            yaml_dict = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(f"Something wrong with your config file: '{file_path}.'")
            print("Please check if it's valid YAML.")
            raise
    return yaml_dict


def write_yaml(yaml_dict, file_path):
    """
    Writes a dictionary as a YAML file.

    args:
        yaml_dict: dict
            Python dictionary containing values to be dumped.
        file_path: str
            Path of the output file.
    """
    try:
        with open(file_path, "w") as f:
            yaml.dump(yaml_dict, f, default_flow_style=False, sort_keys=False)

        print(f"Wrote YAML to '{file_path}'.")
    except:
        print(f"Something wrong when writing YAML to file '{file_path}'.")
        print("Please check the YAML dict or the file path.")
        raise


def confirm_message(cores, mem_mb):
    """
    Prints a confirmation message requiring user input to continue.

    Used in the metaphor executable commands to confirm number of cores and RAM.

    args:
        cores: int
            Number of cores to be used in the workflow.
        mem_mb: int
            Number of MB RAM to be used in the workflow.
    """
    print(
        "You can suppress this confirmation message by running the Metaphor command with the `-y` flag.\n"
    )
    yn = input(
        f"Snakemake will start with {cores} cores and {mem_mb} MB RAM. Ok to continue? [y/N]\n"
    )
    if yn.lower() != "y":
        print("Metaphor execution cancelled.")
        sys.exit()
