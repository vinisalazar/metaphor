import yaml


def get_successful_completion(bool, msg):
    if not bool:
        raise Exception(
            "The workflow did not complete successfully, "
            "or there is nothing to be done (all expected outputs have been produced)."
        )
    else:
        print(msg)


def load_yaml(file_path):
    with open(file_path) as f:
        try:
            yaml_dict = yaml.safe_load(f)
        except yaml.YAMLError as exc:
            print(f"Something wrong with your config file: '{file_path}.'")
            print("Please check if it's valid YAML.")
            raise
    return yaml_dict


def write_yaml(yaml_dict, file_path):
    try:
        with open(file_path, "w") as f:
            yaml.dump(yaml_dict, f, default_flow_style=False)

        print(f"Wrote YAML to '{file_path}'.")
    except:
        print(f"Something wrong when writing YAML to file '{file_path}'.")
        print("Please check the YAML dict or the file path.")
        raise
