from pathlib import Path

data_dir = Path(__file__).parent.joinpath("data/")
test_config = Path(__file__).parent.joinpath("test-config.yaml")
default_config = Path(__file__).parent.joinpath("default-config.yaml")
example_input = Path(__file__).parent.joinpath("samples.csv")
conda_prefix = Path(__file__).parent.joinpath("conda/")
