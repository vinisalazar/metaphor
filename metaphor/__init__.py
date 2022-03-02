from .workflow import snakefile
from .config import test_config, default_config

# Snakemake wrapper version
# This should match the latest released tag on: https://github.com/snakemake/snakemake-wrappers
wrapper_version = "v1.2.0"
wrapper_prefix = "https://github.com/snakemake/snakemake-wrappers/raw/"


def get_successful_completion(bool, msg):
    if not bool:
        raise Exception("The workflow did not complete successfully.")
    else:
        print(msg)


ascii_art = """
   __  ___      __               __            
  /  |/  /___  / /_ ___ _ ___   / /  ___   ____
 / /|_/ // -_)/ __// _ `// _ \ / _ \/ _ \ / __/
/_/  /_/ \__/ \__/ \_,_// .__//_//_/\___//_/   
                       /_/                     

Metaphor - Metagenomic Pipeline for Short Reads

© The University of Melbourne 2022
"""
