__author__ = "Vini Salazar"
__email__ = "17276653+vinisalazar@users.noreply.github.com"
__license__ = "MIT"
__copyright__ = "2022, The University of Melbourne"
__version__ = "1.1.1"

__doc__ = """
Metaphor init module. This defines the workflow version as a dunder variable, as well as other helper variables.
"""

# Snakemake wrapper version
# This should match the latest released tag on: https://github.com/snakemake/snakemake-wrappers
wrapper_version = "v1.2.0"
wrapper_prefix = "https://github.com/snakemake/snakemake-wrappers/raw/"

ascii_art = f"""
   __  ___      __               __            
  /  |/  /___  / /_ ___ _ ___   / /  ___   ____
 / /|_/ // -_)/ __// _ `// _ \ / _ \/ _ \ / __/
/_/  /_/ \__/ \__/ \_,_// .__//_//_/\___//_/   
                       /_/                     

Metaphor v{__version__} - Metagenomic Pipeline for Short Reads

© The University of Melbourne 2022
"""
