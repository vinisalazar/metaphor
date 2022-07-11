__author__ = "Vini Salazar"
__email__ = "17276653+vinisalazar@users.noreply.github.com"
__license__ = "MIT"
__copyright__ = "The University of Melbourne 2022"
__version__ = "1.5.2"

__doc__ = """
Metaphor top-level module. This is the entrypoint for the `cli`, `config`, and `workflow` packages.

This module contains package-level dunder variables such as licence, copyright, and version,
and also variables for the Snakemake wrapper versions to be used, and the package ASCII art.
"""

# Snakemake wrapper version
# This should match the latest released tag on: https://github.com/snakemake/snakemake-wrappers
wrapper_version = "v1.7.0"
wrapper_prefix = "https://github.com/snakemake/snakemake-wrappers/raw/"

ascii_art = f"""
   __  ___      __               __            
  /  |/  /___  / /_ ___ _ ___   / /  ___   ____
 / /|_/ // -_)/ __// _ `// _ \ / _ \/ _ \ / __/
/_/  /_/ \__/ \__/ \_,_// .__//_//_/\___//_/   
                       /_/                     

Metaphor v{__version__} - Metagenomic Pipeline for Short Reads

© {__copyright__}
{__license__} licence
"""
