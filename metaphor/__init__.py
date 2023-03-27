__author__ = "Vini Salazar"
__email__ = "17276653+vinisalazar@users.noreply.github.com"
__license__ = "MIT"
__copyright__ = "The University of Melbourne 2022"
__version__ = "1.7.3"

__doc__ = """
Metaphor top-level module. This is the entrypoint for the `cli`, `config`, and `workflow` packages.

This module contains package-level dunder variables such as licence, copyright, and version,
and also variables for the Snakemake wrapper versions to be used, and the package ASCII art.
"""

# Snakemake wrapper version
# This should match the latest released tag on: https://github.com/snakemake/snakemake-wrappers
wrapper_version = "v1.23.3"
wrapper_prefix = "https://github.com/snakemake/snakemake-wrappers/raw/"

github_url = "https://github.com/vinisalazar/metaphor"
docs_url = "https://metaphor-workflow.readthedocs.io/"

ascii_art = f"""
   __  ___      __               __            
  /  |/  /___  / /_ ___ _ ___   / /  ___   ____
 / /|_/ // -_)/ __// _ `// _ \ / _ \/ _ \ / __/
/_/  /_/ \__/ \__/ \_,_// .__//_//_/\___//_/   
                       /_/                     

Metaphor v{__version__} - Metagenomic Pipeline for Short Reads

© {__copyright__}
{__license__} licence

Documentation available at '{docs_url}'.

To give feedback and report bugs or issues, please go to '{github_url}'.
"""
