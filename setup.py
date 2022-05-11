from setuptools import setup
from metaphor import __author__, __email__, __copyright__, __license__, __version__

with open("README.md", "r") as f:
    readme = f.read()

setup(
    name="metaphor",
    version=__version__,
    author=__author__,
    author_email=__email__,
    license=__license__,
    url="https://github.com/vinisalazar/metaphor",
    description="Metaphor - Metagenomic Pipeline for Short Reads",
    long_description=readme,
    long_description_content_type="text/markdown",
    packages=["metaphor", "metaphor.cli", "metaphor.workflow", "metaphor.config"],
    data_files=[
        (
            ".",
            [
                "README.md",
            ],
        )
    ],
    include_package_data=True,
    package_data={"metaphor": ["workflow/*", "workflow/*/*", "config/*"]},
    entry_points={"console_scripts": ["metaphor = metaphor.cli.cli:main"]},
    install_requires=[  # PyPI dependencies only
        "jinja2",
        "networkx",
        "pandas",
        "pyyaml",
        "requests",
        "snakemake",
        "snakemake-wrapper-utils",
        "tqdm",
    ],
    extras_require={"docs": "myst_parser"},
    keywords="metagenomics binning assembly snakemake workflow pipeline",
    python_requires=">=3.6",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.7",
    ],
)
