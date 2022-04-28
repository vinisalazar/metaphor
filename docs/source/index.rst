.. Metaphor documentation master file, created by
   sphinx-quickstart on Tue Feb  8 08:44:38 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Metaphor's documentation!
====================================

Metaphor -- Metagenomic Pipeline for Short Reads -- is a `Snakemake <https://snakemake.readthedocs.io/>`_-based workflow for
assembly and binning of metagenomes. It also performs quality control, and basic functional and taxonomic annotation of the
assembled contigs, using the `NCBI COG <https://www.ncbi.nlm.nih.gov/research/cog/>`_ database.

Please refer to this website to learn how to use Metaphor. If you have any questions or comments, or would like to report a
problem, don't hesitate to `open an issue in the GitHub repo <https://github.com/vinisalazar/metaphor/issues/new/choose>`_.
Any and all feedback is appreciated!

**Quickstart**

.. code-block:: console

   $ conda install metaphor -c bioconda
   $ metaphor execute -i path/to/fastq/

Please read through the `Introduction <usage/introduction.md>`_ for more detailed usage instructions.

**Table of Contents**

.. toctree::
   :maxdepth: 2

   usage/introduction
   usage/output
   usage/advanced
   usage/configuration
   usage/troubleshooting
   usage/contributing
   usage/reference

..
   The following section is commented.

.. Indices and tables
.. ==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
