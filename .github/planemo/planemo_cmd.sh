#!/bin/bash
planemo tool_init --id 'metaphor' \
       --name 'Metaphor - Assembly and binning of metagenomes'  \
       --description 'Metaphor is a workflow for assembly and binning of metagenomes. Use it to generate MAGs and obtain contig-level functional taxonomic annotation.' \
       --example_command 'metaphor execute -f metaphor_settings.yaml -c 12' \
       --input metaphor_settings.yaml \
       --output output \
       --doi '10.1101/2023.02.09.527784' \
       --cite_url 'https://www.biorxiv.org/content/10.1101/2023.02.09.527784' \
       --help_from_command 'metaphor execute --help' \
       --requirement 'metaphor@1.7.7' \
       --force
