#!/bin/bash
# Usage: interleave_fastq.sh f.fastq r.fastq > interleaved.fastq
# 
# Interleaves the reads of two FASTQ files specified on the
# command line and outputs a single FASTQ file of STDOUT.
# 
# Can interleave 100 million paired reads (200 million total
# reads; a 2 x 22Gbyte files), in memory (/dev/shm), in 6m54s (414s)
# 
# Latest code: https://gist.github.com/4544979
# Also see my deinterleaving script: https://gist.github.com/3521724

paste $1 $2 | paste - - - - | awk -v OFS="\n" -v FS="\t" '{print($1,$3,$5,$7,$2,$4,$6,$8)}'
