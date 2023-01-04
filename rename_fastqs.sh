#!/bin/bash

# Run from within sample directory with fastqs
for f in *.R1.fastq.gz; do mv -v "$f" "${f/.R1/_R1}"; done;
for f in *.R2.fastq.gz; do mv -v "$f" "${f/.R2/_R2}"; done;