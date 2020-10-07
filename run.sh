#!/usr/bin/env bash

export NXF_SINGULARITY_CACHEDIR=/home/bs/Git/nf-tnseq/singularity
NEXTFLOW=/home/bs/Downloads/nextflow
$NEXTFLOW run tn_seq.nf \
  --with-dag flowchart.png \
  --with-timeline timeline.html \
  --with-report nextflow_report.html
