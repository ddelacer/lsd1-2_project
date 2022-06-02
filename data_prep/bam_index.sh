#!/bin/bash
################################################################################
# Data preparation lsd1/lsd2 samples published at: https://doi.org/10.3390/cells9040955
# 
# Created by David de la Cerda
# script name: bam_index.sh
# 
# 
# 
# input: bam files generated from fq_to_bam.sh
# output:  bam files
# required software: samtools 1.10
################################################################################
module load rhel7/samtools/1.10
for sample in *bam
do
    echo $sample
    samtools index $sample
done