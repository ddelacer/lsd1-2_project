#!/bin/bash
################################################################################
# Data preparation lsd1/lsd2 samples published at: https://doi.org/10.3390/cells9040955
# 
# Created by David de la Cerda
# script name: bam_coverage_cpm.sh
# 
# 
# 
# input: bam files generated from fq_to_bam.sh
# output:  bigwig files for visualization on genome browser
# required software: python 3.7 samtools 1.10
################################################################################
module load rhel7/python/3.7.0
module load rhel7/samtools/1.10
export PATH="$PATH:/home/ddelacer/.local/bin/"
source ~/.bashrc
#CPM normalization 
for sample in *bam
do
    echo $sample
    describer=$(echo ${sample} |sed 's/Aligned.sortedByCoord.out.bam//')
    echo $describer
    bamCoverage -b $sample -bs 1 -p 10 --normalizeUsing CPM -o ${describer}_cpm.bw
done
mv *cpm* coverage/
#generate mean files
bigwigCompare -b1 C1_cpm.bw -b2 C1_cpm.bw --operation sum -p 10 -of bigwig -bs 1 -o C1C2_sum.bw
bigwigCompare -b1 C1C2_sum.bw -b2 C3_cpm.bw --operation mean --scaleFactors.5:1 -p 10 -of bigwig -bs 1 -o C1C2C3_mean.bw
rm C1C2_sum.bw
bigwigCompare -b1 T7_cpm.bw -b2 T8_cpm.bw --operation sum -p 10 -of bigwig -bs 1 -o T7T8_sum.bw
bigwigCompare -b1 T7T8_sum.bw -b2 T9_cpm.bw --operation mean --scaleFactors.5:1 -p 10 -of bigwig -bs 1 -o T7T8T9_mean.bw
rm T7T8_sum.bw
bigwigCompare -b1 T10_cpm.bw -b2 T11_cpm.bw --operation sum -p 10 -of bigwig -bs 1 -o T10T11_sum.bw
bigwigCompare -b1 T10T11_sum.bw -b2 T12_cpm.bw --operation mean --scaleFactors.5:1 -p 10 -of bigwig -bs 1 -o T10T11T12_mean.bw
rm T10T11_sum.bw