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
# output:  gene count quantification file
# required software: featureCounts 1.6.3
################################################################################
for sample in *bam
do
    echo $sample
    describer=$(echo ${sample} |sed 's/Aligned.sortedByCoord.out.bam//')
    echo $describer
    #gene
    /deac/bio/peaseGrp/ddelacer/subread-1.6.3-source/bin/featureCounts -t gene -g gene_id -a Schizosaccharomyces_pombe.ASM294v2.43.gtf -G Schizosaccharomyces_pombe_all_chromosomes_.fa -O -M --fraction -T 8 -p -o ${describer}_gene.txt $sample
    cut -f1,7-8 ${describer}_gene.txt | sed -e '1d' - > ${describer}_gene_final.txt
done
mv *gene* gene_counts/
cd gene_counts
paste *final.txt > Spombe_gene_withO_fraction.txt