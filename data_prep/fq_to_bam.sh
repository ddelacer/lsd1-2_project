#!/bin/bash
################################################################################
# Data preparation lsd1/lsd2 samples published at: https://doi.org/10.3390/cells9040955
# 
# Created by David de la Cerda
# script name: fq_to_bam.sh
# 
# 
# 
# input: fq files uploaded to GEO GSE148191
# output: mapped bam files
# required software: star 2.7.2d
################################################################################
module load rhel7/compilers/intel-2018-lp64
module load rhel7/star/2.7.2d
STAR --runMode genomeGenerate --genomeDir /deac/bio/peaseGrp/ddelacer/pombe/ --genomeFastaFiles Schizosaccharomyces_pombe_all_chromosomes_.fa --sjdbGTFfile Schizosaccharomyces_pombe.ASM294v2.43.gtf --sjdbOverhang 49 --genomeSAindexNbases 10
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAAQRAAPEI-30_1.fq --outFileNamePrefix C1 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAARRAAPEI-32_1.fq --outFileNamePrefix C2 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAAARAAPEI-14_1.fq --outFileNamePrefix C3 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAAHRAAPEI-21_1.fq --outFileNamePrefix T10 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAAIRAAPEI-22_1.fq --outFileNamePrefix T11 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAAJRAAPEI-23_1.fq --outFileNamePrefix T12 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAAKRAAPEI-24_1.fq --outFileNamePrefix T13 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAALRAAPEI-25_1.fq --outFileNamePrefix T14 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAAMRAAPEI-26_1.fq --outFileNamePrefix T15 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAANRAAPEI-27_1.fq --outFileNamePrefix T16 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAAORAAPEI-28_1.fq --outFileNamePrefix T17 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAAPRAAPEI-29_1.fq --outFileNamePrefix T18 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180515_I91_CL100071967_L2_HK500SCHorhRAABRAAPEI-15_1.fq,180506_I91_CL100071673_L2_HK500SCHorhRAABRAAPEI-15_1.fq --outFileNamePrefix T4 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAACRAAPEI-16_1.fq --outFileNamePrefix T5 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAADRAAPEI-17_1.fq --outFileNamePrefix T6 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAAERAAPEI-18_1.fq --outFileNamePrefix T7 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAAFRAAPEI-19_1.fq --outFileNamePrefix T8 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded
STAR --genomeDir /pombe_genome/pombe-genome-starIndex/ --readFilesIn 180506_I91_CL100071673_L2_HK500SCHorhRAAGRAAPEI-20_1.fq --outFileNamePrefix T9 --outSAMtype BAM Unsorted SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded