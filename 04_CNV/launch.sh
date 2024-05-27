#!/usr/bin/env bash

cd /scratch/$USER/trabajo-final-bioinformatics/04_CNV

# variant calling for each of the bams
bams_path=/scratch/$USER/trabajo-final-bioinformatics/02_align/dna_bams
bams=($(ls $bams_path/*_recalibrated.bam))

mkdir -p logs

# launch fastqc
jID=$(sbatch --array=0-$((${#bams[@]} - 1))%5 pipeline1.sbs $bams_path)
# extract number from jID
jID=$(echo $jID | cut -d ' ' -f 4)


