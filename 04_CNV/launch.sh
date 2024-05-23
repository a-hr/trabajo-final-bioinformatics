#!/usr/bin/env bash

REFERENCE=/scratch/heral/TRABAJO_FINAL/samples/references/GRCh38.primary_assembly.genome.fa
DICT_REFERENCE=/scratch/aberbelgara/NGS22/samples/Homo_sapiens_assembly38.dict
OUTPUT_DIR=/scratch/aberbelgara/NGS22/Output
INTERVALS=/scratch/aberbelgara/NGS22/samples/tutorial_11682/targets_C.interval_list
TUMOR=/scratch/aberbelgara/NGS22/samples/tutorial_11682/tumor.bam
GATK=/scratch/aberbelgara/NGS22/gatk-4.3.0.0/gatk


reference_fa=/scratch/heral/TRABAJO_FINAL/samples/references/GRCh38.primary_assembly.genome.fa

# variant calling for each of the bams
bams_path=/scratch/heral/TRABAJO_FINAL/02_align/dna_bams
bams=($(ls $bams_path/*.bam))

mkdir -p logs

# base recalibration
jID2=$(sbatch --array=0-$(( ${#bams[@]} - 1 )) --dependency=afterok:$jID1 base_recalibrator.sbs $bams_path ./vcf $reference_fa no)
jID2=$(echo $jID2 | cut -d ' ' -f 4)  # extract number from jID


