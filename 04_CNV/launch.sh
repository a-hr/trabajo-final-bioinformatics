#!/usr/bin/env bash



REFERENCE=/scratch/$USER/trabajo-final-bioinformatics/samples/references/GRCh38.primary_assembly.genome.fa
DICT_REFERENCE=/scratch/$USER/trabajo-final-bioinformatics/samples/references/Homo_sapiens_assembly38.dict
OUTPUT_DIR=/scratch/$USER/trabajo-final-bioinformatics/04_CNV/Output
INTERVALS=/scratch/$USER/trabajo-final-bioinformatics/samples/tutorial_11682/targets_C.interval_list
TUMOR=/scratch/$USER/trabajo-final-bioinformatics/samples/tutorial_11682/tumor.bam #este va a ser el bam que en el sbs lo meta como TUMOR
GATK=/scratch/$USER/trabajo-final-bioinformatics/bin/gatk-4.3.0.0/gatk


reference_fa=/scratch/$USER/trabajo-final-bioinformatics/samples/references/GRCh38.primary_assembly.genome.fa

# variant calling for each of the bams
bams_path=/scratch/$USER/trabajo-final-bioinformatics/02_align/dna_bams
bams=($(ls $bams_path/*_recalibrated.bam))

mkdir -p logs

# base recalibration
jID2=$(sbatch --array=0-$(( ${#bams[@]} - 1 )) --dependency=afterok:$jID1 base_recalibrator.sbs $bams_path ./vcf $reference_fa no)
jID2=$(echo $jID2 | cut -d ' ' -f 4)  # extract number from jID


