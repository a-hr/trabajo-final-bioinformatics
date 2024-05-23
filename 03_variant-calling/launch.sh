#!/usr/bin/env bash

reference_fa=/scratch/$USER/trabajo-final-bioinformatics/samples/references/GRCh38.primary_assembly.genome.fa

# variant calling for each of the bams
bams_path=/scratch/$USER/trabajo-final-bioinformatics/02_align/dna_bams
bams=($(ls $bams_path/*.bam))

mkdir -p logs

# index creation
if [ ! -f $(dirname $reference_fa)/$(basename $reference_fa .fa).dict ]; then
    jID0=$(sbatch index.sbs $reference_fa)
    jID0=$(echo $jID0 | cut -d ' ' -f 4)  # extract number from jID

    # variant calling
    jID1=$(sbatch --array=0-$(( ${#bams[@]} - 1 )) --dependency=afterok:$jID0 haplotype_caller.sbs $bams_path $reference_fa no)
    jID1=$(echo $jID1 | cut -d ' ' -f 4)  # extract number from jID

else
    # variant calling
    jID1=$(sbatch --array=0-$(( ${#bams[@]} - 1 )) haplotype_caller.sbs $bams_path $reference_fa no)
    jID1=$(echo $jID1 | cut -d ' ' -f 4)  # extract number from jID

fi

# base recalibration
jID2=$(sbatch --array=0-$(( ${#bams[@]} - 1 )) --dependency=afterok:$jID1 base_recalibrator.sbs $bams_path ./vcf $reference_fa no)
jID2=$(echo $jID2 | cut -d ' ' -f 4)  # extract number from jID

# apply recalibration
jID3=$(sbatch --array=0-$(( ${#bams[@]} - 1 )) --dependency=afterok:$jID2 apply_recalibrator.sbs $bams_path ./stats_reports)
jID3=$(echo $jID3 | cut -d ' ' -f 4)  # extract number from jID


## vcf recalibrated files
jID4=$(sbatch --array=0-$(( ${#bams[@]} - 1 )) --dependency=afterok:$jID3 haplotype_caller.sbs $bams_path $reference_fa yes)
jID4=$(echo $jID4 | cut -d ' ' -f 4)  # extract number from jID

jID5=$(sbatch --array=0-$(( ${#bams[@]} - 1 )) --dependency=afterok:$jID4 base_recalibrator.sbs $bams_path ./vcf $reference_fa yes)
jID5=$(echo $jID5 | cut -d ' ' -f 4)  # extract number from jID

## Report Covariates construction
sbatch --dependency=afterok:$jID5 analyze_covariates.sbs ./stats_reports $samples_dir/dna 
