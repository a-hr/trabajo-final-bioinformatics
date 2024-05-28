#!/usr/bin/env bash

#SBATCH --job-name=NGS_CNV_2
#SBATCH --output=NGS_CNV_2.out
#SBATCH --time=00:20:00
#SBATCH --mem=3G
#SBATCH --cpus-per-task=1

# ---- Inputs ----
GATK=bin/gatk-4.3.0.0/gatk

SAMPLES_DIR=/scratch/heral/bioinfo/CNV
REFERENCE=$SAMPLES_DIR/data/GRCh38.primary_assembly.genome.fa
DICT_REFERENCE=$SAMPLES_DIR/data/GRCh38.primary_assembly.genome.dict
OUTPUT_DIR=NGS_CNV_2
TUMOR_FOLDER_1=$SAMPLES_DIR/data/tutorial_11682
TUMOR_FOLDER_2=$SAMPLES_DIR/data/tutorial_11683
INTERVALS=$TUMOR_FOLDER_1/targets_C.interval_list
TUMOR=$TUMOR_FOLDER_2/tumor.bam
NORMAL=$TUMOR_FOLDER_2/normal.bam


## STEPS:
#1 Count ref and alt alleles at common germline variant sites

$GATK --java-options "-Xmx3g" CollectAllelicCounts \
    -L $INTERVALS \
    -I $NORMAL \
    -R $REFERENCE \
    -O $OUTPUT_DIR/hcc1143_N_clean.allelicCounts.tsv


$GATK --java-options "-Xmx3g" CollectAllelicCounts \
    -L $INTERVALS \
    -I $TUMOR \
    -R $REFERENCE \
    -O $OUTPUT_DIR/hcc1143_T_clean.allelicCounts.tsv

#2 Group contiguous copy ratios into segments

$GATK --java-options "-Xmx4g" ModelSegments \
    --denoised-copy-ratios $TUMOR_FOLDER_2/hcc1143_T_clean.denoisedCR.tsv \
    --allelic-counts $OUTPUT_DIR/hcc1143_T_clean.allelicCounts.tsv \
    --normal-allelic-counts $OUTPUT_DIR/hcc1143_N_clean.allelicCounts.tsv \
    --output $OUTPUT_DIR \
    --output-prefix hcc1143_T_clean


#3 Call copy-neutral (Amplified and deleted segments)

# This step is not required for plotting segmentation results.
# The resulting called.seg data adds the sixth column to the provided copy
# ratio segmentation table.
# The tool denotes amplifications with a + plus sign, deletions with a - minus sign and neutral segments with a 0 zero.

$GATK CallCopyRatioSegments \
    --input $OUTPUT_DIR/hcc1143_T_clean.cr.seg \
    --output $OUTPUT_DIR/hcc1143_T_clean.called.seg


#4 Plot modeled copy ratio and allelic fraction segments
## (PlotDenoisedCopyRatios requires packages optparse and data.table)

# $GATK PlotModeledSegments \
#     --denoised-copy-ratios tutorial_11683/hcc1143_T_clean.denoisedCR.tsv \
#     --allelic-counts $OUTPUT_DIR/hcc1143_T_clean.hets.tsv \
#     --segments $OUTPUT_DIR/hcc1143_T_clean.modelFinal.seg \
#     --sequence-dictionary $DICT_REFERENCE \
#     --minimum-contig-length 46709983 \
#     --output $OUTPUT_DIR/plots \
#     --output-prefix hcc1143_T_clean
