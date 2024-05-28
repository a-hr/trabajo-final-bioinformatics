#!/usr/bin/env bash

#SBATCH --job-name=NGS_CNV_1
#SBATCH --output=NGS_CNV_1.out
#SBATCH --time=00:15:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

mkdir -p NGS_CNV_1

# ---- Inputs ----
SAMPLES_DIR=/scratch/heral/bioinfo/CNV
REFERENCE=$SAMPLES_DIR/data/GRCh38.primary_assembly.genome.fa
DICT_REFERENCE=$SAMPLES_DIR/data/GRCh38.primary_assembly.genome.dict
OUTPUT_DIR=NGS_CNV_1
TUMOR_FOLDER_1=$SAMPLES_DIR/data/tutorial_11682
INTERVALS=$TUMOR_FOLDER_1/targets_C.interval_list
TUMOR=$TUMOR_FOLDER_1/tumor.bam

GATK=bin/gatk-4.3.0.0/gatk

# ---- Step 1: Collect raw counts data ----
module load Java/11.0.2

$GATK PreprocessIntervals \
    -L $INTERVALS \
    -R $REFERENCE \
    --bin-length 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O $OUTPUT_DIR/targets_C.preprocessed.interval_list

$GATK CollectReadCounts \
    -I $TUMOR \
    -L $INTERVALS \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O $OUTPUT_DIR/tumor.counts.hdf5

# ---- Step 2: Create a CNV panel of normals ----
# $GATK CreateReadCountPanelOfNormals \
#     -I data/tutorial_11682/HG00133.alt_bwamem_GRCh38DH.20150826.GBR.exome.counts.hdf5 \
#     -I data/tutorial_11682/HG00733.alt_bwamem_GRCh38DH.20150826.PUR.exome.counts.hdf5 \
#     -I data/tutorial_11682/NA19654.alt_bwamem_GRCh38DH.20150826.MXL.exome.counts.hdf5 \
#     --minimum-interval-median-percentile 5.0 \
#     -O $OUTPUT_DIR/cnvponC.pon.hdf5

#---- Step 3: Standardize and denoise case read counts against the PoN ----
$GATK DenoiseReadCounts \
    -I data/tutorial_11682/hcc1143_T_clean.counts.hdf5 \
    --count-panel-of-normals data/tutorial_11682/cnvponC.pon.hdf5 \
    --standardized-copy-ratios $OUTPUT_DIR/hcc1143_T_clean.standardizedCR.tsv \
    --denoised-copy-ratios $OUTPUT_DIR/hcc1143_T_clean.denoisedCR.tsv

#4 Plot standardized and denoised copy ratios
## (PlotDenoisedCopyRatios requires packages optparse and data.table)
# $GATK PlotDenoisedCopyRatios \
#     --standardized-copy-ratios $OUTPUT_DIR/hcc1143_T_clean.standardizedCR.tsv \
#     --denoised-copy-ratios $OUTPUT_DIR/hcc1143_T_clean.denoisedCR.tsv \
#     --sequence-dictionary $DICT_REFERENCE \
#     --minimum-contig-length 46709983 \
#     --output $OUTPUT_DIR/plots \
#     --output-prefix hcc1143_T_clean
