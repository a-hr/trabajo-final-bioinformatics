#!/usr/bin/env bash

cd /scratch/heral/TRABAJO_FINAL/01_QC

mkdir -p logs

samples_path=/scratch/heral/TRABAJO_FINAL/samples

# ---- RNA QC ----
samples=( $(ls $samples_path/rna/*.gz) )
fastqc_output=$samples_path/rna/fastqc_output

mkdir -p $fastqc_output

# launch fastqc
jID=$(sbatch --array=0-$((${#samples[@]} - 1))%5 run_fastqc.sbs $samples_path/rna $fastqc_output)
# extract number from jID
jID=$(echo $jID | cut -d ' ' -f 4)

# launch multiqc
sbatch --dependency=afterok:$jID run_multiqc.sbs $fastqc_output multiqc_RNA

# ---- DNA QC ----
samples=( $(ls $samples_path/dna/*.gz) )
fastqc_output=$samples_path/dna/fastqc_output

mkdir -p $fastqc_output

# launch fastqc
jID=$(sbatch --array=0-$((${#samples[@]} - 1))%5 run_fastqc.sbs $samples_path/dna $fastqc_output)
# extract number from jID
jID=$(echo $jID | cut -d ' ' -f 4)

# launch multiqc
sbatch --dependency=afterok:$jID run_multiqc.sbs $fastqc_output multiqc_DNA