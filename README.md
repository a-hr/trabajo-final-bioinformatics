# Bioinformatics final project

> Andoni Berbel & Álvaro Herrero, 2024/05

- [Bioinformatics final project](#bioinformatics-final-project)
  - [Introduction](#introduction)
  - [Module 01 - Quality control](#module-01---quality-control)
  - [Module 02 - Read alignment](#module-02---read-alignment)
  - [Module 03 - Variant calling](#module-03---variant-calling)
  - [Module 04 - CNV](#module-04---cnv)


## Introduction

This repository contains the final project by Álvaro Herrero and Andoni Berbel for the High Performance Computing and Bioinformatics subjects of our Master's degree. 

The project involves a comprehensive analysis of RNAseq and whole genome DNAseq data. Our pipeline includes several key bioinformatics procedures:

1. **Quality Control (QC)**: We first ensure the quality of our data using FastQC and MultiQC for a comprehensive report.

2. **Read Alignment**: We align the RNA sequences to the reference genome using Kallisto, and the DNA sequences using BWA (Burrows-Wheeler Aligner).

3. **Variant Calling**: We identify variants in the DNA sequences using the Genome Analysis Toolkit (GATK).

4. **Copy Number Variation (CNV) Analysis**: Finally, we analyze copy number variations in the DNA sequences.

Each of these steps is encapsulated in a module, as detailed in the following sections.

## Module 01 - Quality control

Starting from the FASTQ files, parallelized `FASTQC`is performed to both RNA and DNA samples. The results are then merged using `MultiQC` to generate a comprehensive report.

All the subprocesses are launched from the `launch.sh` script, which manages their dependency relationships through SLURM job dependencies:

```
#!/usr/bin/env bash

cd /scratch/$USER/trabajo-final-bioinformatics/01_QC

mkdir -p logs

samples_path=/scratch/$USER/trabajo-final-bioinformatics/samples

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
```
The script starts by changing to the working directory and creating a logs directory. It sets the path for the sample files and begins the RNA quality control (QC) process by listing all gzipped RNA samples and creating an output directory for FastQC results. It submits a batch job array to run FastQC on the RNA samples, extracts the job ID, and then submits a MultiQC job to aggregate the FastQC results, ensuring it runs after all FastQC jobs complete.

The script then performs the same steps for DNA samples: listing the gzipped DNA samples, creating the output directory, submitting a batch job array for FastQC, extracting the job ID, and submitting a MultiQC job to aggregate the FastQC results after all FastQC jobs are finished.

SLURM job to run `FASTQC`:

```
#!/bin/bash

#SBATCH --job-name=fastqc
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --time=00:25:00
#SBATCH -o ./logs/%x_%A_%a.out

samples=$1
fastqc_output=$2

files=($(realpath $samples/*.fastq.gz))
file=${files[$SLURM_ARRAY_TASK_ID]}

module load Python
conda activate /scratch/$USER/envs/NGS

fastqc -o $2 $file
```

The script configures a SLURM job named "fastqc" with specific resource allocations and logs output. It takes the sample directory and FastQC output directory as arguments, creates an array of FASTQ files, and selects one based on the SLURM array task ID. It loads the Python module, activates a conda environment, and runs FastQC on the selected file, saving the results to the specified output directory.

SLURM job to run `MULTIQC`:

```
#!/bin/bash

#SBATCH --job-name=multiqc
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --time=00:25:00
#SBATCH -o ./logs/%x_%A_%a.outb

module load Python
conda activate /scratch/$USER/envs/NGS

multiqc $1 -n $2
```

The script sets up a SLURM job named "multiqc" with specific resource allocations and logs output. It loads the Python module, activates a conda environment, and runs MultiQC on the specified input directory, naming the report with the provided output name.

## Module 02 - Read alignment

Starting from the `FASTQ` files, RNA and DNA samples are alligned using Kallisto and BWA respectively.

Again all the subprocesses are launched from the `launch.sh` script:

```
#!/usr/bin/env bash

cd /scratch/$USER/trabajo-final-bioinformatics/02_align
mkdir -p logs

samples_path=/scratch/$USER/trabajo-final-bioinformatics/samples
reference_fa=/scratch/$USER/trabajo-final-bioinformatics/samples/references/GRCh38.primary_assembly.genome.fa
reference_fa_transcriptome=$(dirname $reference_fa)/gencode.v38.transcripts.fa.gz
reference_gtf=/scratch/$USER/trabajo-final-bioinformatics/samples/references/gencode.v41.primary_assembly.annotation.gtf

bootstraps=40

# ---- Align RNA: Kallisto ----

samples_rna=( $(ls $samples_path/rna/*_1.fastq.gz) )
outdir=./rna_bams

# download transcriptome if not present
if [ ! -f $reference_fa_transcriptome ]
then
    curl http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz \
        -o $reference_fa_transcriptome
    # chmod 777 human.fa.gz
fi


# check if kallisto index exists
kallisto_index=$(dirname $reference_fa_transcriptome)/kallisto_index.idx

if [ -f $(dirname $reference_fa_transcriptome)/kallisto_index.idx ]
then
    # kallisto index already exists
    sbatch --array=0-$((${#samples_rna[@]} - 1))%5 \
        align_RNA.sbs \
        $kallisto_index \
        $samples_path/rna \
        $outdir \
        $bootstraps
else
    # create index
    jID1=$(sbatch create_index.sbs $reference_fa_transcriptome kallisto)
    jID1=$(echo $jID1 | cut -d ' ' -f 4)  # extract number from jID

    # align after index is created
    kallisto_index=$(dirname $reference_fa_transcriptome)/kallisto_index.idx
    sbatch --array=0-$((${#samples_rna[@]} - 1))%5 \
        --dependency=afterok:$jID1 \
        align_RNA.sbs \
        $kallisto_index \
        $samples_path/rna \
        $outdir \
        $bootstraps
fi

# ---- Align DNA: BWA ----

samples_dna=( $(ls $samples_path/dna/*_1.fastq.gz) )

# check if bwa index exists
if [ -f ${reference_fa}.bwt ]
then
    # bwa index already exists -> align
    sbatch --array=0-$((${#samples_dna[@]} - 1))%5 \
        align_DNA.sbs \
        $samples_path/dna \
        $reference_fa
else
    # create index
    jID2=$(sbatch create_index.sbs $reference_fa bwa)
    jID2=$(echo $jID2 | cut -d ' ' -f 4)  # extract number from jID

    # align
    sbatch --dependency=afterok:$jID2 \
        --array=0-$((${#samples_dna[@]} - 1))%5 \
        align_DNA.sbs \
        $samples_path/dna \
        $reference_fa
fi
```

This script automates the alignment of RNA and DNA samples to their respective reference genomes. It sets up the environment, defines paths for data and reference genomes, and specifies parameters like bootstraps.

For RNA alignment using Kallisto, it checks for the transcriptome file and the Kallisto index. If absent, it downloads the transcriptome and creates the index before aligning the samples.

For DNA alignment using BWA, it checks for the BWA index and creates it if necessary before aligning the samples.

SLURM job to create indexes:

```
#!/bin/bash

#SBATCH --job-name=index
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH -o ./logs/%x_%A_%a.out

module load Python
conda activate /scratch/$USER/envs/NGS

reference_fa=$1
aligner=$2 # kallisto, bwa

if [ $aligner == "kallisto" ]; then
    kallisto index -i $(dirname $reference_fa)/kallisto_index.idx $reference_fa
elif [ $aligner == "bwa" ]; then
    bwa index $reference_fa
else
    echo "Invalid aligner"
fi
```

This script generates indices for either Kallisto or BWA aligners based on the provided reference genome file and aligner type. It sets up SLURM job configurations, including resource allocations and log output location. Then, it loads the Python module and activates a Conda environment. If the aligner is "kallisto", it generates a Kallisto index using the reference genome file and stores it in the same directory with the name kallisto_index.idx. Otherwise, it generates a BWA index directly from the reference genome file.

SLURM job to align DNA:

```
#!/bin/bash

#SBATCH --job-name=BWA
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH -o ./logs/%x_%A_%a.out

module load Python
conda activate /scratch/$USER/envs/NGS

samples_path=$1
samples1=($(realpath $samples_path/*_1.fastq.gz))

reference_fa=$2 # reference fasta file, the index will be along side it

i=$SLURM_ARRAY_TASK_ID
fq1=${samples1[$i]}
fq2=${fq1/_1/_2}

mkdir -p dna_bams stats_reports_DNA

samplename=$(basename $fq1 _1.fastq.gz).bam

bwa mem -M -t 8 \
	-R "@RG\tID:ind1\tSM:ind1\tLB:ind1\tPU:ind1\tPL:Illumina" \
	$reference_fa \
	$fq1 $fq2 \
	| samblaster -M | samtools fixmate - - | samtools sort -O bam \
	-o dna_bams/$samplename

samtools index dna_bams/$samplename
samtools stats dna_bams/$samplename | grep ^SN \
    | cut -f 2- \
    > stats_reports_DNA/$(basename $samplename .bam).stats
```

This script aligns DNA samples using BWA. It configures SLURM job settings, loads required modules, and activates a Conda environment. It takes input parameters for sample data and the reference fasta file. After listing sample files, it aligns reads using BWA mem, processes the output with samblaster and samtools, saves the resulting BAM file, and generates alignment statistics.

SLURM job to align RNA:

```
#!/bin/bash

#SBATCH --job-name=kallisto
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --time=00:15:00
#SBATCH -o ./logs/%x_%A_%a.out

module load Python
conda activate /scratch/$USER/envs/NGS

index=$1
fastqs_path=$2
outdir=$3
bootstraps=$4

# Select the fastq files
fastqs=($(realpath $fastqs_path/*_1.fastq.gz))
fq1=${fastqs[$SLURM_ARRAY_TASK_ID]}
fq2=${fq1/_1/_2}

# Create output directory
out=$outdir/$(basename $fq1 _1.fastq.gz)    
mkdir -p $out

# Run kallisto
kallisto quant -t 10 -i $index -b $bootstraps -o $out $fq1 $fq2
```

This script configures the SLURM job with specific resource allocations and log output settings. After loading necessary modules and activating a Conda environment, the script takes input parameters, including the Kallisto index path, the directory containing the fastq files, the output directory where results will be saved, and the number of bootstraps for statistical analysis. It then selects the forward read fastq files and creates an output directory for each sample. Finally, it executes Kallisto with the provided options, including multi-threading, using the index, and saving the quantification results in the designated output directory.

## Module 03 - Variant calling

```
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
```

This script manages variant calling and base recalibration processes for DNA samples. It begins by defining the reference genome file path and the path to the BAM files. If the index file for the reference genome does not exist, it is created first. Variant calling is then performed on each BAM file with dependency on the index creation job. Base recalibration follows variant calling, and recalibration is applied to the VCF files generated earlier. Additional variant calling and base recalibration are conducted on the recalibrated files. Finally, a job is submitted to analyze covariates for constructed reports, dependent on the completion of the recalibration jobs.

SLURM job to construct indexes:

```
#!/bin/bash

#SBATCH --job-name=SeqDict
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH -o ./logs/%x_%A_%a.out

: '
Extra dependencies:

wget https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar
wget https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
unzip gatk-4.3.0.0.zip

'

PICARD="/scratch/$USER/trabajo-final-bioinformatics/bin/picard.jar"

module load Java/11.0.2

reference_fa=$1

parent_dir=$(dirname $reference_fa)

java -jar $PICARD CreateSequenceDictionary \
    -R $reference_fa \
    -O $parent_dir/$(basename $reference_fa .fa).dict

# create index file (.fai) for reference fasta
module load Python
conda activate /scratch/$USER/envs/NGS

samtools faidx $reference_fa
```

This script generates a sequence dictionary and index file for a reference genome in fasta format. It is configured as a SLURM job with resource allocations and log output settings. First, it defines the path to the Picard tool and loads the Java module. Then, it takes the reference genome file path as an argument. The script extracts the directory path and file name from the reference genome path. Using Picard, it creates a sequence dictionary (.dict) file in the parent directory of the reference genome file. Additionally, it activates a Conda environment, loads the Python module, and uses samtools to create an index file (.fai) for the reference fasta.

SLURM job to perform base recalibration:

```
#!/bin/bash

#SBATCH --job-name=base_recalibrator
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH -o ./logs/%x_%A_%a.out

module load Java/11.0.2

GATK=/scratch/$USER/trabajo-final-bioinformatics/bin/gatk-4.3.0.0/gatk

bams_path=$1
vcfs_path=$2
reference_fa=$3
recalibrated=$4  # yes or no

if [ $recalibrated == "yes" ]; then
    bam_files=($(ls $bams_path/*recalibrated.bam))
    i=$SLURM_ARRAY_TASK_ID
    bam=${bam_files[i]}

    outname=$(basename $bam .bam).table

elif [ $recalibrated == "no" ]; then
    bam_files=($(ls $bams_path/*.bam | grep -v recalibrated))
    i=$SLURM_ARRAY_TASK_ID
    bam=${bam_files[i]}

    outname=$(basename $bam .bam).table

fi

vcf_files=($(ls $vcfs_path/*.vcf))
vcf=${vcf_files[i]}


mkdir -p stats_reports

$GATK BaseRecalibrator \
    -I $bam \
	-R $reference_fa \
	--known-sites $vcf \
	-O stats_reports/$outname
```

This script conducts base recalibration using GATK (Genome Analysis Toolkit) for BAM files. It is configured as a SLURM job with resource allocations and log output settings. The script loads the Java module and specifies the path to GATK. Input parameters include paths to BAM files, VCF files, reference fasta file, and a flag indicating whether the BAM files have been recalibrated ("yes" or "no"). Depending on the value of the "recalibrated" flag, the script selects BAM files accordingly. It also selects VCF files based on the current SLURM task ID. Base recalibration is then performed using GATK's BaseRecalibrator tool, utilizing the selected BAM and VCF files, and saving the results in the "stats_reports" directory. 

SLURM job to apply recalibration:

```
#!/bin/bash

#SBATCH --job-name=apply_recalibrator
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH -o ./logs/%x_%A_%a.out

module load Java/11.0.2

GATK=/scratch/$USER/trabajo-final-bioinformatics/bin/gatk-4.3.0.0/gatk

bams_path=$1
tables_path=$2

bam_files=($(ls $bams_path/*.bam))
i=$SLURM_ARRAY_TASK_ID
bam=${bam_files[i]}

table_files=($(ls $tables_path/*.table))
table=${table_files[i]}

outname=$(basename $bam .bam)_recalibrated.bam

$GATK ApplyBQSR \
    -I $bam \
	--bqsr $table \
	-O $bams_path/$outname
```

This script applies base quality score recalibration (BQSR) using GATK (Genome Analysis Toolkit) to BAM files. It is configured as a SLURM job with resource allocations and log output settings. The script loads the Java module and specifies the path to GATK. Input parameters include paths to BAM files and tables generated from base recalibration. The script selects BAM and table files based on the current SLURM task ID. Base recalibration is then applied using GATK's ApplyBQSR tool, utilizing the selected BAM and table files, and saving the recalibrated BAM files in the specified directory.

SLURM job to haplotype caller:

```
#!/bin/bash

#SBATCH --job-name=VCF
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH -o ./logs/%x_%A_%a.out

module load Java/11.0.2

GATK=/scratch/$USER/trabajo-final-bioinformatics/bin/gatk-4.3.0.0/gatk

bams_path=$1
reference_fa=$2
recalibrated=$3  # yes or no

if [ $recalibrated == "yes" ]; then
    bam_files=($(ls $bams_path/*recalibrated.bam))
    i=$SLURM_ARRAY_TASK_ID
    bam=${bam_files[i]}
    outname=$(basename $bam .bam).vcf

elif [ $recalibrated == "no" ]; then
    bam_files=($(ls $bams_path/*.bam | grep -v recalibrated))
    i=$SLURM_ARRAY_TASK_ID
    bam=${bam_files[i]}
    outname=$(basename $bam .bam).vcf
fi

mkdir -p vcf
$GATK HaplotypeCaller \
    -I $bam \
    -R $reference_fa \
    -O vcf/$outname
```

This script conducts variant calling using GATK's HaplotypeCaller for BAM files. It is configured as a SLURM job with resource allocations and log output settings. The script loads the Java module and specifies the path to GATK. Input parameters include paths to BAM files, reference fasta file, and a flag indicating whether the BAM files have been recalibrated ("yes" or "no"). Depending on the value of the "recalibrated" flag, the script selects BAM files accordingly. Variant calling is then performed using GATK's HaplotypeCaller tool, utilizing the selected BAM files and saving the VCF files in the "vcf" directory.

SLURM job to analyze covariates:

```
#!/bin/bash

#SBATCH --job-name=report
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH -o ./logs/%x_%A_%a.out

module load Java/11.0.2 R

GATK=/scratch/$USER/trabajo-final-bioinformatics/bin/gatk-4.3.0.0/gatk

tables_dir=$1
samples_dir=$2

fastqs=($(ls $samples_dir/*_1.fastq.gz))
sample=${fastqs[$SLURM_ARRAY_TASK_ID]}
sample=$(basename $sample _1.fastq.gz)

table_before=$(ls $tables_dir/$sample*.table | grep -v recalibrated)
table_after=$(ls $tables_dir/$sample*.table  | grep recalibrated)

mkdir -p reports

$GATK AnalyzeCovariates \
	--before $table_before \
	--after $table_after \
	--plots reports/${sample}_report.pdf
```
This script generates reports on base quality score recalibration using GATK's AnalyzeCovariates tool. It is configured as a SLURM job with resource allocations and log output settings. The script loads the Java and R modules and specifies the path to GATK. Input parameters include paths to directories containing tables generated before and after recalibration and directories containing sample fastq files. The script selects the fastq files based on the current SLURM task ID and extracts the sample name. It then identifies the tables corresponding to the sample before and after recalibration. Reports are generated using GATK's AnalyzeCovariates tool, comparing the tables before and after recalibration, and saving the plots in the "reports" directory

## Module 04 - CNV

2 different pipelines were used to perform CNV analysis over the bams. Both pipelines are called by a launch script, and each of them is in a different sbash script.

The launch script:

```
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

# launch fastqc
jID_=$(sbatch --array=0-$((${#bams[@]} - 1))%5 pipeline2.sbs $bams_path)
# extract number from jID
jID_=$(echo $jID_ | cut -d ' ' -f 4)
```

SLURM job of first pipeline:

```
#!/bin/bash

#SBATCH --job-name=report
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH -o ./logs/%x_%A_%a.out

REFERENCE=/scratch/$USER/trabajo-final-bioinformatics/samples/references/GRCh38.primary_assembly.genome.fa
DICT_REFERENCE=/scratch/$USER/trabajo-final-bioinformatics/samples/references/GRCh38.primary_assembly.genome.dict
OUTPUT_DIR=/scratch/$USER/trabajo-final-bioinformatics/04_CNV/Output
INTERVALS=/scratch/$USER/trabajo-final-bioinformatics/samples/tutorial_11682/targets_C.interval_list
GATK=/scratch/$USER/trabajo-final-bioinformatics/bin/gatk-4.3.0.0/gatk

bams_path=$1
bams=($(ls $bams_path/*_recalibrated.bam))
TUMOR=${bams[$SLURM_ARRAY_TASK_ID]}
name=$(basename $TUMOR _recalibrated.bam)

## STEPS:

module load Java/11.0.2
module load R/3.6.2-intel-2019a

#1 Collect raw counts data

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
    -O $OUTPUT_DIR/$name.counts.hdf5

#2 Generate a CNV panel of normals
# (warning)

#$GATK CreateReadCountPanelOfNormals \
#    -I tutorial_11682/HG00133.alt_bwamem_GRCh38DH.20150826.GBR.exome.counts.hdf5 \
#    -I tutorial_11682/HG00733.alt_bwamem_GRCh38DH.20150826.PUR.exome.counts.hdf5 \
#    -I tutorial_11682/NA19654.alt_bwamem_GRCh38DH.20150826.MXL.exome.counts.hdf5 \
#    --minimum-interval-median-percentile 5.0 \
#    -O $OUTPUT_DIR/cnvponC.pon.hdf5

#3 Standardize and denoise case read counts against the PoN

$GATK DenoiseReadCounts \
    -I $OUTPUT_DIR/$name.counts.hdf5 \
    --count-panel-of-normals /scratch/$USER/trabajo-final-bioinformatics/samples/tutorial_11682/cnvponC.pon.hdf5 \
    --standardized-copy-ratios $OUTPUT_DIR/$name.standardizedCR.tsv \
    --number-of-eigensamples 10 \
    --denoised-copy-ratios $OUTPUT_DIR/$name.denoisedCR.tsv


#4 Plot standardized and denoised copy ratios
## (PlotDenoisedCopyRatios requires packages optparse and data.table)



$GATK PlotDenoisedCopyRatios \
    --standardized-copy-ratios $OUTPUT_DIR/$name.standardizedCR.tsv \
    --denoised-copy-ratios $OUTPUT_DIR/$name.denoisedCR.tsv \
    --sequence-dictionary $DICT_REFERENCE \
    --minimum-contig-length 46709983 \
    --output $OUTPUT_DIR/plots \
    --output-prefix $name
```

This script performs CNV analysis using GATK. It collects raw counts, denoises and standardizes them against a panel of normals, and plots the resulting copy ratios. The script is configured as a SLURM job with resource allocations. It includes commented-out code for generating a panel of normals.

SLURM job of second pipeline:

```
#!/bin/bash

#SBATCH --job-name=report
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=16G
#SBATCH --time=10:00:00
#SBATCH -o ./logs/%x_%A_%a.out

## INPUTS:


REFERENCE=/scratch/$USER/trabajo-final-bioinformatics/samples/references/GRCh38.primary_assembly.genome.fa
DICT_REFERENCE=/scratch/$USER/trabajo-final-bioinformatics/samples/references/GRCh38.primary_assembly.genome.dict
INTERVALS=/scratch/$USER/trabajo-final-bioinformatics/samples/tutorial_11683/chr17_theta_snps.interval_list
#TUMOR=/scratch/aberbelgara/NGS22/samples/tutorial_11683/tumor.bam
GATK=/scratch/$USER/trabajo-final-bioinformatics/bin/gatk-4.3.0.0/gatk
NORMAL=/scratch/$USER/trabajo-final-bioinformatics/samples/tutorial_11683/normal.bam
OUTPUT_DIR=/scratch/$USER/trabajo-final-bioinformatics/04_CNV/Output2

bams_path=$1
bams=($(ls $bams_path/*_recalibrated.bam))
TUMOR=${bams[$SLURM_ARRAY_TASK_ID]}
name=$(basename $TUMOR _recalibrated.bam)

## STEPS:
module load Java/11.0.2
module load R/3.6.2-intel-2019a

#1 Count ref and alt alleles at common germline variant sites

$GATK --java-options "-Xmx3g" CollectAllelicCounts \
    -L $INTERVALS \
    -I $NORMAL \
    -R $REFERENCE \
    -O $OUTPUT_DIR/$name _N.allelicCounts.tsv


$GATK --java-options "-Xmx3g" CollectAllelicCounts \
    -L $INTERVALS \
    -I $TUMOR \
    -R $REFERENCE \
    -O $OUTPUT_DIR/$name _T.allelicCounts.tsv

#2 Group contiguous copy ratios into segments

$GATK --java-options "-Xmx4g" ModelSegments \
--denoised-copy-ratios /scratch/$USER/trabajo-final-bioinformatics/04_CNV/Output/$name.denoisedCR.tsv \
    --allelic-counts $OUTPUT_DIR/$name _T.allelicCounts.tsv \
    --normal-allelic-counts $OUTPUT_DIR/$name _N.allelicCounts.tsv \
    --output $OUTPUT_DIR \
    --output-prefix $name


#3 Call copy-neutral (Amplified and deleted segments)

# This step is not required for plotting segmentation results.
# The resulting called.seg data adds the sixth column to the provided copy
# ratio segmentation table.
# The tool denotes amplifications with a + plus sign, deletions with a - minus sign and neutral segments with a 0 zero.

$GATK CallCopyRatioSegments \
    --input $OUTPUT_DIR/$name _T.cr.seg \
    --output $OUTPUT_DIR/$name _T.called.seg


#4 Plot modeled copy ratio and allelic fraction segments
## (PlotDenoisedCopyRatios requires packages optparse and data.table)

$GATK PlotModeledSegments \
    --denoised-copy-ratios /scratch/$USER/trabajo-final-bioinformatics/04_CNV/Output/$name.denoisedCR.tsv \
    --allelic-counts $OUTPUT_DIR/$name.hets.tsv \
    --segments $OUTPUT_DIR/$name.modelFinal.seg \
    --sequence-dictionary $DICT_REFERENCE \
    --minimum-contig-length 46709983 \
    --output $OUTPUT_DIR/plots \
    --output-prefix $name
```