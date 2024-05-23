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



