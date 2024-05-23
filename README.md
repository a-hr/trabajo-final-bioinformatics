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

All the subprocesses are launched from the `launch.sh` script, which manages their dependency relationships through SLURM job dependencies.

## Module 02 - Read alignment


## Module 03 - Variant calling

## Module 04 - CNV

