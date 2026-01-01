# SSTDhunter v1.0.0
A tool for investigating steroid-17,20-desmolase (SSTD) encoding genes in metagenomic data

## Introduction
    Steroid-17,20-desmolase (SSTD) activity is believed to be responsible for the biotransformation of cortisol into 11-oxy-androgen. It is also endowed with the ability to cleave the side chain of other endogenous steroids and pharmaceutical glucocorticoids, resulting in similar 11-oxy-androgens. For instance, the side-chain cleavage product of prednisone can significantly promote the proliferation of prostate cancer cells.
    Here, we proposed SSTDhunter, a specialized gene database and associated tool designed for investigating SSTD-encoding genes in metagenomic data. The database is constructed upon large-scale genomic surveys and analyses of SSTD operons. SSTDhunter relies on Bowtie2 or similar tools (such as BWA) to generate alignment files in SAM format. The associated tool, "SSTDhunter.pl", is written in Perl and generates the abundance of SSTD-encoding genes in RPKM.

## Recommended pipeline:
Here is the recommended pipeline for metagenomes investigation under linux environment.

## Requirements
- Linux env
- Bowtie2 v2.5.4
- Perl v5.26.2

### Step 1: Build the index using bowtie2-build
$ bowtie2-build $path/SSTDhunterDB.fa SSTDhunterDB

### Step 2: Run a simple alignment
$ bowtie2 -p 30 -x $path/SSTDhunterDB -1 $path/[prefix]_R1.fq.gz -2 $path/[prefix]_R2.fq.gz -S [prefix].sam

### Step 3: Calculate the abundance of SSTD-encoding genes
$ perl $path/SSTDhunter.pl [prefix].sam

### Tips:
The "$path" represents the readable path of the corresponding file.
The "[prefix]" represents the prefix of the metagenome sample.

### Output files
[prefix]_SSTDhunter_abundance.txt: Contains the abundance (RPKM) of each SSTD-encoding gene.
[prefix]_SSTDhunter_reads.txt: Contains the reads of SSTD-encoding genes in FASTA format.
