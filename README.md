# aDNA_data_processing
Repository for data processing scripts common to all LClab projects. 

List of scripts and purpose:

## Alignment 
`seqPro_v7.py`  
Script to trim adapters, align fastqs to a reference and do initial duplicate removal and sample-level merges, as well as calculating alignment statistics.  
Written in python 2.7; assumes cutadapt, AdapterRemoval, samtools, bwa and mapDamage are all accessible in $PATH.

## Alignment Processing
`BamBams_v6.py`  
Script to do final BAM merges, duplicate removal, indel realignment, length and MQ filtering and softclipping.  
`filter_bam.py`  
Script to filter BAM file for length and MQ as well as softclipping.
`softclip.py`  
old script for softclipping.

## Metagenomic analysis

`metascreen.py`  
Script to filter unaligned BAM files to remove duplicates, fish for exact index matches and profile using kraken2 + bracken. Assumes kraken is accessible in $PATH  
`edit_dist.R`   
Script to assess edit distance distribution in microbial alignments.

## Pseudohaploid SNP Calling
`HapSNP_v9.py`  
Script to call pseudohaploid SNPs in ancient DNA data. Assumes plink1.9 accessible in $PATH.
