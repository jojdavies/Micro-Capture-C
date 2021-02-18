# Micro-Capture-C

This is the code used to make a pipeline for analysing Micro Capture-C (MCC) data.

There are three main scripts for analysing the data, which can be used independently to generate output files and the whole process can be run manually.

1. MCC_BLATfa.pl
2. MCC_fastq_splitter.pl
3. MCC_analyser.pl

These scripts are subject to academic only use and can be downloaded from the Oxford University Software Store https://process.innovation.ox.ac.uk/software/p/16529a/micro-capture-c-academic/1


Prior to running the analysis a reference file of the sequence surrounding the oligo positions needs to be generated using the MCC_BLATfa.pl script 

The basic scheme of analysis is as follows:
1. Trim adaptors (Trim_galore Babraham Institute https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
2. Merge reads into a single read using the central area of overlap (FLASh https://ccb.jhu.edu/software/FLASH/)
3. Convert FASTQ to FA (I used sed)
4. Use BLAT (from ucsctools) to map reads to the 800bp around the capture site
5. Use MCC_fastq_splitter.pl to split the fastqs into multiple files, one for each view point
6. Bowtie to align these files
7. MCC_analyser.pl to analyse these files

Available here are two wrapper scripts (imaginatively called MCC_pipe1.pl and MCC_pipe2.pl) that you can use to automate analysis of the data but you may need to modify these to get them to run on your system

This pipeline has a number of dependencies including the following:
1. Trim_galore
2. FLASh
3. BLAT
4. Bowtie2
5. Samtools
6. MACS2
7. WigToBigWig (ucsctools) 
8. Sun Grid engine

The pipeline will also make you a trackhub to view the data. 
