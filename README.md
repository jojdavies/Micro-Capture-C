# Micro-Capture-C

This is the code used to make a pipeline for analysing Micro Capture-C (MCC) data.

There are three main scripts for analysing the data, which can be used independently to generate output files and the whole process can be run manually.

1. MCC_oligo_fa.pl
2. MCC_fastq_splitter.pl
3. MCC_analyser.pl

These scripts are subject to academic only use and can be downloaded from the Oxford University Software Store https://process.innovation.ox.ac.uk/software/p/16529a/micro-capture-c-academic/1
Please download these scripts and the MCC_pipe scripts all into the same folder in order for it to run properly.

Prior to running the analysis a reference file of the sequence surrounding the oligo positions needs to be generated using the MCC_oligo_fa.pl script 
It also creates a bed file of the oligos and you can alter this to specify the colours of the tracks in the track hub at the end if you like.

The basic scheme of analysis is as follows:
1. Trim adaptors (Trim_galore Babraham Institute https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
2. Merge reads into a single read using the central area of overlap (FLASh https://ccb.jhu.edu/software/FLASH/)
3. Convert FASTQ to FA (I used sed)
4. Use BLAT (from ucsctools) to map reads to the 800bp around the capture site. This uses the fasta file generated by MCC_oligo_fa.pl as the reference to map sequences to
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

To run the pipeline initially run the first script MCC_pipe1.pl to make a config file called MCC_pipe_config.txt in the directory which contains the fastq files you want to analyse

> perl /path_to_scripts_folder/MCC_pipe1.pl -config
 
This will generated a file called MCC_pipe_config.txt, which you need to edit to specify the, paths to the genome files for bowtie; the names of the fastq files and the project name and the paths to the oligo files. 

It also allows you to specify all of the different parts of the pipeline you want to run.
Initially I would recommend running the first part of the pipe and then running the second part once this has completed.

To run the first part of the pipe change the following in the MCC_pipe_config.txt file

>trim_galore=Y

>flash=Y

>fq2fa=Y

>blat=Y

>MCC_splitter=Y

but leave 
>Pipe2=N

To run the pipeline type 
>perl /path_to/MCC_pipe1.pl

This will generate a tmp.sh file with the commands to execute. 
It is worthwhile checking that this seems sensible before executing this!

If you run the script with the -q flag it will execute the tmp.sh script on slurm (you may need to alter the exact command in the MCC_pipe1.pl and MCC_pipe2.pl scripts to adapt it to your local server). 

Note that if you have multiple experiments the tmp.sh file is overwritten each time and so you will only see the command run for the final files.

The MCC_pipe_config.txt has a small bash script at the top, which allows you to run the pipeline (with the -q flag) by typing:

>bash MCC_pipe_config.txt

If this has run correctly it should create a folder which has a fastq file in it for each of your targets. 
You can then run the second part of the pipe by changing the initial parts from "Y" to "N" in the MCC_pipe_config.txt and changing the following:

>Pipe2=Y

>bowtie=Y

>sort=Y

>MCCanal=Y

>track_hub=Y

The options gzip, postgzip, clean_wig and samtobam allow you to compress and clean up the files quickly

The macs2 option will call MACS2 for simple peak calling
