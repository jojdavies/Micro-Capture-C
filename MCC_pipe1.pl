#!/usr/bin/perl -w

use strict;
use Cwd; use Cwd 'abs_path';
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

=head1 NAME
 
 perl MCC_pipe1.pl 

=head1 SYNOPSIS

This script runs the MCC analysis pipeline - it outputs a tmp.sh file for each sample and executes it with slurm.
The pipeline is in two parts but MCC_pipe1.pl will run MCC_pipe2.pl.
Be careful running MCC_pipe2.pl because it will submit hundreds of jobs onto your cluster (one for each of your targets).

The files must be named {filename}_R1.fastq and {filename}_R2.fastq at the beginning

It is generally better to specify the parameters in the MCC_pipe_config.txt file rather than trying to make a really long command.

The MCC_pipe_config.txt can be created by typing:
perl MCC_pipe1.pl-config

You must specify -q if you want the command to be submitted to qsub.

The pipeline has two components. The first pipeline MCC_pipe1.pl processes the large fastq files and splits these into smaller fastq files in a separate directory depending on which targets they BLAT to.
The second pipeline MCC_pipe2.pl processes the smaller files in the subdirectories.

The pipeline has the following dependencies:
Trim_galore https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/
FLASh https://ccb.jhu.edu/software/FLASH/
BLAT
Bowtie2
Samtools
MACS2
WigToBigWig https://www.encodeproject.org/software/wigtobigwig/
gunzip

=head1 EXAMPLE

perl MCC_pipe1.pl -q

=head1 OPTIONS

 -i input config file if different from MCC_pipe_config.txt - generally it is better to specify everything in the MCC_pipe_config.txt file than specify in the command line.
 -f file names
 -q qsub the commands
 -genome genome build to use
 -pf specify public folder
 -pu specify public url
 -h or -help help!
 -man manual
 -config outputs a blank version of the config file

=head1 AUTHOR

 Written by James Davies 2019

=cut


# Specifications
# Default parameters stored in %para hash
my $pipescript_path=abs_path($0);
my ($pipescript_filename, $pipescript_directory) = fileparse($pipescript_path);
my %param;
$param{"master_folder"} = substr($pipescript_directory,0,-1); # The folder where you plan to store the scripts - defaults to the dirctory this script is stored in
$param{"public_folder"} = ""; # The public folder base directory
$param{"public_url"} = ""; # The name of your public url e.g. sara.molbiol.ox.ac.uk/public/user
$param{"bigwig_folder"} = ""; # The path to your bigwig files - it will expect it to be structured /databank/raw/bigwig for the name of the root followed by /genome for the genome followed by genome_size.txt (e.g. /genomes/bigwig/mm9/mm9_sizes.txt - specifiy here /genomes/bigwig )
$param{"email"} = ''; # Your email - for the UCSC track hub

$param{"reference"}="MCC"; #This is specified in the config file 
$param{'oligo_file'}=""; #Path to the oligo file 
$param{'colour_file'}=""; #Path to a color scheme file if you want one
$param{"genome"}="mm9";
$param{"trim_galore"}="N";
$param{"flash"}="N";
$param{"gunzip"}="N";
$param{"gzip"}="N";
$param{"fq2fa"}="N";
$param{"blat"}="N";
$param{"MCC_splitter"}="N";
$param{"MNase_multi"}="N";
$param{"Align_pipe"}="N";
$param{"Pipe2"}="N";
$param{"qsub"}=0;
$param{"combine"}=0;
$param{"clean_wig"}=0;
$param{"MNase_Multi_limit"}=0;
# Strings
# NB. The following strings should be changed to to your public account, server url and bigwig folder.

my $path = getcwd;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(); #defines the start time of the script (printed out into the report file later on)

my $help=0;
my $man=0;
my $config=0;
my $stringent = 0;
my $input_file="MCC_pipe_config.txt";

# Arrays
my @files;

# The GetOptions from the command line
&GetOptions
(
    "i=s"=>\ $input_file, #input config file
    "f=s"=>\ @files,	# -f1		Input filename 
	"q"=>\ $param{"qsub"},			# -pf		Your public folder (e.g. /hts/data0/public/username)
	"genome=s"=>\ $param{"genome"},
    "pf=s"=>\ $param{"public_folder"},			# -pf		Your public folder (e.g. /hts/data0/public/username)
	"pu=s"=>\ $param{"public_url"},				# -pu		Your public url (e.g. sara.molbiol.ox.ac.uk/public/username)
	'h|help'=>\$help,				# -h or -help 	Help - prints the manual
	'man'=>\$man,					# -man  	prints the manual
    'config'=>\$config, 
    'combine'=>\$param{"combine"}, 
);

pod2usage(1) if $help;
pod2usage(-verbose=>2) if $man;
pod2usage(2) unless ($input_file);

if ($config ==1){config_file()} #outputs the standard config file if requested and exits

# This loads the config file and puts it into the hash
open FH, $input_file or die "Couldn't open input file $input_file";
while (my $line = <FH>)
{
chomp $line;
    if ($line =~ /^#/){next}
    elsif($line =~ /(.*=)(.*)#(.*)/){if(length($2)>5){$line = $1.$2} else{next}} #Removes hased comments from lines
    
    if (($line =~ /^perl/) or ($line =~ /^exit/)){}
    elsif ($line =~ /^\s*$/){}
    elsif (($line =~ /(.*)\s*=\s*(.*)/) or ($line =~ /(.*)\s*=\s*(.*)\s*#.*/))
        {
            if ($1 eq "file"){push @files, $2}
            else{$param{$1}=$2}
        }
    else {print "line of input file not read: $line\n"}
}


#Makes a subdirectory in the public folder
$param{"public_folder_expt"} = $param{"public_folder"}."/".$param{"reference"};
if (-d $param{"public_folder"}){}
else {mkdir $param{"public_folder"}};

print_parameters(\%param);


foreach my $filename (@files)
{
print "Writing tmp.sh for $filename\n";
open TMP, ">tmp.sh" or die "Couldn't open input file tmp.sh";
print TMP "#!/bin/bash
# Script started at $hour:$min:$sec $mday/$mon/$year
# Filename $filename\n";
output_parameters(\%param, \*TMP);

########################################################################
#Makes the output path for MNAse_psl_multi and cleans it of the wig files
my $output_path= "$path/".$param{'reference'}.$filename;
if (-d $output_path){}
else {mkdir $output_path};

# Deletes all of the previous wig files in the folders
if (($param{"clean_wig"} =~ /^[Yy1]/) and ($param{"qsub"} =~ /^[Yy1]/))
{
chdir($output_path);
system("rm *.wig -rf");
chdir($path);
}

# To unzip the files
if ($param{"gunzip"} =~ /^[Yy1]/)
{
    if (-e "$filename\_R1.fq.gz"){print TMP "mv $filename\_R1.fq.gz $filename\_R1.fastq.gz\n"}
    if (-e "$filename\_R2.fq.gz"){print TMP "mv $filename\_R2.fq.gz $filename\_R2.fastq.gz\n"}    
    if (-e "$filename\_R1.fastq.gz"){print TMP "gunzip $filename\_R1.fastq.gz\n"}
    if (-e "$filename\_R2.fastq.gz"){print TMP "gunzip $filename\_R2.fastq.gz\n"}
    if (-e "$filename\_ext.fastq.gz"){print TMP "gunzip $filename\_ext.fastq.gz\n"}
    if (-e "$filename\_ext.fa.gz"){print TMP "gunzip $filename\_ext.fa.gz\n"}
    if (-e "$filename\_ext.psl.gz"){print TMP "gunzip $filename\_ext.psl.gz\n"}
}

# To run trim_galore to cut the reads
if ($param{"trim_galore"} =~ /^[Yy1]/)
{
print TMP "
module load trim_galore
trim_galore --paired $filename\_R1.fastq $filename\_R2.fastq
wait\n"
}


# To run flash to merge the paired end reads into single end reads
if ($param{"flash"} =~ /^[Yy1]/)
{
print TMP "module load flash
flash $filename\_R1_val_1.fq $filename\_R2_val_2.fq -o $filename
mv $filename.extendedFrags.fastq $filename\_ext.fastq\n";
if ($param{"gzip"} eq "Y")
    {
    print TMP "gzip $filename\_R1.fastq\n";
    print TMP "gzip $filename\_R2.fastq\n";
    }
}

# To unzip files if you have previously zipped them
if ($param{"gunzip"} =~ /^[Yy1]/)
{
    if (-e "$filename\_ext.fastq.gz"){print TMP "gunzip $filename\_ext.fastq.gz\n"}
    if (-e "$filename\_ext.fa.gz"){print TMP "gunzip $filename\_ext.fa.gz\n"}
    if (-e "$filename\_ext.psl.gz"){print TMP "gunzip $filename\_ext.psl.gz\n"}
}

# Converts fastq to fa
if ($param{"fq2fa"} =~ /^[Yy1]/)
{
print TMP "sed -n '1~4s/^@/>/p;2~4p' $filename\_ext.fastq > $filename\_ext.fa\n";
}

# Blats sequence
if ($param{"blat"} =~ /^[Yy1]/)
{
#print TMP "module load blat\n$param{'blat_command'} -minScore=20 -minIdentity=5 -maxIntron=10000 -tileSize=11 $param{'oligo_file'} $filename\_ext.fa $filename\_ext.psl\nwait\n"
print TMP "module load ucsctools\nblat -minScore=20 -minIdentity=5 -maxIntron=10000 -tileSize=11 $param{'oligo_file'} $filename\_ext.fa $filename\_ext.psl\nwait\n"
}

# Runs the perl script to split the file into multiple files
if (($param{"MCC_splitter"} =~ /^[Yy1]/) or ($param{"MNase_multi"} =~ /^[Yy1]/))
{
print TMP "perl $param{'master_folder'}/MCC_splitter.pl -f $filename\_ext.fastq -p $filename\_ext.psl -r $param{'reference'}$filename -limit $param{'MNase_Multi_limit'} -all\n"
}

#Gzip to zip the fastq and other files if you have finished using them - it also deletes intermediate files
if ($param{"gzip"} =~ /^[Yy1]/)
{
if (-e "$filename\_R1.fastq"){print TMP "gzip $filename\_R1.fastq\n"} 
if (-e "$filename\_R2.fastq"){print TMP "gzip $filename\_R2.fastq\n"}
if (-e "$filename\_R1_val_1.fq"){print TMP "rm $filename\_R1_val_1.fq"}
if (-e "$filename\_R2_val_2.fq"){print TMP "rm $filename\_R2_val_2.fq"}
if (-e "$filename\_ext.fastq"){print TMP "gzip $filename\_ext.fastq\n"} 
if (-e "$filename\_ext.fa"){print TMP "rm $filename\_ext.fa\n"} 
if (-e "$filename\_ext.psl"){print TMP "gzip $filename\_ext.psl\n\n"} 
}
if (($param{"Align_pipe"} =~ /^[Yy1]/) or ($param{"Pipe2"} =~ /^[Yy1]/))
{
print TMP "perl $param{'master_folder'}/MCC_pipe2.pl -p $path/$param{'reference'}$filename -n $param{'reference'}$filename -i $input_file -q\n";
}

close TMP;

# To qsub the command to sun grid engine.
if ($param{"qsub"} =~ /^[Yy1]/)
{
# For sungrid engine
#system ("qsub -cwd -o qsub.out -e qsub.err -N  $filename < ./tmp.sh");

# For slurm
system ("sbatch -o sb.out -e sb.err tmp.sh");

# To make a log file
system ("cat tmp.sh >> log.txt")
}
}

#This ouputs a 2column hash to a file
sub output_parameters
{
    my ($hashref, $filehandleout_ref) = @_;
    foreach my $value (sort keys %$hashref)
    {
    print $filehandleout_ref "# $value\t".$$hashref{$value}."\n";
    }        
}

sub print_parameters
{
    print "Script parameters\n";
    my ($hashref) = @_;
    foreach my $value (sort keys %$hashref)
    {
    if (defined  $$hashref{$value}){print "$value\t".$$hashref{$value}."\n";}
    else {print "This value wasn't defined properly: $value\n"}
    }        
}


# This is the subroutine that generates the initial config file - it might be worthwhile editing this so that you hard code some of the key parameters for your system
sub config_file
{
if (-e "MCC_pipe_config.txt"){print "MCC_pipe_config.txt already exists. Please delete the previous version and run again to overwrite\n"; exit;}
else
    {
        open CONFIG, ">MCC_pipe_config.txt" or die "couldn't open MCC_pipe_config.txt to output";
        print CONFIG '##########################################################################################
#!bin/bash
# This part of the script allows you to run the pipe by executing the config.txt file
perl '.$pipescript_path.' -q -i ${0##*/}
exit

##########################################################################################
# Specifications that are unlikely to change much - please feel free to change these in the script above and delete from the text below
##########################################################################################
# please note the paths need to be in the format /folder1/folder2 (withoutht a final "/")
# master_folder is the path to the scripts /user/scripts
master_folder='.$param{"master_folder"}.'

# insert path public folder e.g. /public/user and public url for the track hub. The email address is for the UCSC track hub only
public_folder=
public_url=
email=

# Insert path that leads to genome_sizes.txt files for WigToBigWig Please note it expects it to be in the format /path/genome/genome_sizes.txt
bigwig_folder=

# Insert the path for bowtie 2 to use for the reference genome
bowtie_genome_path=

# Insert the build of the genome you are using 
genome=

##########################################################################################
# File specifications
##########################################################################################

# NB the files must be named (filename)_R1.fastq   (filename)_R2.fastq - specify the (filename) - without the _R1.fastq
# Please use fastq not fq
# Please avoid special characters in the filenames such as additional "."s or spaces

file=insert_file1_name_here_in_MCC_config.txt
file=insert_file2_name_here_in_MCC_config.txt
file=insert_file3_name_here_in_MCC_config.txt
file=insert_file4_name_here_in_MCC_config.txt

# Change the reference to the name of the experiment
reference=MCC

# The oligo file is generated by the MCC_oligo_fa.pl script
oligo_file= # Oligos.fa

# Bed file specifying the color of the targets - must be in 8 column bed format with colour in the 8th column (r,g,b format)
colour_file=Color_scheme.bed

##########################################################################################
# Parts of the first pipeline to run
##########################################################################################
trim_galore=N
flash=N
fq2fa=N
blat=N

MCC_splitter=N

clean_wig=N
gunzip=N
gzip=N

##########################################################################################
# Parts of the second pipe parameters
##########################################################################################
Pipe2=N

bowtie=N
sort=N
MCCanal=N
postgzip=N
track_hub=N
samtobam=N
macs2=N';

exit;
    }
}
