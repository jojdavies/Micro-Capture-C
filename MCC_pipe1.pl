#!/usr/bin/perl -w

use strict;
use Cwd; use Cwd 'abs_path';
#use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use lib dirname (__FILE__);
use MCC;

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
my $pipescript_path=abs_path($0);
my ($pipescript_filename, $pipescript_directory) = fileparse($pipescript_path);


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

#
my %param;
MCC::default_parameters(\%param); #inserts the default parameters into the %param hash

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

if ($config ==1){MCC::config_file()} #outputs the standard config file if requested and exits

# This loads the config file and puts it into the hash
# Loads the MCC config file into the %param hash and @files array
load_config_file("$path/$input_file", \%param, \@files);


#Makes a subdirectory in the public folder
$param{"public_folder_expt"} = $param{"public_folder"}."/".$param{"reference"};
if (-d $param{"public_folder"}){}
else {mkdir $param{"public_folder"}};

system ("test -e log.txt && rm log.txt");

MCC::print_parameters(\%param); 


foreach my $filename (@files)
{
print "Writing tmp.sh for $filename\n";
open TMP, ">tmp.sh" or die "Couldn't open input file tmp.sh";
print TMP "#!/bin/bash
set -e
set -u
set -o pipefail

# Script started at $hour:$min:$sec $mday/$mon/$year
# Filename $filename\n";
MCC::output_parameters(\%param, \*TMP);

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
    print TMP 'if wait; then echo "'."$filename gunzip completed successfully".'">>log.txt;fi'."\n";
}

# To run trim_galore to cut the reads
if ($param{"trim_galore"} =~ /^[Yy1]/)
{
print TMP "
module load trim_galore
test -e $filename\_R1.fastq && trim_galore --paired $filename\_R1.fastq $filename\_R2.fastq || exit 0\n";
print TMP 'if wait; then echo "trim_galore completed successfully">>log.txt;fi'."\n";
# trim_galore --paired $filename\_R1.fastq $filename\_R2.fastq
}


# To run flash to merge the paired end reads into single end reads
if ($param{"flash"} =~ /^[Yy1]/)
{
print TMP "module load flash
test -e $filename\_R1_val_1.fq && flash $filename\_R1_val_1.fq $filename\_R2_val_2.fq -o $filename || exit 0\n";
print TMP 'if wait; then echo "flash completed">>log.txt;fi'."\n";
print TMP "mv $filename.extendedFrags.fastq $filename\_ext.fastq\n";
if ($param{"gzip"} eq "Y")
    {
    if (-e "$filename\_R1.fastq"){print TMP "gzip $filename\_R1.fastq\n"};
    if (-e "$filename\_R2.fastq"){print TMP "gzip $filename\_R2.fastq\n"};
    }
print TMP 'if wait; then echo "'."$filename flash completed successfully".'">>log.txt;fi'."\n";
}


# Converts fastq to fa
if ($param{"fq2fa"} =~ /^[Yy1]/)
{
print TMP "test -e $filename\_ext.fastq && sed -n '1~4s/^@/>/p;2~4p' $filename\_ext.fastq > $filename\_ext.fa || exit 0\n";
print TMP 'if wait; then echo "'."$filename FASTA completed successfully".'">>log.txt;fi'."\n";
}

# Blats sequence
if ($param{"blat"} =~ /^[Yy1]/)
{
#print TMP "module load blat\n$param{'blat_command'} -minScore=20 -minIdentity=5 -maxIntron=10000 -tileSize=11 $param{'oligo_file'} $filename\_ext.fa $filename\_ext.psl\nwait\n"
print TMP "module load ucsctools\nblat -minScore=20 -minIdentity=5 -maxIntron=10000 -tileSize=11 $param{'oligo_file'} $filename\_ext.fa $filename\_ext.psl\n";
print TMP 'if wait; then echo "BLAT completed successfully">>log.txt;fi'."\n";
}

# Runs the perl script to split the file into multiple files
if (($param{"MCC_splitter"} =~ /^[Yy1]/) or ($param{"MNase_multi"} =~ /^[Yy1]/))
{
print TMP "perl $param{'master_folder'}/MCC_splitter.pl -f $filename\_ext.fastq -p $filename\_ext.psl -r $param{'reference'}$filename -limit $param{'MNase_Multi_limit'} -all\n";
print TMP 'if wait; then echo "'."$filename MCC_splitter completed successfully".'">>log.txt;fi'."\n";
}

#Gzip to zip the fastq and other files if you have finished using them - it also deletes intermediate files
if ($param{"gzip"} =~ /^[Yy1]/)
{
if (-e "$filename\_R1.fastq"){print TMP "gzip $filename\_R1.fastq\n"} 
if (-e "$filename\_R2.fastq"){print TMP "gzip $filename\_R2.fastq\n"}
#if (-e "$filename\_R1_val_1.fq"){print TMP "rm $filename\_R1_val_1.fq"}
#if (-e "$filename\_R2_val_2.fq"){print TMP "rm $filename\_R2_val_2.fq"}
#if (-e "$filename\_ext.fastq"){print TMP "gzip $filename\_ext.fastq\n"} 
#if (-e "$filename\_ext.fa"){print TMP "rm $filename\_ext.fa\n"} 
#if (-e "$filename\_ext.psl"){print TMP "gzip $filename\_ext.psl\n\n"} 
}
if (($param{"Align_pipe"} =~ /^[Yy1]/) or ($param{"Pipe2"} =~ /^[Yy1]/))
{

print TMP "perl $param{'master_folder'}/MCC_pipe2.pl -p $path/$param{'reference'}$filename -n $param{'reference'}$filename -i $path/$input_file -q\n";
print TMP 'if wait; then echo "'."#########################################################################################
$filename Pipe2 started".'">>log.txt;fi'."\n";
}

close TMP;

# To qsub the command to sun grid engine.
if ($param{"qsub"} =~ /^[Yy1]/)
{
MCC::submit("tmp.sh", $filename);

# To make a log file
system ('cat tmp.sh >> tmplog.txt
echo "
#########################################################################################


">>tmplog.txt')
}
}

sub load_config_file
{
my ($config_file, $param, $files) =  @_;

open (FH, $config_file) or die "Couldn't open input file $config_file";

while (my $line = <FH>)
{
chomp $line;
    if ($line =~ /^#/){}
    elsif ($line =~ /^\s*$/){}
    elsif (($line =~ /^perl/) or ($line =~ /^exit/)){}
    elsif (($line =~ /(.*)\s*=\s*(.*)/) or ($line =~ /(.*)\s*=\s*(.*)\s*#.*/))
        {
            if ($1 eq "file"){push @$files, $2}
            else{$$param{$1}=$2}
        }
    else {print "line of input file not read: $line\n"}
}

}
