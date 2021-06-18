#!/usr/bin/perl -w

use strict;
use Cwd;
#use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use lib dirname (__FILE__);
use MCC;

=head1 NAME
 
MCC_pipe2.pl

=head1 SYNOPSIS

This script runs the second part of the MCC_pipeline
It should be called by the MCC_pipe1.pl
Be careful running this script it will potentially execute hundreds of jobs on your server simultaneously

=head1 EXAMPLE

perl MCC_pipe2.pl

=head1 OPTIONS

It is best to use the MCC_config.txt file to enter the parameters

=head1 AUTHOR

 Written by James Davies 2019

=cut


# Specifications

my $sleep_time =10;

# Default parameters stored in %para hash
my %param;
my $input_file="MCC_pipe_config.txt";
MCC::default_parameters(\%param); #inserts the default parameters into the %param hash

# Strings
# NB. The following strings should be changed to to your public account, server url and bigwig folder.

my $path = getcwd;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(); #defines the start time of the script (printed out into the report file later on)
my $uniq_id = "M$min$sec";
my $help=0;
my $man=0;

# Arrays
my @files;

# The GetOptions from the command line
&GetOptions
(
    "i=s"=>\ $param{"input_file"},
    "p=s"=>\ $path,
	"q"=>\ $param{"qsub"},			# -pf		Your public folder (e.g. /hts/data0/public/username)
	"pf=s"=>\ $param{"public_folder"},			# -pf		Your public folder (e.g. /hts/data0/public/username)
	"pu=s"=>\ $param{"public_url"},				# -pu		Your public url (e.g. sara.molbiol.ox.ac.uk/public/username)
    "genome=s"=>\ $param{"genome"},             # -genome   Genome
    "n=s"=>\ $param{"name"},             # -genome   Genome
    "b"=>\ $param{"bowtie"},                    # -b        Run bowtie
    "sam"=>\ $param{"samtools"},                # -sam      Run samtools
    "m"=>\ $param{"MCCanal"},
    "s"=>\ $param{"sort"},
    "gzip"=>\ $param{"postgzip"},
    "hub"=>\ $param{"track_hub"},
	'h|help'=>\$help,				# -h or -help 	Help - prints the manual
	'man'=>\$man,					# -man  	prints the manual
);

pod2usage(1) if $help;
pod2usage(-verbose=>2) if $man;


# Pulls in the MCC_pipe_config file and puts the data into a hash / array of the files
load_config_file("$input_file", \%param, \@files);

# Makes a subdirectory in the public folder
$param{"public_folder_expt"} = $param{"public_folder"}."/".$param{"reference"};
if (-d $param{"public_folder_expt"}){}
else {mkdir $param{"public_folder_expt"}};

# Opens the directory created by the first part of the pipe and pulls in the file names
opendir DIR, $path or die "cannot open dir $path: $!";
@files= readdir DIR;
closedir DIR;

chdir $path;
my $target = '';
if (exists ($param{"single_target"})){$target = $param{"single_target"}}

MCC::print_parameters(\%param);

# Loops through all the files in the directory and prints the commands into a temporary file
foreach my $file (@files)
{
    #if ($file =~ /(.*).fastq/)
    if ($file =~ /(.*$target.*).fastq/)
    {
my $filename=$1;
print "Analysing file $filename\n";
open TMP, ">align_tmp.sh" or die "Couldn't open input file tmp.sh";
print TMP "#!/bin/bash
# Script started at $hour:$min:$sec $mday/$mon/$year
# Filename $filename\n";

MCC::output_parameters(\%param, \*TMP);

# To compress the files afterwards
if ($param{"postgunzip"} =~ /^[Yy1]/)
{
    if (-e "$filename.fastq.gz"){print TMP "gunzip $filename.fastq.gz\n"}
}


# Runs bowtie2
if ($param{"bowtie"} =~ /^[Yy1]/)
{
print TMP "
module load bowtie2
test -e $filename.fastq && bowtie2 -p 1 -X 1000 -x $param{'bowtie_genome_path'} $filename.fastq -S $filename.sam || exit 0\n";
print TMP 'if wait; then echo "Bowtie on'.$filename.'fastq completed successfully">>log.txt;fi'."\n";
}

# Converts to bam back to sam if it has previously been compressed
if ($param{"bamtosam"} =~ /^[Yy1]/) #Reason for elsif is that you don't want to run both of these options - the first trumps the second.
{
print TMP "samtools view -h $filename.bam > $filename.sam\n";
print TMP 'if wait; then echo "BAM to SAM on'.$filename.' completed successfully">>log.txt;fi'."\n";
}

# Convert from sam to bam visualising and trouble shooting
if ($param{"samtools"} =~ /^[Yy1]/)
{
print TMP "samtools view -S -b -o $filename.bam $filename.sam
samtools sort $filename.bam  -o $filename.sorted.bam
samtools index $filename.sorted.bam";
print TMP 'if wait; then echo "SAM to BAM on'.$filename.' completed">>log.txt;fi'."\n";
print TMP "cp $filename.sorted.bam $param{'public_folder_expt'}
cp $filename.sorted.bam.bai $param{'public_folder_expt'}
echo track type=bam name=\"$filename\" bigDataUrl=$param{'public_url'}$param{'public_folder_expt'}/$filename.sorted.bam >> UCSC.txt\n";
print TMP 'if wait; then echo "BAMs on'.$filename.' copied and sorted completed successfully">>log.txt;fi'."\n";
}
elsif ($param{"samtobam"} =~ /^[Yy1]/) #Reason for elsif is that you don't want to run both of these options - the first trumps the second.
{
print TMP "samtools view -S -b -o $filename.bam $filename.sam
samtools sort $filename.bam  -o $filename.sorted.bam\n";
print TMP 'if wait; then echo "SAM to BAM on'.$filename.' completed successfully">>log.txt;fi'."\n";
}

# Sorts the sam file - this is important for the next step in the analysis
if ($param{"sort"} =~ /^[Yy1]/)
{
print TMP "module load samtools
test -e $filename.sam && samtools sort -n -o $filename\_sort.sam $filename.sam || exit 0
mv $filename\_sort.sam $filename.sam\n";
print TMP 'if wait; then echo "Sort on'.$filename.' completed successfully">>log.txt;fi'."\n";
}

# Runs the MCC_analyser.pl script which is the main script which pulls out the junctions
if ($param{"MCCanal"} =~ /^[Yy1]/)
{
print TMP "test -e $filename.sam && perl $param{'master_folder'}/MCC_analyser.pl -f $filename.sam -pf $param{'public_folder_expt'} -bf $param{'bigwig_folder'} -genome $param{'genome'} -o $param{'oligo_file'} || exit 0\n";
print TMP 'if wait; then echo "MCC_analyser.pl on'.$filename.' completed successfully">>log.txt;fi'."\n";
}

# To compress the files afterwards
if ($param{"postgzip"} =~ /^[Yy1]/)
{
    if (-e "$filename.fastq"){print TMP "gzip $filename.fastq\n"}
    if ((-e "$filename.sam") and (-e "$filename.bam")){print TMP "rm $filename.sam\n"}
    elsif (-e "$filename.sam"){print TMP "samtools view -S -b -o $filename.bam $filename.sam\n"}
    print TMP 'if wait; then echo "Gzip / samtobam on '.$filename.' completed successfully">>log.txt;fi'."\n";
}

# To analyse with macs2
if ($param{"macs2"} =~ /^[Yy1]/)
{
if (-d "$path/macs2"){}
else {mkdir "$path/macs2"};

print TMP "module load macs2
test -e $filename.bam && macs2 callpeak -t $filename.bam -f BAM -g1.87e9 --nomodel --extsize 100 -n  $filename --outdir $path/macs2 || exit 0\n";
print TMP 'if wait; then echo "MACS2 on'.$filename.' completed successfully">>log.txt;fi'."\n";
}

# To analyse with macs2
if ($param{"macs2junction"} =~ /^[Yy1]/)
{
if (-d "$path/macs2_junction"){}
else {mkdir "$path/macs2_junction"};

print TMP "module load macs2
test -e $filename\_junction.bam && macs2 callpeak -t $filename\_junction.bam -f BAM -g1.87e9 --nomodel --extsize 100 -n  $filename\_junction --outdir $path/macs2 || exit 0\n";
print TMP 'if wait; then echo "MACS2 junction on'.$filename.' completed successfully">>log.txt;fi'."\n";
}

if ($param{"double_norm"} =~ /^[Yy1]/)
{
    if ((-e "$filename\_de_norm_rep.wig") and (exists($param{"normalisation_dhs_bed"}))){print TMP "perl /t1-data/project/fgenomics/jdavies/00Scripts/MCC_normaliser.pl $filename\_de_norm_rep.wig $param{'normalisation_dhs_bed'}\n"}
}


close TMP;

# To submit to the Sun Grid Engine
if ($param{"qsub"} =~ /^[Yy1]/)
{
# Sun grid engine
#system ("qsub -cwd -o qsub.out -e qsub.err -N  $uniq_id < ./align_tmp.sh");

# For Slurm
system ("chmod 755 align_tmp.sh
sbatch -o sb.out -e sb.err -J $uniq_id align_tmp.sh");

system ('cat align_tmp.sh >> tmplog.txt
        echo "
#########################################################################################



">>tmplog.txt')
}    

    
    }
}

if ($param{"track_hub"} =~ /^[Yy1]/)
{
# Waits until the qstat commands have finished
my $flag=0;
until($flag==1)
    {
    #my $qstat = `qstat -u jdavies`;
    
    my $qstat = `squeue -u jdavies`;
    #print "$uniq_id\t $qstat\n";
    if ($qstat !~ /.*$uniq_id*/){$flag=1; last}
    else {sleep $sleep_time;}
    }

# To make the track hub
print "Making track hub....\n";
opendir DIR, $path or die "cannot open dir $path: $!";
@files= readdir DIR;
closedir DIR;


# Allows the parameters to be specified in the input file
my %colours;
if (-e $param{"colour_file"})
{
    open FH, $param{"colour_file"} or die "Couldn't open input file".$param{"colour_file"};

    while (my $line = <FH>)
    {
    chomp $line;
    my ($chr, $start, $end, $name, $score, $dot, $start2, $end2, $colour) = split /\t/, $line; 
        if ($colour =~/\d++,\d++,\d++/){$colours{$name}=$colour}
        else {print "colour file error $line\n"}
    }
}

#print Dumper (\%colours); exit;

open UCSC, ">>UCSC.txt" or die "Couldn't open input file UCSC.txt";
print UCSC "$param{'public_url'}$param{'public_folder_expt'}/$param{'name'}_hub.txt\n";
close UCSC;

open HUB, ">$param{'public_folder_expt'}/$param{'name'}_hub.txt" or die "couldn't make hub.txt file";
print HUB "hub $param{'name'}
shortLabel $param{'name'}
longLabel $param{'name'}
genomesFile $param{'public_url'}$param{'public_folder_expt'}/$param{'name'}_genomes.txt
email $param{'email'}";
close HUB;

open GENOMES, ">$param{'public_folder_expt'}/$param{'name'}_genomes.txt" or die "Couldn't open input file genomes.txt";
print GENOMES "genome $param{'genome'}
trackDb $param{'public_url'}$param{'public_folder_expt'}/$param{'name'}_tracks.txt";
close GENOMES;

open TRACKS, ">$param{'public_folder_expt'}/$param{'name'}_tracks.txt" or die "Couldn't open input file tracks.txt";


my $colour;
my $i=0;

foreach my $file (@files)
{
    if ($file =~ /(.*).wig/)
    {
    my $file_short=$1;
    my $filesize = -s  $file;  # checks the filesize
    if ($filesize >1000)
    {
$colour = colour_find($file_short, \%colours, $i);
print TRACKS "track $file_short
type bigWig 
longLabel $file_short
shortLabel $file_short
bigDataUrl $param{'public_url'}$param{'public_folder_expt'}/$file_short.bw
visibility hide
autoScale off
viewLimits 0:100
priority 200
color $colour
alwaysZero on

"
    }
    }
}
}

# This searches the %colour_hash to see if a colour has been specified for the target. If it hasn't it returns a colour from the array
sub colour_find
{
    my ($file, $colour_hash, $i) = @_;
    my @colours =("19,0,24", "194,0,11", "52,153,128", "255,128,0", "130,2,126", "0,160,198", "83,42,133", "135,54,0", "11,83,69", "27,38,49");
    foreach my $target (keys %$colour_hash)
    {
        #print "$file\t$target\n"; 
        if ($file =~ /.*$target.*/){return $$colour_hash{$target}}
    }
    return $colours[$i]
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
