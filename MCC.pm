package MCC;
# MCC.pm

use strict;
use warnings;
use File::Basename;
use Cwd 'abs_path';

my $LEVEL = 1;

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
my $pipescript_path=abs_path($0);
my ($pipescript_filename, $pipescript_directory) = fileparse($pipescript_path);
if (-e "MCC_pipe_config.txt"){print "MCC_pipe_config.txt already exists. Please delete the previous version and run again to overwrite\n"; exit;}
else
    {
        open CONFIG, ">MCC_pipe_config.txt" or die "couldn't open MCC_pipe_config.txt to output";
        print CONFIG '##########################################################################################
#!bin/bash
set -e
set -u
set -o pipefail

# This part of the script allows you to run the pipe by executing the config.txt file
perl '.$pipescript_path.' -q -i ${0##*/}
exit

##########################################################################################
# Specifications that are unlikely to change much - please feel free to change these in the script above and delete from the text below
##########################################################################################
# please note the paths need to be in the format /folder1/folder2 (withoutht a final "/")
# master_folder is the path to the scripts /user/scripts
master_folder='.$pipescript_path.'

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

reference=Please_change_the_reference_in_MCC_config.txt

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




sub default_parameters
{
my ($param) =  @_;
# Default parameters stored in %para hash
# Default parameters stored in %para hash
my $pipescript_path=abs_path($0);
my ($pipescript_filename, $pipescript_directory) = fileparse($pipescript_path);
$$param{"master_folder"} = substr($pipescript_directory,0,-1); # The folder where you plan to store the scripts - defaults to the dirctory this script is stored in
$$param{"public_folder"} = ""; # The public folder base directory
$$param{"public_url"} = ""; # The name of your public url e.g. sara.molbiol.ox.ac.uk/public/user
$$param{"bigwig_folder"} = ""; # The path to your bigwig files - it will expect it to be structured /databank/raw/bigwig for the name of the root followed by /genome for the genome followed by genome_size.txt (e.g. if the full path to the file is /genomes/bigwig/mm9/mm9_sizes.txt - specifiy here /genomes/bigwig )
$$param{"blat_command"} = ""; # The command to use BLAT on the server for the Oxford CBRG server use "/package/blat/35/bin/blat"
$$param{"reference"}="MCC"; #This is specified in the config file 
$$param{'oligo_file'}=""; #Path to the oligo file 
$$param{'colour_file'}=""; #Path to a color scheme file if you want one
$$param{"genome"}="mm9";
$$param{"trim_galore"}="N";
$$param{"flash"}="N";
$$param{"gunzip"}="N";
$$param{"gzip"}="N";
$$param{"fq2fa"}="N";
$$param{"blat"}="N";
$$param{"MCC_splitter"}="N";
$$param{"MNase_multi"}="N";
$$param{"Align_pipe"}="N";
$$param{"Pipe2"}="N";
$$param{"qsub"}=0;
$$param{"bowtie"}=0;
$$param{"samtools"}=0; #not normally necessary to run unless trouble shooting 
$$param{"samtobam"}=0; #not normally necessary but can be useful
$$param{"bamtosam"}=0;
$$param{"sort"}=0;
$$param{"MCCanal"}=0;
$$param{"postgzip"}=0;
$$param{"postgunzip"}=0;
$$param{"track_hub"}=0;
$$param{"double_norm"}=0;
$$param{"combine"}=0;
$$param{"clean_wig"}=0;
$$param{"MNase_Multi_limit"}=0;
$$param{"name"}=$$param{"reference"}; #Required to make a unque name for the track hub
$$param{"macs2"}=0;
$$param{"macs2junction"}=0;

$$param{"normalisation_dhs_bed"}='';
}




sub submit
{
my ($filename, $ref) =  @_;
if ($ref eq ""){$ref = "UNDEF"}

if (-e $filename){system ("sbatch -o sb.out -e sb.err -J $ref $filename")}
else {print STDERR "######################## Unable to submit job for $filename #######################\n"}
}
