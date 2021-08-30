#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Pod::Usage;

=head1 NAME

MCC_BLATfa.pl

=head1 SYNOPSIS
 
This script generates a reference 'genome' of the target oligoes for BLAT to map the sequences to. 

It can also take in files in the format
chr:start-stop_NameofTarget
or BED format
chr /t start /t stop /t NameofTarget

It requires a path to a reference genome sequence (same as the one used by bowtie) in FASTA format

Please don't use strange characters in the names - just numbers and letters it may mess things up later on. 

=head1 EXAMPLE

perl MCC_BLATFA.pl -f input_file -g path_to_bowtie_fasta_of_the_genome

=head1 OPTIONS

 -f	Input filename 
 -g	Specify the genome file - this uses the bowtie reference file in fasta format for the genome you want to use
 -l     Target lenght (Default 800)

=head1 AUTHOR

 Written by James Davies 2020

=cut


my $target_length = 800;

my $help=0;
my $man=0;



# The GetOptions from the command line
&GetOptions
(
	"f=s"=>\ my $full_filename,	# -f		Input filename 
	"g=s"=>\ my $genome_file,	# -o		File of the genome in FASTA format - the bowtie index file for example
	"l=i"=>\ $target_length,	# -l		Length of the desired target sequence for BLAT -default 800bp
	'h|help'=>\$help, 		# -h or -help 	Help - prints the manual
	'man'=>\$man,			# -man  	prints the manual
);

pod2usage(1) if $help;
pod2usage(-verbose=>2) if $man;
pod2usage(2) unless ($full_filename); 


# Opens the file and splits the filename from the path - this ouputs the file to the same path as the script
unless ($full_filename =~ /(.*)\.(.*)/) {die"filename does not match fasta/fa format"};
my $file_name=$1;
my $file_path="";
if ($file_name =~ /(.*\/)(\V++)/) {$file_path = $1; $file_name = $2};

my $output_filename = $file_name."_split_$target_length"."bp.fa";

# Uploads the genome into a hash for quick acccess
unless (open(FH, $genome_file)) {print "Cannot open file $genome_file\n"; exit;}
my $chr="undefined";
my %genomehash;
while (my $line= <FH>)
{
if ($line =~ />chr(.*)/)
    {
        $chr = $1; chomp $chr;
        next
    }
else
    {
        chomp $line;
        $genomehash{$chr}.=$line;  
    }
}


open(FHOUT, ">$output_filename") or die "Cannot open file $output_filename $!\n";

my %counters;

open(INFH, "$full_filename") or die "can't open input filename";

my %gene_hash;


while (my $line = <INFH>)
{
    chomp $line;
    $counters{"lines searched"}++;
    my $gene; my $chr; my $start; my $end; my $flag=0;
    if ($line =~ /^chr(.*):(\d++)-(\d++)_(.*)/){ $gene = $4; $chr = $1; $start = $2; $end = $3; $flag=1}
    elsif($line =~ /chr(.*)\t(\d++)\t(\d++)\t(.*)/){ $gene = $4; $chr = $1; $start = $2; $end = $3; $flag=1}
    if ($flag==1)
	{
		
		$gene_hash{$gene}++; #generates a number if there is more than one oligo per gene
		
		if ($end<$start) {$counters{"00 Error end coord before start"}++; next}
		else
		{
		my $oligo_mid_point = $start + ($end - $start)/2;
		my $search_start = int($oligo_mid_point -  ($target_length/2));
        my $search_end = $search_start + $target_length;
        my $sequence;
        my $j = 0;

        $sequence = make_oligo ($chr, $search_start, $search_end, $gene, \*FHOUT);
        print "$gene\t$chr\t$search_start\t$search_end\n";
		}
	}
    else{$counters{"lines not parsed"}++}
    
  $counters{"Lines_searched"}++;
}

sub make_oligo
{
    my ($chr, $start, $end, $name, $FH) = @_;
    my $target_length = $end-$start;
    my $sequence = substr($genomehash{$chr},$start ,$target_length);
    print $FH ">chr$chr:$start-$end\_$name\n$sequence\n";
    return $sequence
}
