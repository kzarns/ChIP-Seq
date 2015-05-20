#!/usr/bin/perl

#Created by Adam Burkholder, National Institute of Environmental Health Sciences, 2010-2015

use warnings;
use strict 'vars';
use Getopt::Long;
use DBI;

my $first_end;
my $second_end;
my $first_location;
my $second_location;
my $seq;
my @entries;
my $length;
my $first_id;
my $second_id;
my $min=100;
my $max=200;
my $out_type="b";
my $binsize=25;
my $chr;
my @sorted;
my $result;
my %bins;
my $bin;
my $length_file="";
my $ref_genome="";
my %length_table;
my $end;
my $dbh;
my $rh;
my @result;
my $Dflag=0;
my $mit="mit";

sub usage;

$result=GetOptions("o=s" => \$out_type,"min=i" => \$min,"max=i" => \$max,"b=i" => \$binsize,"l=s" => \$length_file,"r=s" => \$ref_genome, "D" => \$Dflag, "M=s" => \$mit);

if(@ARGV<3 || $result==0) {
	if(@ARGV==1) {
		print "Error: An output file prefix must be specified\n";
	}
	if(@ARGV==2) {
		print "Error: miseq or hiseq must be specified\n";
	}
	usage();
}

if($ref_genome ne "") {																					#if the user specified a UCSC genome, look up chromosome lengths via mysql query
	$dbh=DBI->connect("DBI:mysql:database=".$ref_genome.";host=genome-mysql.cse.ucsc.edu","genome");
	die "Error: could not fetch list of chromosome lengths for reference genome \"$ref_genome\"\n" if $dbh->{Active}==0;
	$rh=$dbh->prepare("select chrom, size from chromInfo");
	$rh->execute();						
	while(@result=$rh->fetchrow_array()) {												#store each row of the UCSC-derived list in the length_table hash
		$length_table{$result[0]}=$result[1];
	}
	die "Error: could not fetch list of chromosome lengths for reference genome \"$ref_genome\"\n" if keys(%length_table)==0;
}
elsif($length_file ne "") {																					#if the user specified a chromosome lengths file, store in length_table hash
	open(LENGTHFILE,"$length_file") || die "Error: Cannot open chromosome length file \"$length_file\"\n";
	while(<LENGTHFILE>) {
			/(\S+)\s+(\S+)/;
			$length_table{$1}=$2;
	}
	close(LENGTHFILE);
}

open(INFILE,"$ARGV[0]") || die "Error: Cannot open input file \"$ARGV[0]\"\n";					#extract left-most mapping position, ID, and chromosome for end1 and end2 reads
while(<INFILE>) {
	$first_end=$_;
	($first_id,$chr,$first_location)=($first_end=~/^([^\t]+)\t\S+\t(\S+)\t(\S+)/);
	$second_end=<INFILE>;
	($second_id,$second_location,$seq)=($second_end=~/^([^\t]+)\t\S+\t\S+\t(\S+)\t(\S+)/);
	if($ARGV[2] eq "miseq") {																	#miseq/nextseq instruments include end information, separated by a space from primary read ID
		$first_id=~s/\s\S+$//;
		$second_id=~s/\s\S+$//;
		$first_id=~s/\/\d$//;																	#strip /1 and /2 from read ID
		$second_id=~s/\/\d$//;
	}
	else {
		$first_id=~s/\/\d$//;
		$second_id=~s/\/\d$//;	
	}
	$first_location++;																			#convert form 0-based coordinates to 1-based
	$second_location++;
	$chr="chrM" if $chr=~m/$mit/i && $Dflag==0;													#unless disabled by user, append "chr" to beginning of chromosome and rename mitochondrial chromosome chrM
	$chr="chr".$chr if $chr!~/chr/ && $Dflag==0;
	if($first_id ne $second_id) {																#confirm end1 and end2 read IDs match
		print "Identifiers of current pair do not match, skipping: $first_id, $second_id\n";
		next;
	}
	$length=$second_location-$first_location+length($seq);
	if($length>=$min && $length<=$max) {														#exclude pairs shorter/longer than user-specified limits
		if($out_type ne "g") {																	#if bed output is requested, store fragment start and end in array
			push(@entries,[$chr,$first_location,$first_location+$length-1]);
		}
		if($out_type ne "b") {																	#if bedgraph output is requested, determine genomic bin of fragment center, add single count to bin
			$bin=int(int(((2*$first_location)+$length-1)/2)/$binsize)*$binsize;
			$bin=1 if $bin==0;
			$bins{$chr}={} if !exists($bins{$chr});
			$bins{$chr}->{$bin}++;
		}
	}
}
close(INFILE);

if($out_type ne "g") {																							#write bed file, sorted by chromosome and fragment start
	@sorted=sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @entries;
	open(BEDFILE,">${ARGV[1]}.bed") || die "Error: Could not create output file \"${ARGV[1]}.bed\"\n";
	print BEDFILE "track name=\"$ARGV[1]_bed\" description=\"$ARGV[1]\" visibility=full color=179,27,27\n";
	foreach(@sorted) {
		print BEDFILE "$_->[0]\t$_->[1]\t$_->[2]\n";
	}
	close(BEDFILE);
}

if($out_type ne "b") {																												#write bedgraph file, sorted by chromosome and bin start
	open(GRAPHFILE,">${ARGV[1]}.bedgraph") || die "Error: Could not create output file \"${ARGV[1]}.bedgraph\"\n";
	print GRAPHFILE "track type=bedGraph name=\"$ARGV[1]_bedgraph\" description=\"$ARGV[1]\" visibility=full color=179,27,27\n";
	foreach $chr(sort(keys(%bins))) {
		foreach(sort {$a <=> $b} keys(%{$bins{$chr}})) {
			$end=$_+$binsize-1;
			$end-- if $_==1;
			$end=$length_table{$chr} if ($length_file ne "" && $end>$length_table{$chr});
			print GRAPHFILE "$chr\t$_\t$end\t$bins{$chr}->{$_}\n";
		}
	}
	close(GRAPHFILE);
}

sub usage() {
	die "Usage: extract_fragments.pl [options] [input file] [output prefix] [miseq|hiseq]\nAvailable Options:\n\n\t-o [b,g,a]\t\tSpecify output format: BED (b), BEDGRAPH (g),\n\t\t\t\tor both (a) (default=b).\n\t-min [integer]\tSpecify minimum fragment size (default=100).\n\t-max [integer]\tSpecify maximum fragment size (default=200).\n\t-b [integer]\tSpecify bin size to be used in creation of\n\t\t\t\tBEDGRAPH files (default=25).\n\t-l [string]\t\tprovide list of chromosome names and lengths,\n\t\t\t\tto prevent each final bin from extending beyond\n\t\t\t\tthe chromosome's end\n\t-r [string]\t\trather than specifying a chromosome length list,\n\t\t\t\tspecify a UCSC reference genome identifier\n\t\t\t\t(e.g. mm9),to have chromosome lengths fetched\n\t\t\t\tautomatically\n\t-D\t\t\tdisable appending \"chr\" to beginning of\n\t\t\t\tchromosome names and renaming of mitochondrial\n\t\t\t\tchromosome\n\t-M [string]\t\tstring used to match mitochondrial chromosome\n\t\t\t\tname, default: \"mit\"\n";
}
