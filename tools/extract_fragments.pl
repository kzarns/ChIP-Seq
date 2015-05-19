#!/usr/bin/perl

use warnings;
use strict 'vars';
use Getopt::Long;
#use DBI;

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
#my $dbh;
my $rh;
my @result;

sub usage;

$result=GetOptions("o=s" => \$out_type,"min=i" => \$min,"max=i" => \$max,"b=i" => \$binsize,"l=s" => \$length_file,"r=s" => \$ref_genome);

if(@ARGV<3 || $result==0) {
	if(@ARGV==1) {
		print "Error: An output file prefix must be specified\n";
	}
	if(@ARGV==2) {
		print "Error: miseq or nisc must be specified\n";
	}
	usage();
}

#if($ref_genome ne "") {
#	$dbh=DBI->connect("DBI:mysql:database=".$ref_genome.";host=genome-mysql.cse.ucsc.edu","genome");	#connect the the user-specified UCSC database
#	die "Error: could not fetch list of chromosome lengths for reference genome \"$ref_genome\"\n" if $dbh->{Active}==0;
#	$rh=$dbh->prepare("select chrom, size from chromInfo");								#prepare and execute a request for a list of chromosomes and their lengths
#	$rh->execute();						
#	while(@result=$rh->fetchrow_array())												#store each row of the list in the length_table hash
#		{
#		$length_table{$result[0]}=$result[1];
#		}
#	die "Error: could not fetch list of chromosome lengths for reference genome \"$ref_genome\"\n" if keys(%length_table)==0;
#	}
#elsif($length_file ne "") {
	open(LENGTHFILE,"$length_file") || die "Error: Cannot open chromosome length file \"$length_file\"\n";
	while(<LENGTHFILE>) {
			/(\S+)\s+(\S+)/;
			$length_table{$1}=$2;
	}
	close(LENGTHFILE);
#}


open(INFILE,"$ARGV[0]") || die "Error: Cannot open input file \"$ARGV[0]\"\n";

while(<INFILE>) {
	
	$first_end=$_;
	($first_id,$chr,$first_location)=($first_end=~/^([^\t]+)\t\S+\t(\S+)\t(\S+)/);
	$second_end=<INFILE>;
	($second_id,$second_location,$seq)=($second_end=~/^([^\t]+)\t\S+\t\S+\t(\S+)\t(\S+)/);
	if($ARGV[2] eq "miseq") {
		$first_id=~s/\s\S+$//;
		$second_id=~s/\s\S+$//;
		$first_id=~s/\/\d$//;
		$second_id=~s/\/\d$//;
	}
	else {
		$first_id=~s/\/\d$//;
		$second_id=~s/\/\d$//;	
	}
	$first_location++;
	$second_location++;
	$chr="chr".$chr if($chr!~/chr/);
	if($first_id ne $second_id) {
		print "Identifiers of current pair do not match, skipping: $first_id, $second_id\n";
		next;
	}
	$length=$second_location-$first_location+length($seq);
	if($length>=$min && $length<=$max) {
		if($out_type ne "g") {
			push(@entries,[$chr,$first_location,$first_location+$length-1]);
		}
		if($out_type ne "b") {
			$bin=int(int(((2*$first_location)+$length-1)/2)/$binsize)*$binsize;
			$bin=1 if $bin==0;
			$bins{$chr}={} if !exists($bins{$chr});
			$bins{$chr}->{$bin}++;
		}
	}
}
close(INFILE);

if($out_type ne "g") {
	@sorted=sort {$a->[0] cmp $b->[0] || $a->[1] <=> $b->[1]} @entries;
	open(BEDFILE,">${ARGV[1]}.bed") || die "Error: Could not create output file \"${ARGV[1]}.bed\"\n";
	print BEDFILE "track name=\"$ARGV[1]_bed\" description=\"$ARGV[1]\" visibility=full color=179,27,27\n";
	foreach(@sorted) {
		print BEDFILE "$_->[0]\t$_->[1]\t$_->[2]\n";
	}
	close(BEDFILE);
}

if($out_type ne "b") {
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
	die "Usage: extract_fragments.pl [options] [input file name] [output file prefix] [miseq|nisc]\nAvailable Options:\n\n\t-o [b,g,a]\t\tSpecify output format: BED (b), BEDGRAPH (g),\n\t\t\t\tor both (a) (default=b).\n\t-min [numeric value]\tSpecify minimum fragment size (default=100).\n\t-max [numeric value]\tSpecify maximum fragment size (default=200).\n\t-b [numeric value]\tSpecify bin size to be used in creation of\n\t\t\t\tBEDGRAPH files (default=25).\n\t-l [file name]\t\tprovide list of chromosome names and lengths,\n\t\t\t\tto prevent each final bin from extending beyond\n\t\t\t\tthe chromosome's end\n\t-r [genome]\t\trather than specifying a chromosome length list,\n\t\t\t\tspecify a UCSC reference genome identifier\n\t\t\t\t(e.g. mm9),to have chromosome lengths fetched\n\t\t\t\tautomatically\n";
}
