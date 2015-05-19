#!/usr/bin/perl

#normalize_bedgraph.pl
#
#This script was written by kris zarns
#
#This script takes a bedgraph file and normalizes it in one of several ways
#
#Normalize by number of genes
#Normalize by number of non-zero genes
#Normalize by number of hits
#?

#generate warnings
use warnings;
#generate warnings when undefined variables are referenced
#use strict 'vars';
use strict;
#allow use of threads
use threads;
#use get options package
use Getopt::Std;
use Time::HiRes qw(time);
use POSIX;

#initialize variables
my $in_file_name_1;
my $in_file_name_2;
my $out_dir = "./";
my $type_policy = "up";
my $normalize_policy = "hits";
my $multiplier = 1;

my %options = ();

#extract options from the command line argument
getopts('o:t:n:m:', \%options) || print_usage();

if(exists($ARGV[0])){
	$in_file_name_1 = $ARGV[0];
}else{
	print_usage();
}

if(exists($ARGV[1])){
	$in_file_name_2 = $ARGV[1]
}else{
	print_usage();
}

if(exists($options{'o'})){
	$out_dir = $options{'o'};
}
my $out_file_name = "$out_dir/tmp.bedgraph";

if(exists($options{'t'})){
	$type_policy = $options{'t'};
}

if(exists($options{'n'})){
	$normalize_policy = $options{'n'};
}

if(exists($options{'m'})){
	$multiplier = $options{'m'};
}

write_normalized_file($in_file_name_1, $in_file_name_2, $out_file_name, $type_policy, $normalize_policy, $multiplier);

sub read_bedgraph_line{
	my($line) = @_;

	my $values = [];
	@$values = split("\t", $line);
	#chromasome, start, end, hits
	if(scalar(@$values) != 4){
		return undef;
	}
	return $values;
}

sub read_bedgraph_file{
	my ($file_name) = @_;

	open(DATA, "<", $file_name) || die "Could not open file $file_name";
	my $lines = [];
	@$lines = <DATA>;

	close DATA;
	return $lines;
}

sub normalize_line{
	my ($line, $factor, $type_policy) = @_;
	my $values = read_bedgraph_line($line);

	if(defined $values){
		my $normalized = -1;
		my $hits = $values->[3];
		if($type_policy eq "up"){
			$normalized = int(ceil($hits / $factor));
		}elsif($type_policy eq "down"){
			$normalized = int(floor($hits / $factor));
		}elsif($type_policy eq "float"){
			$normalized = $hits / $factor;
		}else{
			die "Bad type_policy: $type_policy";
		}
		die "normalize failedfor line: $line, result: $normalized, type_policy: $type_policy" unless $normalized >= 0;

		$values->[3] = $normalized;
	}
	return $values;
}

sub normalize_lines{
	my ($lines, $factor, $type_policy) = @_;

	my $norm_lines = {};
	foreach my $line (@$lines){
		my $values = normalize_line($line, $factor, $type_policy);
		if(defined $values){
			push(@{$norm_lines->{$values->[0]}}, $values);
		}
	}
	return $norm_lines;
}

sub determine_factor{
	my ($lines, $normalize_policy) = @_;

	my $factor = 0;
	
	if($normalize_policy eq "genes"){
		my $size = @$lines;
		$factor = $size;
	}elsif($normalize_policy eq "nz_genes"){
		foreach my $line (@$lines){
			my $values = read_bedgraph_line($line);
			if(defined $values){
				my $hits = $values->[3];

				if($hits > 0){
					$factor += 1;
				}
			}
		}
	}elsif($normalize_policy eq "hits"){
		foreach my $line (@$lines){
			my $values = read_bedgraph_line($line);
			if(defined $values){
				my $hits = $values->[3];
				if($hits > 0){
					$factor += $hits;
				}
			}
		}
	}

	die "bad factor or normalize_policy: factor: $factor, normalize_policy: $normalize_policy" unless $factor > 0 and $multiplier > 0;
	print "factor: $factor\n";
	$factor /= $multiplier;
	print "multiplied factor: $factor\n";

	return $factor;
}


sub normalize_file{
	my ($in_file_name, $type_policy, $normalize_policy) = @_;

	my $lines = read_bedgraph_file($in_file_name);
	my $factor = determine_factor($lines, $normalize_policy, $multiplier);
	print "lines: " . scalar(@$lines) . ", factor: $factor, normalize_policy: $normalize_policy, type_policy: $type_policy\n";
	my $norm_lines = normalize_lines($lines, $factor, $type_policy);
	return $norm_lines;
}

sub write_normalized_file{
	my ($in_file_name_1, $in_file_name_2, $out_file_name, $type_policy, $normalize_policy, $multiplier) = @_;

	#don't do any work unless we have an outfile
	open(OUT_FILE, ">", $out_file_name) || die "Could not open $out_file_name for writing";
	
	my $norm_lines_1 = normalize_file($in_file_name_1, $type_policy, $normalize_policy, $multiplier);
	my $norm_lines_2 = normalize_file($in_file_name_2, $type_policy, $normalize_policy, $multiplier);

	#TEST
	print OUT_FILE "test 1: $in_file_name_1\n";
	foreach my $chrom (keys %$norm_lines_1){
		foreach my $values (@{$norm_lines_1->{$chrom}}){
			print OUT_FILE join("\t", @$values) . "\n";
		}
	}
	print OUT_FILE "test 2: $in_file_name_2\n";
	foreach my $chrom (keys %$norm_lines_2){
		foreach my $values (@{$norm_lines_2->{$chrom}}){
			print OUT_FILE join("\t", @$values) . "\n";
		}
	}
	close OUT_FILE;
}

#subroutine called when an error in command line input is detected
sub print_usage{
	die join("\n",
	  	"Usage: bowtie2bedgraph.pl [options] [input file1] [input file2]",
		"\t-o [output_directory]\tspecify where the results should be sent. Default is the current directory.",
		"\t-t [up|down|float]\tthe type policy. Should results be int rounded up, int rounded down, or float? Default is up.",
		"\t-n [genes|nz_genes|hits]\tthe normalize policy. Should we normalize by the number of genes, number of non-zero genes, or number of hits? Default is hits.",
	);
}

1;
