#!/bin/perl
use strict;
use warnings;

use Data::Dumper;
use Getopt::Std;
use File::Basename;

#Written by Kris Zarns
#
#This script takes a tss file with lines of form 
#name, chrom, strand, txstart, txend, geneSymbol
#forms a interval around the transcription start site according to strandedness
#reduces duplicate symbols occuring with overlapping intervals
#reorders to name, genSymbol, crhom, txstart, txend, strand for the cppmatch utility
#
#The script also provides methods for limiting output based on a list of genes or cpp results
#The script also provides an option to generate gene_list files for use in draw_heatmap
#gene list files can also be used by match_genes to generate a fasta file
my %option;
getopts('i:o:f:r:s:c', \%option) || print_usage();

$Data::Dumper::Indent = 1;
$Data::Dumper::Terse = 1;

#represents our expected starting order
my $enum = {
	gene_id      => 0,
	chrom        => 1,
	strand       => 2,
	start        => 3,
	end          => 4,
	gene_symbol  => 5,
};


my $file_name = $option{i};
my $forward_interval = $option{f};
my $reverse_interval = $option{r};
my $out_dir = $option{o};
my $sub_name = $option{s};

my $subset = {};
my $use_subset = 0;
#read in a subset file to use for limiting output to a subset of gene symbols
if($sub_name){
	$use_subset = 1;
	open(my $sub_fh, "<", $sub_name) or die "can't open subset file: $sub_name";

	if($option{c}){
		my $header = readline($sub_fh);
		while(my $line = readline($sub_fh)){
			my $source = read_cpp_result_total_line($line);
			$subset->{$source->{"symbol"}} = 1;
		}
	}else{
		while(my $line = readline($sub_fh)){
			chomp($line);
			$subset->{$line} = 1;
		}
	}
	close $sub_fh;
}

sub read_cpp_result_total_line{
	my $line = shift;
	chomp($line);

	my ($id, $symbol_tss, $hits) = split("\t", $line);

	$symbol_tss =~ /^(.*):([^:]+)$/;
	my $symbol = $1;
	my $tss = $2;

	return {
		"gene_id" => $id,
		"symbol" => $symbol,
		"tss" => $tss,
		"hits" => $hits,		
	};
}

my $usage = join("\n",
	"USAGE: perl tss2db_file.pl -c -i tss_file_name -f interval_forward -r interval_reverse -s [subset_file_name] -o output_dir",
	"    -c using subset file and file is of form cpp_result_total",
	"    -i the tss input file",
	"    -f interval forward",
	"    -r interval reverse",
	"    -s optional subset file default form is one gene symbol per line",
	"    -o output dir",
);

die $usage unless defined $file_name and defined $out_dir;

#switch between using user defined intervals or strictly start/end sites from the tss file
#if no interval is given this will output a gene_listand the ordering of start and end will not be swapped according to strand
my $use_intervals = 1;
unless(defined $forward_interval and defined $reverse_interval){
	$use_intervals = 0;
}

my ($short_name, $path, $suffix) = fileparse($file_name, qr/\.[^.]*/);

#generate output file name
my $out_name = "$out_dir/" . $short_name;
$out_name .= "_+$forward_interval" . "_-$reverse_interval" if $use_intervals;

if($sub_name){
	my ($sub_short_name, $sub_path, $sub_suffix) = fileparse($sub_name, qr/\.[^.]*/);
	$out_name = $out_name . "_$sub_short_name";
}
if($use_intervals){
	$out_name .= ".db_file";
}else{
	$out_name .= ".gene_list";
}

#read the tss file
open(my $in_fh, "<", $file_name) or die "can't open $file_name";

my $ordered_results = [];
my $seen = {};

while(my $line = readline($in_fh)){
	#print "base: $line";
	my $values = [];
	@$values = split(/\s+/, $line);

	#copy the header forward
	if($line =~ /^#/){
		push(@$ordered_results, $values);
		next;
	}

	#FLIP as needed to account for negative strand start coming after end
	#flip when using intervals
	if($values->[$enum->{strand}] eq '-' and $use_intervals){
		$values->[$enum->{start}] = $values->[$enum->{end}];
	}

	if($use_subset and not defined $subset->{$values->[$enum->{gene_symbol}]}){
		next;
	}

	#don't record nearby duplicates
	my $found = 0;
	if($seen->{$values->[$enum->{gene_symbol}]}){

		#don't remember values we've seen already if they overlap the start site of known values 
		foreach my $set (@{$seen->{$values->[$enum->{gene_symbol}]}}){
			if($set->[$enum->{start}] <= $values->[$enum->{start}] and $set->[$enum->{end}] >= $values->[$enum->{start}]){
				$found  = 1;
			}
		}
	}
	if(!$found){

		my $initial_gene_symbol = $values->[$enum->{gene_symbol}];
		if($values->[$enum->{strand}] eq '+'){
			#concatenate the start to the gene_symbol
			$values->[$enum->{gene_symbol}] = join(":", $values->[$enum->{gene_symbol}], $values->[$enum->{start}]);
			if($use_intervals){
				$values->[$enum->{end}] = $values->[$enum->{start}] + $forward_interval;
				$values->[$enum->{start}] = $values->[$enum->{start}] - $reverse_interval;
			}
		}else{
			#minus strand, concatenate the end to the gene_symbol
			$values->[$enum->{gene_symbol}] = join(":", $values->[$enum->{gene_symbol}], $values->[$enum->{end}]);
			if($use_intervals){
				$values->[$enum->{end}] = $values->[$enum->{start}] + $reverse_interval;
				$values->[$enum->{start}] = $values->[$enum->{start}] - $forward_interval;
			}
		}
		#remember what we've seen
		push(@{$seen->{$initial_gene_symbol}}, $values);
		#remember what order it was in
		push(@$ordered_results, $values);
	}
}

close $in_fh;

my $strand_map = {"+" => "plus", "-" => "minus"};

open(my $out_fh, ">", $out_name) or die "can't open $out_name";
foreach my $values (@$ordered_results){
	unless($use_intervals){
		$values->[$enum->{strand}] = $strand_map->{$values->[$enum->{strand}]};
	}
	#reorders to name, genSymbol, crhom, txstart, txend, strand for the cppmatch utility
	print $out_fh join("\t",
		$values->[$enum->{gene_id}],
		$values->[$enum->{gene_symbol}],
		$values->[$enum->{chrom}],
		$values->[$enum->{start}],
		$values->[$enum->{end}],
		$values->[$enum->{strand}],
	) . "\n"
}
close $out_fh;
#foreach my $key (keys %$seen){
#	foreach my $set (@{$seen->{$key}}){
#		print "$set->[$enum->{gene_id}]\t$set->[$enum->{strand}]\t$set->[$enum->{gene_symbol}]\t" . (($set->[$enum->{chrom}] eq "-" and $set->[$enum->{end}] > $set->[$enum->{start}]) ? "$set->[$enum->{end}]\t$set->[$enum->{start}]" : "$set->[$enum->{start}]\t$set->[$enum->{end}]") . "\t$set->[$enum->{chrom}]\n";
#	}
#}
