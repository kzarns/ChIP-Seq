#!/usr/bin/perl

#This script was written by Kris Zarns

use strict;
use warnings;
use Getopt::Std;
use File::Basename;

my $usage = join("\n",
	"perl sort_heatmap.pl <heatmap_file> [cpp_total_file]",
	"The cpp_total_file is optional.  If that file is present the heatmap will be sorted by it.  If not, it will be sorted by by the sum of the values on each line"
);

my $file_name = $ARGV[0];
open(my $fh, "<", $file_name) or die "$file_name is missing or can't be opened";
my @lines = <$fh>;
close($fh);

my @sorted;
my $count = 0;
if(@ARGV == 1){
	print "sort by heatmap magnitude";
	@sorted = sort {line_val($b) <=> line_val($a)} @lines[9..$#lines];
}elsif(@ARGV == 2){
	print "sort by cpp";
	my $master = read_cpp_result_total_file($ARGV[1]);
	@sorted = sort {cpp_val($b, $master) <=> cpp_val($a, $master)} @lines[9..$#lines];
}else{
	print "Too many or too few args\n$usage";
	exit 1;
}


my ($short_name, $path, $suffix) = fileparse($file_name, qr/\.[^.]*/);
my $out_name = "$path/$short_name" . "_sorted$suffix";

open(my $fho, ">", $out_name) or die "$out_name can't be opened";
print $fho @lines[0..8];
print $fho @sorted;
close($fh);

sub line_val{
	my ($line) = @_;
	chomp($line);
	my (@vals) = split("\t",$line);
	
	my $total = 0;
	$total += $_ for (@vals[7..$#vals]);
	return $total
}

sub cpp_val{
	my ($line, $master) = @_;
	my (@vals) = split("\t",$line);
	my $out = $master->{$vals[1]};
	if(not defined $master->{$vals[1]}){
		#gene was not in the cpp_result_total file
		#gene_list is more inclusive than db_file which excludes anything within an interval of the first gene it finds
		$out = 0;
	}
	return $out;
}

sub read_cpp_result_total_file{

	my ($cpp_name) = @_;
	my $cpp_result_total_master = {};

	open(my $fhcpp, "<", $cpp_name) or die "$ARGV[1] not readable";
	my $header = <$fhcpp>;
	while(my $line = <$fhcpp>){

		my $source = read_cpp_result_total_line($line);
		if(not defined $cpp_result_total_master->{$source->{"symbol_tss"}} or $source->{"hits"} > $cpp_result_total_master->{$source->{"symbol_tss"}}){
			$cpp_result_total_master->{$source->{"symbol_tss"}} = $source->{"hits"};
		}
	}
	close $fhcpp;
	return $cpp_result_total_master;
}

sub read_cpp_result_total_line{
	my $line = shift;
	chomp($line);

	my ($id, $symbol_tss, $hits) = split("\t", $line);
	return {
		"symbol_tss" => $symbol_tss,
		"hits" => $hits,
	};
}

1;
#
#
#
#                                                 1;

#
#print STDERR "cpp long: $cpp_result_total_file_name\n";
#my ($cpp_short_name, $path, $suffix) = fileparse($cpp_result_total_file_name, qr/\.[^.]*/);
#print STDERR "cpp short: $cpp_short_name\n";
#
#
#my $out_cpp_result_total_name = $ARGV[1] . "/$cpp_short_name" . "_deduplicated" . "$suffix";
##print STDERR $out_cpp_result_total_name;
#open(CPP_RESULT_TOTAL_OUT, ">", $out_cpp_result_total_name) or die "can't open out cpp_result_total file: $out_cpp_result_total_name";
#
#print CPP_RESULT_TOTAL_OUT $cpp_result_total_header;
#foreach my $gene (values %$cpp_result_total_master){
#	print CPP_RESULT_TOTAL_OUT join("\t", $gene->{"gene_id"} , join(":", $gene->{"symbol"}, $gene->{"tss"}), $gene->{"hits"}) . "\n";
#}
#
#close(CPP_RESULT_TOTAL_OUT);
#
#

#
#if($both){
##process cpp_result
#	my $cpp_result_file_name = $path . "/" . $cpp_short_name . ".cpp_result";
#
#	open(CPP_RESULT_FILE, $cpp_result_file_name) or die "can't open cpp_result file: $cpp_result_file_name";
#	my $cpp_result_header = <CPP_RESULT_FILE>;
#
#	my $cpp_result_master = {};
#
#	$count = 0;
#	while(my $line = <CPP_RESULT_FILE>){
#
#		my $source = read_cpp_result_line($line);
#		if($source->{"hits"} ne "" and (not defined $cpp_result_master->{$source->{"symbol"}} or $source->{"hits"} > $cpp_result_master->{$source->{"symbol"}}->{"hits"})){
#			$cpp_result_master->{$source->{"symbol"}} = $source;
#			$count++;
#		}
#		last if defined $limit and $count >= $limit;
#	}
#	close CPP_RESULT_FILE;
#
#	my $out_cpp_result_name = $ARGV[1] . "/$cpp_short_name" . "_deduplicated.cpp_result";
##print STDERR $out_cpp_result_name;
#	open(CPP_RESULT_OUT, ">", $out_cpp_result_name) or die "can't open cpp file: $out_cpp_result_name";
#
#	print CPP_RESULT_OUT $cpp_result_header;
#	foreach my $gene (values %$cpp_result_master){
#		print CPP_RESULT_OUT join("\t",
#			$gene->{"gene_id"} ,
#			join(":", $gene->{"symbol"}, $gene->{"tss"}),
#			$gene->{"chr"},
#			$gene->{"start"},
#			$gene->{"end"},
#			$gene->{"q_id"},
#			$gene->{"hits"},
#			$gene->{"q_chr"},
#			$gene->{"q_start"},
#			$gene->{"q_end"},
#			$gene->{"q_match"},
#		) . "\n";
#	}
#
#	close(CPP_RESULT_OUT);
#}
#
#sub read_cpp_result_line{
#	my $line = shift;
#	chomp($line);
#
#	#my ($id, $symbol_tss, $hits) = split("\t", $line);
#	my ($id, $symbol_tss, $chr, $start, $end, $q_id, $hits, $q_chr, $q_start, $q_end, $q_match) = split("\t", $line);
#
#	$symbol_tss =~ /^(.*):([^:]+)$/;
#	my $symbol = $1;
#	my $tss = $2;
#
#	return {
#		"gene_id" => $id,
#		"symbol" => $symbol,
#		"tss" => $tss,
#		"chr" => $chr,
#		"start" => $start,
#		"end" => $end,
#		"q_id" => $q_id,
#		"hits" => $hits,
#		"q_chr" => $q_chr,
#		"q_start" => $q_start,
#		"q_end" => $q_end,
#		"q_match" => $q_match		
#	};
#}

1;
