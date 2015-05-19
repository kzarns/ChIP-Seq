#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use File::Basename;

my $usage = join("\n",
	"perl deduplicate.pl args <cpp_result_total_file> <output directory>",
	"choose the highest value of duplicate gene symbols",
	"ARGS:",
	"    -b process both cpp_result and cpp_result_total",
	"    -l <num> limit the output to num genes"
);

use Data::Dumper;
$Data::Dumper::Terse = 1;
$Data::Dumper::Purity = 1;

my %option = ();
getopts('bl:', \%option) || print $usage;

my $both = (exists($option{b}) ? 1 : 0);
my $limit = (exists($option{l}) ? $option{l} : undef);
print "opt l: " . $option{l} . "\n";

print "both: $both, limit: $limit\n";
print "argv: $#ARGV\n";

#prcess cpp_result_total
my $cpp_result_total_file_name = $ARGV[0];

open(CPP_RESULT_TOTAL_FILE, $cpp_result_total_file_name) or die "can't open cpp_result_total file: $cpp_result_total_file_name";
my $cpp_result_total_header = <CPP_RESULT_TOTAL_FILE>;

my $cpp_result_total_master = {};

my $count = 0;
while(my $line = <CPP_RESULT_TOTAL_FILE>){

	my $source = read_cpp_result_total_line($line);
	if(not defined $cpp_result_total_master->{$source->{"symbol"}} or $source->{"hits"} > $cpp_result_total_master->{$source->{"symbol"}}->{"hits"}){
		$cpp_result_total_master->{$source->{"symbol"}} = $source;
		$count++;
	}
	last if defined $limit and $count >= $limit;
}
close CPP_RESULT_TOTAL_FILE;

print STDERR "cpp long: $cpp_result_total_file_name\n";
my ($cpp_short_name, $path, $suffix) = fileparse($cpp_result_total_file_name, qr/\.[^.]*/);
print STDERR "cpp short: $cpp_short_name\n";


my $out_cpp_result_total_name = $ARGV[1] . "/$cpp_short_name" . "_deduplicated" . "$suffix";
#print STDERR $out_cpp_result_total_name;
open(CPP_RESULT_TOTAL_OUT, ">", $out_cpp_result_total_name) or die "can't open out cpp_result_total file: $out_cpp_result_total_name";

print CPP_RESULT_TOTAL_OUT $cpp_result_total_header;
foreach my $gene (values %$cpp_result_total_master){
	print CPP_RESULT_TOTAL_OUT join("\t", $gene->{"gene_id"} , join(":", $gene->{"symbol"}, $gene->{"tss"}), $gene->{"hits"}) . "\n";
}

close(CPP_RESULT_TOTAL_OUT);


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

if($both){
#process cpp_result
	my $cpp_result_file_name = $path . "/" . $cpp_short_name . ".cpp_result";

	open(CPP_RESULT_FILE, $cpp_result_file_name) or die "can't open cpp_result file: $cpp_result_file_name";
	my $cpp_result_header = <CPP_RESULT_FILE>;

	my $cpp_result_master = {};

	$count = 0;
	while(my $line = <CPP_RESULT_FILE>){

		my $source = read_cpp_result_line($line);
		if($source->{"hits"} ne "" and (not defined $cpp_result_master->{$source->{"symbol"}} or $source->{"hits"} > $cpp_result_master->{$source->{"symbol"}}->{"hits"})){
			$cpp_result_master->{$source->{"symbol"}} = $source;
			$count++;
		}
		last if defined $limit and $count >= $limit;
	}
	close CPP_RESULT_FILE;

	my $out_cpp_result_name = $ARGV[1] . "/$cpp_short_name" . "_deduplicated.cpp_result";
	#print STDERR $out_cpp_result_name;
	open(CPP_RESULT_OUT, ">", $out_cpp_result_name) or die "can't open cpp file: $out_cpp_result_name";

	print CPP_RESULT_OUT $cpp_result_header;
	foreach my $gene (values %$cpp_result_master){
		print CPP_RESULT_OUT join("\t",
			$gene->{"gene_id"} ,
			join(":", $gene->{"symbol"}, $gene->{"tss"}),
			$gene->{"chr"},
			$gene->{"start"},
			$gene->{"end"},
			$gene->{"q_id"},
			$gene->{"hits"},
			$gene->{"q_chr"},
			$gene->{"q_start"},
			$gene->{"q_end"},
			$gene->{"q_match"},
		) . "\n";
	}

	close(CPP_RESULT_OUT);
}

sub read_cpp_result_line{
	my $line = shift;
	chomp($line);

	my ($id, $symbol_tss, $chr, $start, $end, $q_id, $hits, $q_chr, $q_start, $q_end, $q_match) = split("\t", $line);

	$symbol_tss =~ /^(.*):([^:]+)$/;
	my $symbol = $1;
	my $tss = $2;

	return {
		"gene_id" => $id,
		"symbol" => $symbol,
		"tss" => $tss,
		"chr" => $chr,
		"start" => $start,
		"end" => $end,
		"q_id" => $q_id,
		"hits" => $hits,
		"q_chr" => $q_chr,
		"q_start" => $q_start,
		"q_end" => $q_end,
		"q_match" => $q_match		
	};
}

1;
