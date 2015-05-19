#!/usr/bin/perl

#cpg_enrichment.pl version 2.0
#Written by Kris Zarns

use strict;
use warnings;
use Getopt::Std;
use Time::HiRes qw(time);
use File::Basename;
use Data::Dumper;

$Data::Dumper::Terse = 1;
$Data::Dumper::Indent = 1;

my $usage = join("\n",
	"perl cpg_enrichment.pl <cpg_results file> <cpp_results_total file> <db_file_name>",
	"",
	"This script will compare a cpp_results_total file with a cpg_islands file to determine a correlation between enrichment and cpg islands",
	"The cpp_results_total file will be sorted in order of enrichment.  The corresponding genes will be placed into B buckets.",
	"Genes which have a corresponding cpg island will be counted and the count normalized by the number of genes",
);


die $usage unless $#ARGV eq 2;
my $args = [];
foreach my $key (0 .. $#ARGV){
	push (@$args, $ARGV[$key]);
}
print join("\t", @$args) . "\n";

my $cpg_file_name = $ARGV[0];
my $cpp_file_name = $ARGV[1];
my $db_file_name = $ARGV[2];

#open input file, run print_usage subroutine if file cannot be opened
open(CPGFILE, $cpg_file_name) or die "can't open cpg file";

my $cpgs = {};
my $cpg_header1 = <CPGFILE>;
my $cpg_header2 = <CPGFILE>;

my $known_chroms = {};

#if we track what chromosome the db data is from, we can speed things up with a hash here
while(my $line = <CPGFILE>){
	chomp($line);
	my ($chrom, $gene_id, $junk) = split("\t", $line, 3);
	$known_chroms->{$chrom} = 1;

	$cpgs->{$gene_id} = 1;
}

close CPGFILE;

#print STDERR "cpgs: " . Dumper($cpgs);

#the cpp_results_total file doesn't contain non matching genes, we need to construct a list
#we will remove the cpp genes fromthe list to find the remainder and tack that onto the end
#of the results
open(DBFILE, $db_file_name) or die "can't open db_file";
my $master_genes = {};
my $full_master_genes = {};
while(my $line = <DBFILE>){
	chomp($line);
	my ($gene_id, $symbol_tss, $chr, $start, $end, $strand) = split("\t", $line);
	$master_genes->{$gene_id} = 1;
	my ($gene_symbol, $tss) = split(":", $symbol_tss);
	
	$full_master_genes->{$gene_id} = [$chr, join("\t", $chr, $start, $end, "$gene_id:$gene_symbol", 0, $strand)];
}
#print STDERR "master genes: " . Dumper($master_genes);

my $total_genes = keys %$master_genes;
print STDERR "total_genes: $total_genes\n";
sort_cpp($cpp_file_name);


open(CPPFILE, $cpp_file_name) or die "can't open cpp file";

my $cpp_header = <CPPFILE>;
my $line_num = 0;
my $num_bins = 10;
my $bins = {};
my $bin_size = $total_genes / $num_bins;

#bin our matches
while(my $line = <CPPFILE>){

	my $source = read_cpp_results_line($line);

	my $bin_num = int($line_num / $bin_size);

	$bins->{$bin_num}{"max_enrichment"} = $source->{"hits"} if not exists $bins->{$bin_num}{"max_enrichment"} or $source->{"hits"} > $bins->{$bin_num}{"max_enrichment"};

	push(@{$bins->{$bin_num}{"genes"}}, $source->{"gene_id"});

	$line_num += 1;
	delete $master_genes->{$source->{"gene_id"}};
}

#bin our misses
for my $gene_id (keys %$master_genes){
	my $bin_num = int($line_num / $bin_size);

	push(@{$bins->{$bin_num}{"genes"}}, $gene_id);

	$line_num += 1;
}


my ($cpp_short_name, $dirs, $ext) = fileparse($cpp_file_name, qr/\.[^.]*/);

print STDERR "cpp short: $cpp_short_name\n";
my $cpg_gene_files = {};
open(CPG_STATS, ">", "$cpp_short_name.stats") or die "can't open $cpp_short_name.stats file";

#print STDERR "bins: " . Dumper($bins);

print CPG_STATS "Bin\tpol2\tgenes\tmatches\t%CPG\n";
#print CPG_GENES "bin, chromasome, start, end\n";
for my $bin (0 .. $num_bins -1){
	my $matches = 0;

	$bins->{$bin}{"max_enrichment"} = 0 if not exists $bins->{$bin}{"max_enrichment"};
	for my $gene_id (@{$bins->{$bin}{"genes"}}){
		my $chrom = $full_master_genes->{$gene_id}[0];
		if($bin eq 0){
			unless($cpg_gene_files->{$chrom}){
				open(my $chrom_file, ">", "$cpp_short_name" . "_$chrom.bed") or die "can't open $cpp_short_name" . "_$chrom.bed file";
				$cpg_gene_files->{$chrom} = $chrom_file;
			}

			print {$cpg_gene_files->{$chrom}} "$full_master_genes->{$gene_id}[1]\n";
		}
		$matches++ if $cpgs->{$gene_id };
	}

	print CPG_STATS "$bin\t" . $bins->{$bin}{"max_enrichment"} . "\t" . @{$bins->{$bin}{"genes"}} . "\t$matches\t" . int(10000 * $matches / $bin_size)/100 . "\n";
}


sub sort_cpp{
	my $data_file_name = shift;
	#sort by hits
	`(head -n 1 $data_file_name && tail -n +2 $data_file_name | sort -k 3,3nr -S 10%) > $data_file_name.tmp && mv $data_file_name.tmp $data_file_name`;
}

#this is probably garbage
sub read_cpp_results_line{
	my $line = shift;
	chomp($line);

	my ($id, $symbol_tss, $hits) = split("\t", $line);

	$symbol_tss =~ /^(.*):([^:]+)$/;
	my $symbol = $1;
	my $tss = $2;

	return {
		"gene_id" => $id,
		"hits" => $hits,		
	};
}

foreach my $cpg_gene_file (values %$cpg_gene_files){
	close $cpg_gene_file;
}
close CPPFILE;
1;
