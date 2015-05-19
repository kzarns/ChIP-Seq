#!/usr/bin/perl

#cpgmatch.pl version 1.0
#written by Kris Zarns
#This script is available for use by UND students, staff, and faculty
#proper OSS license to be added later

use strict;
use warnings;
use Getopt::Std;
use Time::HiRes qw(time);

my $usage = join("\n",
	"perl cpg_match <cpg file> <db file> <cpg interval offset> <db interval offset>",
	"",
	"This script will combine a db_file with a cpg_islands file.",
	"Combinations are determined by comparing two intervals for overlap.",
	"The cpg interval.  This interval ranges from cpg_start - 'cpg interval offset' to cpg_end + 'cpg interval offset'",
	"The tss interval.  This interval ranges from +/- 'db interval offset' of the gene tss.",
);


die $usage unless $#ARGV eq 3;
my $args = [];
foreach my $key (0 .. $#ARGV){
	push (@$args, $ARGV[$key]);
}
print join("\t", @$args) . "\n";

my $cpg_file_name = $ARGV[0];
my $db_file_name = $ARGV[1];
my $cpg_offset = $ARGV[2];
my $db_offset = $ARGV[3];

#open input file, run print_usage subroutine if file cannot be opened
open(CPGFILE, $cpg_file_name) or die "can't open cpg file";

my $cpgs = {};
my $cpg_header = <CPGFILE>;

#if we track what chromosome the db data is from, we can speed things up with a hash here
while(my $line = <CPGFILE>){
	chomp($line);
	my ($chrom, $start, $end, $cpg_id) = split("\t", $line);
	#my ($junk, $id) = split("_", $cpg_id);

	push(@{$cpgs->{$chrom}}, {"start" => $start, "end" => $end, "id" => $cpg_id});
}

close CPGFILE;

#could bucket tss to improve speed
open(DBFILE, $db_file_name) or die "can't open db file";

my $db_header = <DBFILE>;

print join("\t",
	"chrom",
	"gene_id",
	"gene_symbol",
	"gene_tss",
	"cpg_id",
	"cpg_start",
	"cpg_end",
) . "\n";

while(my $line = <DBFILE>){

	my $source = read_db_file_line($line);

	if(not exists $cpgs->{$source->{"chrom"}}){
		print STDERR "chrom: $source->{'chrom'} not found in db_file\n";
		next;
	}

	for my $cpg (@{$cpgs->{$source->{"chrom"}}}){
		if($source->{"tss"} - $db_offset < $cpg->{"end"} + $cpg_offset and $cpg->{"start"} - $cpg_offset < $source->{"tss"} + $db_offset){
			#this check covers overlapping
			#  the start
			#  the end
			#  tss completely overlapping cpg
			#  cpg completely overlapping tss
			
			print STDOUT join("\t",
				$source->{'chrom'},
				$source->{"id"},
				$source->{"symbol"},
				$source->{"tss"},
				$cpg->{'id'},
				$cpg->{'start'},
				$cpg->{'end'},
			). "\n";
		}
	}
}



#from the db_file
sub read_db_file_line{
	my $line = shift;
	chomp($line);

	my ($id, $symbol_tss, $chr, $remainder) = split("\t", $line, 4);

	$symbol_tss =~ /^(.*):([^:]+)$/;
	my $symbol = $1;
	my $tss = $2;

	return {
		"id" => $id,
		"symbol" => $symbol,
		"tss" => $tss,
		"chrom" => $chr,
	};
}


#this is probably garbage
sub read_cpp_results_line{
	my $line = shift;

	my ($id, $symbol_tss, $hits) = split("\t", $line);

	$symbol_tss =~ /^(.*):([^:]+)$/;
	my $symbol = $1;
	my $tss = $2;

	return {
		"id" => $id,
		"symbol" => $symbol,
		"tss" => $tss,		
	};
}

close DBFILE;
1;
