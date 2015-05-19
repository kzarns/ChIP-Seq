#!/usr/bin/perl

#This script originated at the NIH
#UND students, faculty, and staff have permission to use it
#The script primarily is used to reorder and optionally reduce bedgraph files
#It may be merged directly into cppmatch in the future
use warnings;
use strict 'vars';

my $strand=1;
my $n=0;

sub print_usage($);

open(INFILE,"$ARGV[0]") || print_usage("i");
open(OUTFILE,">$ARGV[1]") || print_usage("o");

if((scalar @ARGV)<3) {
	$strand=0;
	print "No strand information specfied, strand column will be suppressed.\n";
}

while(<INFILE>) {
	next if /track/ || /^$/;
	/(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/;
	print OUTFILE sprintf("row_%07d\t",$n);
	print OUTFILE "$4\t$1\t$2\t$3";
	if($strand==1) {
		print OUTFILE "\t$ARGV[2]\n";
	}
	else {
		print OUTFILE "\n";
	}
	$n++;
}

sub print_usage($) {
	my $err=shift;
	if($err eq "i") {
		print "Could not open input file \"$ARGV[0]\"\n";
	}
	elsif($err eq "o") {
		print "Could not create output file \"$ARGV[2]\"\n";
	}
	die "Usage:  bedgraph2cppmatch.pl [input file name] [output file name] [strand name (optional)]\n";
}
