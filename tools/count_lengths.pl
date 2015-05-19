#!/usr/bin/perl

use warnings;
use strict 'vars';

my $first_end;
my $second_end;
my $first_location;
my $second_location;
my $seq;
my %counts;
my $length;
my @sorted;
my $first_id;
my $second_id;

if(@ARGV<3) {
	if(@ARGV==1) {
		print "Error: An output file name must be specified\n";
	}
	if(@ARGV==2) {
		print "Error: miseq or nisc must be specified\n";
	}
	die "Usage: count_lengths.pl [input file name] [output file name] [miseq|nisc]\n";
}

open(INFILE,"$ARGV[0]") || die "Error: cannot open input file \"$ARGV[0]\"\n";
open(OUTFILE,">$ARGV[1]") || die "Error: cannot create output file \"$ARGV[1]\"\n";

while(<INFILE>) {
	$first_end=$_;
	($first_id,$first_location)=($first_end=~/^([^\t]+)\t\S+\t\S+\t(\S+)/);
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
	if($first_id ne $second_id) {
		print "Identifiers of current pair do not match, skipping: $first_id, $second_id\n";
		next;
	}
	$length=$second_location-$first_location+length($seq);
	$counts{$length}++;
}

@sorted=sort {$a <=> $b} keys(%counts);
foreach(@sorted) {
	print OUTFILE "$_\t$counts{$_}\n";
}

close(INFILE);
close(OUTFILE);
