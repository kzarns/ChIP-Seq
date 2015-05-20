#!/usr/bin/perl

#Created by Adam Burkholder, National Institute of Environmental Health Sciences, 2010-2015

use warnings;								
use strict 'vars';							
use threads;								
use Getopt::Std;							
use DBI;									

my @cols;															
my (%plus_table,%minus_table,%plus_bin_table,%minus_bin_table,%merge_table,%binned_merge_table,%option,%length_table);
my ($key,$plus_thread,$minus_thread,$plus_bin_thread,$minus_bin_thread,$merge_thread,$binned_merge_thread,$temp);
my ($normalized,$binned,$merged,$binned_merged,$bin_size,$shift_val,$trim,$list,$shifted,$swap);
my $norm_value=0;
my $name;
my $dbh;
my $rh;
my @result;
my $mit;
my $Dflag;

getopts('o:b:s:t:l:r:xCM:',\%option) || print_usage();				
									
if(exists($option{o})) {							
	if($option{o}=~/[^nbmv]/) {
		print_usage("o");					
	}
	$normalized=($option{o}=~/n/i) || 0;				
	$binned=($option{o}=~/b/i) || 0;
	$merged=($option{o}=~/m/i) || 0;
	$binned_merged=($option{o}=~/v/i) || 0;
}

print_usage("bz") if (exists($option{b}) && $option{b}<1);							#bin size must be 1 or greater
print_usage("b") if ((!$binned && !$binned_merged) && exists($option{b}));			#check for contradictions in specified parameters
print_usage("s") if ((!$merged && !$binned_merged) && exists($option{s}));

$bin_size=$option{b} || 25;								
$shift_val=$option{s} || 75;							
$trim=$option{t} || 0;									
$swap=$option{x} || 0;
$Dflag=$option{D} || 0;
$mit=$option{M}	|| "mit";				

if(exists($option{r})) {								
	$list=$option{r};
	$dbh=DBI->connect("DBI:mysql:database=".$list.";host=genome-mysql.cse.ucsc.edu","genome");	#if user specifies UCSC genome, look up chromosome lengths via mysql query
	print_usage("db") if $dbh->{Active}==0;
	$rh=$dbh->prepare("select chrom, size from chromInfo");								
	$rh->execute();						
	while(@result=$rh->fetchrow_array()) {												#store each row of the list in the length_table hash
		$length_table{$result[0]}=$result[1];
	}
	print_usage("db") if keys(%length_table)==0;
}
elsif(exists($option{l})) {
	$list=$option{l};
	open(LISTFILE,$list) || print_usage("lf");
	while(<LISTFILE>) {							#if user specifies chromosome length list, store names/values in length_table hash
		/(\S+)\s+(\S+)/;
		$length_table{$1}=$2;
	}
	close(LISTFILE);
}
else {
	$list="";									#if list file name is not specified, assign blank value
}									

print_usage("") if scalar(@ARGV)==0;
					
open(INFILE,"$ARGV[0]") || print_usage("f");		#open input file, run print_usage subroutine if file cannot be opened

$name=$ARGV[1] || print_usage("t");			#check for output file prefix in argument
$name=~s/\/?(\S*\/)*//;						#strip path from output file prefix

while(<INFILE>) {
	$norm_value++;										#increment overall hit count
	@cols=/^[^\t]+\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t/;		#extract the second through fourth fields of the current line
	$cols[1]="chrM" if $cols[1]=~m/$mit/i && $Dflag==0;				#if user allows chromosome name fixing, rename mitochrondrial chromsome, append chr where absent
	if($cols[1]!~m/chr/i && $Dflag==0) {
		$cols[1]="chr".$cols[1];
	}
	if($cols[0] eq '+') {
		$cols[2]+=(1-$trim);					#convert from 0-based to 1-based coordinates and adjust for trimming
		$key="$cols[1]\t$cols[2]";				#concatenate the chromosome and 5' mapping position, and use it as a hash key
		$plus_table{$key}++;
		if($binned) {	
			$temp=1 if ($temp=int($cols[2]/$bin_size)*$bin_size)<1;		#determine bin, use 1 if calculated value is zero
			$key="$cols[1]\t$temp";										#use chromosome and bin start as hash key
			$plus_bin_table{$key}++;
		}
		$shifted=$cols[2]+$shift_val;									#determine shifted location for merged output
	}
	elsif($cols[0] eq '-') {						
		$cols[2]+=(length($cols[3])+$trim); 			#determine 5' mapping position in 1-based coordinates, adjust for trimming
		$key="$cols[1]\t$cols[2]";						#perform same tasks as above
		$minus_table{$key}++;
		if($binned) {
			$temp=1 if ($temp=int($cols[2]/$bin_size)*$bin_size)<1;
			$key="$cols[1]\t$temp";
			$minus_bin_table{$key}++;
		}
		$shifted=$cols[2]-$shift_val;
	}
	if($merged) {
		$temp=$shifted;
		$temp=1 if $temp<1;
		if(exists($length_table{$cols[1]})) {					
			$temp=$length_table{$cols[1]} if $temp>$length_table{$cols[1]};				
		}
		elsif($list ne "") {																
			die "Error: chromosome $cols[1] could not be found in the chromosome lengths list\n";
		}
		$key="$cols[1]\t$temp";							#use chromosome and shifted 5' mapping position as hash key
		$merge_table{$key}++;							#increment count, for both forward and reverse strand hits
	}
	if($binned_merged) {
		$temp=$shifted;
		$temp=1 if ($temp=int($temp/$bin_size)*$bin_size)<1;
		if(exists($length_table{$cols[1]})) {						
			$temp=int($length_table{$cols[1]}/$bin_size)*$bin_size if $temp>$length_table{$cols[1]};		#determine bin			
		}
		elsif($list ne "") {																
			die "Error: chromosome $cols[1] could not be found in the chromosome lengths list\n";
		}
		$key="$cols[1]\t$temp";					#use chromosome and bin start as hash key				
		$binned_merge_table{$key}++;			#increase bin hit count, do this for both forward and reverse strand hits
	}
}
close(INFILE);								

$norm_value/=1000000;							#divide overall hit count by one million, hit values will be normalized to millions of overall hits

if($swap==0) {
	$plus_thread=threads->create('sort_tables',\%plus_table,"forward");						#create separate threads for sorting the forward and reverse hits
	$minus_thread=threads->create('sort_tables',\%minus_table,"reverse");					#pass a reference to the appropriate hash along with a string to append to the output file
	$plus_bin_thread=threads->create('sort_tables_alt',\%plus_bin_table,"forward_binned") if $binned;
	$minus_bin_thread=threads->create('sort_tables_alt',\%minus_bin_table,"reverse_binned") if $binned;
}
else {																				#swap strands if requested by user
	$plus_thread=threads->create('sort_tables',\%minus_table,"forward");			
	$minus_thread=threads->create('sort_tables',\%plus_table,"reverse");			
	$plus_bin_thread=threads->create('sort_tables_alt',\%minus_bin_table,"forward_binned") if $binned;
	$minus_bin_thread=threads->create('sort_tables_alt',\%plus_bin_table,"reverse_binned") if $binned;
}
$merge_thread=threads->create('sort_tables_alt',\%merge_table,"merged") if $merged;
$binned_merge_thread=threads->create('sort_tables_alt',\%binned_merge_table,"binned_merged") if $binned_merged;	

$plus_thread->join();							#join threads
$minus_thread->join();
$plus_bin_thread->join() if $binned;
$minus_bin_thread->join() if $binned;
$merge_thread->join() if $merged;
$binned_merge_thread->join() if $binned_merged;

sub sort_tables {								#subroutine used by plus and minus threads for sorting
	my $ref=shift;										
	my $strand=shift;
	my @result;
	my $value;
	open(OUTFILE,">$ARGV[1]_${strand}\.bedgraph") || die "Could not create output file: $ARGV[1]_${strand}\.bedgraph\n";			
	@result=sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] }	map { [/(^.*)\t(.*$)/] } keys(%$ref);						#retrieve the hash keys, split them into chromosome and location components, then sort first by chromosome name alphabetically, next by location numerically
	print OUTFILE "track type=bedGraph name=${name}_$strand description=${name}_$strand visibility=full color=179,27,27 altColor=179,27,27 priority=20\n\n";
	if($normalized) {																						
		open(NORMFILE,">$ARGV[1]_norm1m_${strand}\.bedgraph") || die "Could not create output file: $ARGV[1]_norm1m_${strand}\.bedgraph\n";											#create normalized output file if requested by user
		print NORMFILE "track type=bedGraph name=${name}_norm1m_${strand} description=${name}_norm1m_${strand} visibility=full color=179,27,27 altColor=179,27,27 priority=20\n\n";	
		foreach (@result) {													 
			$value=$ref->{"$_->[0]"."\t"."$_->[1]"};				#for each key in the result array, retrieve hit value by rejoining the key components
			print OUTFILE "$_->[0]\t$_->[1]\t$_->[1]\t$value\n";	
			$value/=$norm_value;									#divide hit value by normalization constant
			print NORMFILE "$_->[0]\t$_->[1]\t$_->[1]\t$value\n";	
		}
		close(NORMFILE);						
	}
	else {									#if normalized output is not requested, only write to the standard output file
		foreach (@result) {									 
			$value=$ref->{"$_->[0]"."\t"."$_->[1]"};				#for each key in the result array, retrieve hit value by rejoining the key components
			print OUTFILE "$_->[0]\t$_->[1]\t$_->[1]\t$value\n";	
		}
	}
	close(OUTFILE);
	return;
}
	
sub sort_tables_alt {							#subroutine used by bin and merge threads for sorting
	my $ref=shift;							
	my $strand=shift;
	my @result;
	my ($value,$end);
	open(OUTFILE,">$ARGV[1]_${strand}\.bedgraph") || die "Could not create output file: $ARGV[1]_${strand}\.bedgraph\n";
	@result=sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } map { [/(^.*)\t(.*$)/] } keys(%$ref);
	print OUTFILE "track type=bedGraph name=${name}_$strand description=${name}_$strand visibility=full color=179,27,27 altColor=179,27,27 priority=20\n\n";
	foreach (@result) {													 
		$value=$ref->{"$_->[0]"."\t"."$_->[1]"};
		if($strand=~/bin/) {														#if outfile is binned
			$end=($_->[1]==1 ? ($bin_size-1) : $_->[1]+($bin_size-1));			#calculate end position based on bin size
		}
		else {																	#otherwise use start value as end value
			$end=$_->[1];
		}
		if($list ne "") {
			die "Error: chromosome $_->[0] could not be found in the chromosome lengths list\n" if !exists($length_table{$_->[0]});
			$end=$length_table{$_->[0]} if $end>$length_table{$_->[0]};			#if end value exceeds the current chromosome length and a length list was provided, set it's value to the chromosome's end
		}
		print OUTFILE "$_->[0]\t$_->[1]\t$end\t$value\n";
	}
	close(OUTFILE);
	return;
}
	
sub print_usage {								#subroutine called when an error in command line input is detected
	my $arg=shift;							
	if($arg eq "b")	{						#if -b option is present and binned output is not requested
		print "Error: binned output files must be requested (-o b or -o v) for -b option to be accepted\n";
	}
	elsif($arg eq "bz") {					#if bin size specified is less than 1
		print "Error: bin size must be 1 or larger\n";		
	}
	elsif($arg eq "s") {					#if -s option is present and merged output is not requested
		print "Error: merged output file must be requested (-o m) for -s option to be accepted\n";
	}
	elsif($arg eq "o") {						#if letters other than n, b, v, or m are included in the -o option's argument
		$option{o}=~s/n|b|m|v//g;					
		print "Error: invalid output file(s) requested: $option{o}\n";
	}
	elsif($arg eq "f" && scalar(@ARGV)>0) {						#if input file name is invalid
		print "Error: could not open input file \"$ARGV[0]\"\n";
	}
	elsif($arg eq "t") {						#if track name is not specified
		print "Error: output file prefix must be specified\n";
	}							
	elsif($arg eq "db") {						#if mysql lookup fails
		print "Error: could not fetch list of chromosome lengths for reference genome \"$list\"\n";
	}
	elsif($arg eq "lf") {						#if chromosome lengths file is invalid
		print "Error: could not open chromosome list file \"$list\"\n";
	}							
	die "Usage: bowtie2bedgraph.pl [options] [input file] [output file prefix]\n\t-o [nbmv]\tspecify additional files to be generated: n=normalized,\n\t\t\tb=binned, m=merged, v=binned and merged\n\t-b [integer>=1]\tspecify bin size, requires -o b or -o v\n\t-s [integer>=1]\tspecify number of bps to shift ChIP-seq hits prior to\n\t\t\tmerging, requires -o m or -o v\n\t-l [string]\tprovide list of chromosome names and lengths to prevent\n\t\t\thits from being shifted beyond the end of the\n\t\t\tchromosome, requires -o b, -o m, or -o v\n\t-r [string]\trather than specifying a chromosome length list,\n\t\t\tspecify a UCSC reference genome identifier (e.g. mm9),\n\t\t\tto have chromosome lengths fetched automatically,\n\t\t\trequires -o b, -o m, or -o v\n\t-t [integer>=1]\tspecify number of nt trimmed prior to alignment\n\t-x\t\trequest stranded output files be swapped,\n\t\t\tforward->reverse and reverse->forward\n\t-D\t\tdisable appending \"chr\" to beginning of chromosome names\n\t\t\tand renaming of mitochondrial chromosome\n\t-M [string]\tstring used to match mitochondrial chromosome name,\n\t\t\tdefault: \"mit\"\n";
}
