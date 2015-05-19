#!/usr/bin/perl

#bowtie2bedgraph v2.0
#
#This script originated at the NIH
#UND students, faculty, and staff have permission to use the script
#The script has been refactored for performance and readability by Kris Zarns

#generate warnings
use warnings;
#generate warnings when undefined variables are referenced
use strict;
#allow use of threads
use threads;
#use get options package
use Getopt::Std;

use Fcntl qw(:seek);

#initialize variables
my (%plus_table, %minus_table, %plus_bin_table, %minus_bin_table, %merge_table, %option, %length_table);
my ($normalized, $binned, $merged, $bin_size, $shift_val, $trim, $list, $minus_shift_val);
my $norm_value = 0;

#extract options from the command line argument
getopts('o:b:s:t:l:n:', \%option) || print_usage();

#allow -o, -b, and -s, if other arguments are passed, run print_usage subroutine
#if the -o option is passed
if(exists($option{o})){
	#use -o to select additional output file types
	#check for characters other than n, b, or m in its argument
	if($option{o} =~ /[^nbm]/){
		#if other characters are found, run the print_usage subroutine
		print_usage("o");
	}
	#set flag values for each optional file time,
	$normalized = ($option{o} =~ /n/i) || 0;
	$binned = ($option{o} =~ /b/i) || 0;
	$merged = ($option{o} =~ /m/i) || 0;
}

#if shift value or bin size is less than 1, run the print_usage subroutine
print_usage("bz") if (exists($option{b}) && $option{b} < 1);
print_usage("sz") if (exists($option{s}) && $option{s} < 1);
#if bin size is specified and binned output is not requested, run print_usage
print_usage("b") if (!$binned && !$merged && exists($option{b}));
#if shift value is specified and merged output is not requested, run print_usage
print_usage("s") if (!$merged && exists($option{s}));

#if binned output is requested, set bin size to specified value, or 25, otherwise set it to 1
$bin_size = ($binned ? ($option{b} || 25) : 1);
#if shift value is not specified, set to 75
$shift_val = $option{s} || 75;

#minus only shift val
$minus_shift_val = $option{m} || 0;

#if trim value is not specified, set $trim to 0
$trim = $option{t} || 0;

if(exists($option{l})){

	$list = $option{l};
	#if merged output is not requested, and a list file name is specified, run print_usage
	print_usage("l") if (!$merged);
	#open chromosome list file, run print_usage subroutine if file cannot be opened
	open(LISTFILE, $list) || print_usage("lf");
	#split each line by whitespace and store in the length_table hash, with chromosome name as key and length as value
	while(<LISTFILE>){

		/(\S+)\s+(\S+)/;
		$length_table{$1} = $2;
	}
	close(LISTFILE);
}else{
	#if list file name is not specified, assign blank value
	$list = "";
}

my $source_name = $ARGV[0];
my $sorted_name = $source_name . "_sorted.bowtie";
#for bedgraph conversion, the order of the bowtie file is not important
#however, if we presort by chromosome we can assume that once we've run
#out of a chromosome in the file we will not see another occurance of it
#and can move on to parsing for the next chromosome
#this is a huge performance win since we don't have to keep all
#chromosomes in a hash at once

#always re-sort.  We hit a bug where an empty file created by a failed sort was being reused
my $sort_start = time();
print "sorting first into $sorted_name\n";
`sort -t\$'\t' -k 3,3V -S 10% $source_name.bowtie > $sorted_name`;
print "sorting took: " . (time() - $sort_start) . "\n";

#hit values will be normalized to millions of overall hits
#count number of lines
$norm_value = `wc -l < $source_name.bowtie`;
$norm_value /= 1000000;
print "norm: $norm_value\n";

#check for name in argument, run print_usage subroutine if name is not specified
my $track_name=$ARGV[1] || print_usage("t");
print "name: $track_name\n";

my $found_chrome = {};
my $this_chrome = undef;
my $found_all = 0;

#open input file, run print_usage subroutine if file cannot be opened
open(INFILE, "$sorted_name") || print_usage("f");

#prepare the headers for the various files
if(!$merged){
	initializer($track_name, "forward", $normalized);
	initializer($track_name, "reverse", $normalized);

	if($binned){
		initializer($track_name, "forward_binned", $normalized);
		initializer($track_name, "reverse_binned", $normalized);
	}
}
if($merged){
	initializer($track_name, "merged", $normalized);
}

#just debugging
print "argv: $#ARGV\n";
foreach my $key (0 .. $#ARGV){
	print "arg: $key, val: $ARGV[$key]\n";
}

#track the time for the starting sort, the processing, and the ending sorts separately
my $process_time = time();

while( not $found_all){
#Start doing some work


	#loop over input lines in chunks of one chromosome at a time
	while(my $line = <INFILE>){

		#skip lines starting with @
		next if $line =~ /^@/;

#example line
#SNPSTER3_0697:6:1:3076:1009#0/1	4	*	0	0	*	*	0	0	NATCCCTGTTTGGCTTAAGAAGTGTTCCAAGACCATGGAATACTGTGGAGATCG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	XM:i:0
#SNPSTER3_0697:4:1:7926:1010#0/1 run=101028_SNPSTER3_0697_62N7MAAXX	+	chr11	119103480	NGGATTTGCCTCCTAAACGCATTTTCCAGTTCATTCCCCAGCACAATATGCCTG	BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB	0	0:C>N
#  bowtie
#  0 name of the read SNPSTER... note the space here and account for it -- skip
#  1 strand alignment type, + --
#  2 name of sequence where alignment occurs chr11
#  3 0 based offset in strand forward where leftmost character of the allignment occurs 119103480
#  4 read sequence NGGATTT...
#  5 read qualities BBB...
#  6 ceiling 0
#  7 mismatches 0:C>N
#@cols=($line =~/^\S+\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t/); #extract the second through fourth fields of the current line
#                0       1      2      3      4  bowtie
#                -1      0      1      2      3  cols

		#extract the second through fourth fields of the current line
		my ($alignment, $sequence_name, $forward_offset, $sequence_read) = ($line =~ /^[^\t]+\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t/);

		#HWUSI-EAS1676_0053_FC:6:1:1874:951#0/1	+	chr1	44336077	NCTGGGAGCAGGGAAGCCCAGGTCTGGGGAGGGGAA	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%	0	0:G>N
		die "this didn't match: $line" unless defined $alignment;

		$sequence_name="chrM" if $sequence_name =~ m/mito/i;
		#if the chromosome field does not begin with 'chr', add it
		if($sequence_name !~ m/chr/i){
			$sequence_name = "chr" . $sequence_name;
		}

		#are we searching for our next chromosome?
		#this could probably be optimized out
		if(not defined($this_chrome)){
			#is this a good chrome to start reading for?
			if(not exists($found_chrome->{$sequence_name})){
				$this_chrome = $sequence_name;
				print "start processing a new chrome: $this_chrome\n";
				$found_chrome->{$sequence_name} = 1
			}if($found_chrome->{$sequence_name} eq 0){
				$this_chrome = $sequence_name;
				print "start processing a new chrome: $this_chrome\n";
				$found_chrome->{$sequence_name} = 1
			}else{
				#we've already processed this one, try another line
				next;
			}
		}

		#we have a chrome we are looking for, is this it?
		#again, with a bit more effort we could optimize this out mostly
		if(not $this_chrome eq $sequence_name){

			#it's not, have we seen this one before?
			if(not exists($found_chrome->{$sequence_name})){
				#we haven't, record it
				#this is probably redundent
				$found_chrome->{$sequence_name} = 0;
			}

			#backup to the start of this line so we can read it in again on the next
			#run of this loop
			seek(INFILE, -length($line), SEEK_CUR);
			last;
		}

		#we found what we were looking for, run with it
		#if hit is to forward strand
		#this section is largely unchanged, updated variable names
		if($alignment eq '+'){

			#adjust position from continuum starting at 0 to continuum starting at 1 and adjust for trimming
			$forward_offset += (1 - $trim);
			#concatenate the chromosome and location, and use it as a hash key, allows quick locating of previously hit locations
			my $key = "$sequence_name\t$forward_offset";
			#increase hit count for that chromosome and location
			$plus_table{$key}++;
			
			#create our own length table as best we can
			if($list eq ""){
				if(not exists $length_table{$sequence_name} or $length_table{$sequence_name} < $forward_offset){
					$length_table{$sequence_name} = $forward_offset;
				}
			}
			
			#if binned output is requested
			if($binned){

				#determine alignment's bin, use 1 if calculated value is zero
				my $temp = int($forward_offset / $bin_size) * $bin_size;
				$temp = 1 if $temp < 1;
				#concatenate the chromosome and bin start location, and use it as a hash key
				$key = "$sequence_name\t$temp";
				#increase bin hit count
				$plus_bin_table{$key}++;
			}
			#if merged output is requested
			if($merged){

				#shift alignment specified number of bps downstream

				#determine alignment's bin
				my $temp = int(($forward_offset + $shift_val) / $bin_size) * $bin_size;
				$temp = 1 if $temp < 1;
				#if an entry in the length table exists for the current chromosome
				if(exists($length_table{$sequence_name})){

					#check if the alignment was shifted beyond the end of the chromosome, if so, set position to end
					$temp = int($length_table{$sequence_name} / $bin_size) * $bin_size if $temp > $length_table{$sequence_name};
				}
				#if an entry wasn't found, and the -l option was specified, print an error message and exit
				elsif($list ne ""){

					die "Error: chromosome $sequence_name could not be found in the chromosome list file\n";
				}

				#concatenate the chromosome and bin start location, and use it as a hash key
				#increase bin hit count
				$merge_table{"$sequence_name\t$temp"}++;
			}
		#if hit is to reverse strand
		}elsif($alignment eq '-'){

			#adjust position by adding the length of the aligned read, now specifies 5' location, also add trim value
			$forward_offset += (length($sequence_read) + $trim);

			#create our own length table as best we can
			if($list eq ""){
				if(not exists $length_table{$sequence_name} or $length_table{$sequence_name} < $forward_offset){
					$length_table{$sequence_name} = $forward_offset;
				}
			}

			#perform same tasks as above
			my $key = "$sequence_name\t$forward_offset";
			$minus_table{$key}++;
			if($binned){
				my $temp = int(($forward_offset + $minus_shift_val) / $bin_size) * $bin_size;
				$temp = 1 if $temp < 1;
				$key = "$sequence_name\t$temp";
				$minus_bin_table{$key}++;
			}

			if($merged){
				my $temp = int(($forward_offset - $shift_val - $minus_shift_val) / $bin_size) * $bin_size;
				$temp = 1 if $temp < 1;
				if(exists($length_table{$sequence_name})){

					$temp = int($length_table{$sequence_name} / $bin_size) * $bin_size if $temp > $length_table{$sequence_name};
				}elsif($list ne ""){

					die "Error: chromosome $sequence_name could not be found in the chromosome list file\n";
				}

				#concatenate the chromosome and bin start location, and use it as a hash key
				#increase bin hit count
				$merge_table{"$sequence_name\t$temp"}++;
			}
		}
	}

	#process our results for this chromosome
	print "starting threads for $this_chrome\n";
	my $plus_thread;
	my $minus_thread;
	my $plus_bin_thread;
	my $minus_bin_thread;
	if(!$merged){
		$plus_thread = threads->create('print_table', \%plus_table, $track_name, "forward", $normalized);
		$minus_thread = threads->create('print_table', \%minus_table, $track_name, "reverse", $normalized);
		$plus_bin_thread = threads->create('print_table_alt', \%plus_bin_table, $track_name, "forward_binned") if $binned;
		$minus_bin_thread = threads->create('print_table_alt', \%minus_bin_table, $track_name, "reverse_binned") if $binned;
	}
	my $merge_thread = threads->create('print_table_alt', \%merge_table, $track_name, "merged") if $merged;

	print "threads created, waiting on join\n";
	#join threads once all the work is done for this chromosome
	if(!$merged){
		$plus_thread->join();
		$minus_thread->join();
		$plus_bin_thread->join() if $binned;
		$minus_bin_thread->join() if $binned;
	}
	$merge_thread->join() if $merged;
	print "threads joined\n";

	#cleanup
	#we reinitialize these tables for each chromosome
	#this is a huge memory savings
	undef %plus_table;
	undef %minus_table;
	undef %plus_bin_table;
	undef %minus_bin_table;
	undef %merge_table;
	undef $this_chrome;

	#we hit the end, should we start over?
	$found_all = 1;
	foreach my $chrome (keys(%$found_chrome)){
		if(not $found_chrome->{$chrome}){
			$found_all = 0;
			print "still missing $chrome, run another round\n";
		}
	}
}
#close input file
close(INFILE);

print "process time: " . (time() - $process_time) . "\n";
print "writing results to file and sorting results\n";

my $final_sort_time = time();
#use unix sort to sort the output as well
if(!$merged){
	finalizer($track_name, "forward", $normalized);
	finalizer($track_name, "reverse", $normalized);

	if($binned){
		finalizer($track_name, "forward_binned", $normalized);
		finalizer($track_name, "reverse_binned", $normalized);
	}
}
if($merged){
	finalizer($track_name, "merged", $normalized);
}

print "finalize sort time: " . (time() - $final_sort_time) . "\n";

sub initializer{

	my ($track_name, $strand, $normalized) = @_;

	#clobber
	open(my $out_file, ">", "${track_name}_${strand}\.bedgraph") || die "Could not create output file: ${track_name}_${strand}\.bedgraph\n";

	#create output file for standard hit count
	#print track information to standard output file
	print $out_file "track type=bedGraph name=kz_${track_name}_${strand} description=kz_${track_name}_${strand} visibility=full color=179,27,27 altColor=179,27,27 priority=20\n\n";
	close($out_file);

	if($normalized){

		open(my $norm_file, ">", "${track_name}_${strand}_normal\.bedgraph") || die "Could not create output file: ${track_name}_${strand}_normal\.bedgraph\n";
		#print track information to normalized file
		print $norm_file "track type=bedGraph name=kz_${track_name}_${strand}_normal description=kz_${track_name}_${strand}_normal visibility=full color=179,27,27 altColor=179,27,27 priority=20\n\n";
		close($norm_file);
	}
}

sub finalizer{

	my ($track_name, $strand, $normalized) = @_;

	#clobber
	my $data_file_name = "${track_name}_${strand}\.bedgraph";
	#first two lines are header, preserve those, sort the rest
	#sort on first field, chromosome version sort
	#sort on second field, start point numerically
	#use 10% of ram or less
	`(head -n 2 $data_file_name && tail -n +3 $data_file_name | sort -k 1,1V -k 2,2n -S 10%) > $data_file_name.tmp && mv $data_file_name.tmp $data_file_name`;

	if($normalized){

		my $norm_file_name = "${track_name}_${strand}_normal\.bedgraph";
		`(head -n 2 $norm_file_name && tail -n +3 $norm_file_name | sort -k 1,1V -k 2,2n -S 10%) > $norm_file_name.tmp && mv $norm_file_name.tmp $norm_file_name`;
	}
}

sub print_table{

	my ($ref, $track_name, $strand, $normalized) = @_;

	#append
	open(my $out_file, ">>", "${track_name}_${strand}\.bedgraph") || die "Could not create output file: ${track_name}_${strand}\.bedgraph\n";

	if($normalized){

		open(my $norm_file, ">>", "${track_name}_${strand}_normal\.bedgraph") || die "Could not create output file: ${track_name}_${strand}_normal\.bedgraph\n";
		foreach my $chrome (keys %$ref){
			my ($chrome_name, $chrome_start) = split("\t", $chrome);
			print $out_file "$chrome_name\t$chrome_start\t$chrome_start\t$ref->{$chrome}\n";
			print $norm_file "$chrome_name\t$chrome_start\t$chrome_start\t" . ($ref->{$chrome}/$norm_value) . "\n";
		}
		close($out_file);
	}else{
		foreach my $chrome (keys %$ref){
			my ($chrome_name, $chrome_start) = split("\t", $chrome);
			print $out_file "$chrome_name\t$chrome_start\t$chrome_start\t$ref->{$chrome}\n";
		}
	}
	close($out_file);
}

#subroutine used by bin and merge threads for printing
sub print_table_alt{

	my ($ref, $track_name, $strand) = @_;

	#create binned or merged output file
	open(my $out_file, ">>", "${track_name}_${strand}\.bedgraph") || die "Could not create output file: ${track_name}_${strand}\.bedgraph\n";

	foreach my $chrome (keys %$ref){
		#for each key in the result array, retrieve hit value by rejoining the key components
		my $value = $ref->{$chrome};
		my ($chrome_name, $chrome_value) = split("\t", $chrome);

		my $end;
		#calculate end position based on bin size
		if($chrome_value == 1){
			$end = $chrome_value + ($bin_size - 2);
		}else{
			$end = $chrome_value + ($bin_size - 1);
		}

		#print "track_name: $track_name, strand: $strand, bin_size: $bin_size, chrome_value: $chrome_value, end: $end\n";
		#if end value exceeds the current chromosome length, set it's value to the chromosome's end
		#if($strand eq "merged"){

		#if the length table wasn'tprovided, we created our own length table because it seems binning needs one
		$end = $length_table{$chrome_name} if $end > $length_table{$chrome_name};
		#}

		$end = 1 if($end < 1);
		#print chromosome, location, location, and hit value to output file
		print $out_file "$chrome_name\t$chrome_value\t$end\t$value\n";
	}
	#close output files
	close($out_file);
	return;
}

#subroutine called when an error in command line input is detected
sub print_usage{

	#retrieve argument
	my $arg = shift;
	#if -b option is present and binned output is not requested
	if($arg eq "b"){

		print "Binned output files must be requested (-o b) for -b option to be accepted\n";
		#if -s option is present and merged output is not requested
	}elsif($arg eq "s"){

		print "Merged output file must be requested (-o m) for -s option to be accepted\n";
		#if letters other than n, b, or m are included in the -o option's argument
	}elsif($arg eq "o"){

		#remove valid letters so invalid letters only may be printed
		$option{o} =~ s/n|b|m//g;
		print "Invalid output file(s) requested: $option{o}\n";
		#if bin size specified is less than 1
	}elsif($arg eq "bz"){

		print "Bin size must be greater than 0\n";
		#if shift value specified is less than 1
	}elsif($arg eq "bz"){

		print "Shift value must be greater than 0\n";
		#if input file name is invalid
	}elsif($arg eq "f"){

		print "Could not open input file: $source_name.bowtie\n"
		#if track name is not specified
	}elsif($arg eq "t"){

		print "Track name must be specified\n"
	}elsif($arg eq "l"){

		print "Merged output must be requested for -l option to be accepted\n";
	}elsif($arg eq "lf"){

		print "Could not open chromosome list file: $list\n";
		#print small description of proper command use and exit
	}
	die join("\n",
	  	"Usage: bowtie2bedgraph.pl [options] [input file] [track name]",
		"\t-o [nbm]\tspecify additional files to be generated: n=normalized,",
		"\t\t\tb=binned, m=merged",
		"\t-b [number>=1]\tspecify bin size, requires -o b",
		"\t-s [number>=1]\tspecify number of bps to shift ChIp-seq hits prior to",
		"\t\t\tmerging, requires -o m",
		"\t-l [file name]\tprovide list of chromosome names and lengths to prevent",
		"\t\t\thits from being shifted beyond the end of the",
		"\t\t\tchromosome, requires -o m",
		"\t-t [number>=1]\tspecify number of nt trimmed prior to alignment",
		"\t-n [number>=1]\tspecify number of bps to shift the minus strand",
	);
}

1;
