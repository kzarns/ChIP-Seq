#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

my $usage = join("\n",
	"perl cppmatch_multi <db_file_name> <cpp_file_name_1> .. <cpp_file_name_n>",
	"",
	"map cpp files onto a db file",
);

die "not enough stuff: " . $usage unless $#ARGV > 0;
my $cpp_file_names = [];
my $cpp_files = [];
print STDERR "argv: " . Dumper(\@ARGV);

my $db_file_name = $ARGV[0];
my $file_count = scalar @ARGV;

foreach my $key (@ARGV[1 .. $file_count - 1]){
	push (@$cpp_file_names, $key);
}

print STDERR Dumper $cpp_file_names;

my $master = {};
open(my $db_file, $db_file_name) or die "failed opening $db_file_name";
readline($db_file);

while(my $line = readline($db_file)){
	chomp($line);
	my ($id, $symbol, $chrom, $tss, $tes, $strand) = split("\t", $line);

	#using $id creates a lot of duplicates	
	my $key = $id;
	#my $key = $symbol;
	#take the first
	unless($master->{$chrom}{$key}){
		my $hits;
		@$hits = (0) x ($file_count - 1);
		#print Dumper $hits;
		$master->{$chrom}{$key} = {"id" => $id, "symbol" => $symbol, "tss" => $tss, "tes" => $tes, "hits" => [@$hits]};
	}
}
close($db_file);
print STDERR "read in master\n";

#print STDERR "master genes: " . Dumper($master);

for(my $i = 0; $i < @$cpp_file_names; $i++){
	associate($i);
}

sub associate{
	my $index = shift;
	my $name = $cpp_file_names->[$index];

	open(my $file, $name) or die "failed to open $name: $!";

	while(my $line = readline($file)){
		#print $line;
		chomp($line);

		my ($id, $hits, $chrom, $tss, $tes) = split("\t", $line);
		foreach my $key (keys %{$master->{$chrom}}){
			my $gene = $master->{$chrom}{$key};
			if($gene->{"tss"} <= $tss and $gene->{"tes"} >= $tss){
				$gene->{"hits"}[$index] += $hits;
			}
		}
	}

	close($file);
	print STDERR "opened processed and closed $name\n";
}

my $cpp_file_name = join("_", @$cpp_file_names, ".cpp_result");
open(CPP, $cpp_file_name) or die "failed to open $cpp_file_name: $!";
open(CPP_TOTAL, $cpp_file_name . "_total") or die("failed to open $cpp_file_name" . "_total: $!");

print CPP join("\t",
	"gene_id",
	"gene_symbol:tss",
	join("\t", @{$cpp_file_names}),
) . "\n";

my $fields = [];
foreach my $file_name (@$cpp_file_names){
	my $sub_fields = ["desc1", "hits", "chr", "start", "end"];
	foreach my $sub (@$sub_fields){
		push(@$fields, "$file_name.$sub");
	}
}

print CPP_TOTAL join("\t", "db.desc1", "db.desc2", "db.chr", "db.start", "db.end", @$fields);

foreach my $chrom (keys %$master){
	foreach my $key (keys %{$master->{$chrom}}){
		my $gene = $master->{$chrom}{$key};
		print CPP join("\t",
			$gene->{"id"},
			$gene->{"symbol"},
			join("\t", @{$gene->{"hits"}})
		) . "\n";
	}
}

1;
