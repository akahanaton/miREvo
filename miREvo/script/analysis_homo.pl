#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use File::Basename qw/basename dirname/;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 <lib.name>  <predct.mirna>  <soap_result>
};

getopts('', \%opts);
die($usage) if (@ARGV != 2);
my ($inFile,$outFile) = @ARGV;


open(OUT, ">$outFile") || print "Warning: can not write to $outFile: $!\n";

my $total;
open(IN, $inFile) || print "Warning: can not open $inFile: $!\n" ;
while(my $line = <IN>){
	my @arr = split(/\s+/, $line);
	my(undef,undef,$count) = split(/_/,$arr[0]);
	$predictExp{$arr[8]} += $count;
	$total += $count;
}
close IN;

printf OUT "## Mapped predict miRNA Reaads:\t%d\n", $total;
foreach my $key (keys %predictExp){
	printf OUT "%s\t%d\n", $key,$predictExp{$key};
}

close OUT;
