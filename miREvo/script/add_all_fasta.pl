#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use File::Basename qw/basename dirname/;
use Data::Dumper;
use Bio::SeqIO;
use Bio::AlignIO;

my %opts;
my $usage = qq{
	Usage: $0 <lib.name>  <slim.fas>  <query.fas>
};

getopts('', \%opts);
die($usage) if (@ARGV != 3);
my ($libName,$slim_fas,$query_fas) = @ARGV;


my(@svgFiles,@mirnaTags,%mirnaLen,%mapReadsNum);


my $slim_in=Bio::SeqIO->new(-file=>$slim_fas,-format=>"fasta");
my %mirna_seq;
while(my $seq=$slim_in->next_seq()){
	my ($tmp_acc,$index,$score,$specie)=split(/_/,$seq->display_id());
	my $acc = $tmp_acc.'_'.$index."_".$score;
	$mirna_seq{$acc}{$specie} = $seq->seq();
}

my $query_in=Bio::SeqIO->new(-file=>$query_fas,-format=>"fasta");
while(my $seq=$query_in->next_seq()){
	my @tmpArr=split(/_/,$seq->display_id());
	my $acc = $tmpArr[0].'_'.$tmpArr[1].'_'.$tmpArr[3];
	if(!defined( $mirna_seq{$acc})){
		$mirna_seq{$acc}{"solo"} = $seq->seq();
	}
}
foreach my $acc ( keys %mirna_seq){
	foreach my $specie( keys %{$mirna_seq{$acc}}){
		print ">$acc\_$specie\n";
		print "$mirna_seq{$acc}{$specie}\n";
	}
}

