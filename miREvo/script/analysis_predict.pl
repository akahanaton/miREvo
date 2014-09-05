#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use File::Basename qw/basename dirname/;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 predict_bwt out_file
};

getopts('', \%opts);
die($usage) if (@ARGV != 2);
my ($pred_bwt, $out_file) = @ARGV;

&analysisPredict( $pred_bwt, $out_file);

sub analysisPredict{
	open(OUT, ">$out_file") || print "Warning: can not write to $out_file: $!\n";

	my $pred_bwt = $_[0];

	my (%predictExp,%mappedPredReads);

	open(PRED, $pred_bwt)|| print "Warning: can not open $pred_bwt: $!\n" ;
	while(my $line = <PRED>){
		my @arr = split(/\s+/, $line);
		my(undef,undef,$count) = split(/_/,$arr[0]);
		$count =~ s/x// ;
		$mappedPredReads{$arr[0]} += 1;
		#--------------------------------------------------
		# dmm_62_GACTGGTGATGGTGTGAACGACGCCC_1.2_GTACCGCCGTCGCCAAGTCTG
		#-------------------------------------------------- 
		my($prefix,$index,$seq, $score) = split(/_/,$arr[2]);
		$predictExp{$index}{prefix} = $prefix;
		$predictExp{$index}{seq} = $seq;
		$predictExp{$index}{count} += $count;
		if(defined($score)){
			$predictExp{$index}{score} = $score;
		}
	}
	close PRED;

	my $total;
	foreach my $read (keys %mappedPredReads){
		my(undef,undef,$count) = split(/_/,$read);
		$count =~ s/x// ;
		$total += $count;
	}

	if ( defined($total)) {
		printf OUT "\n## Mapped known miRNA Reaads:\t%d\n\n", $total;
		printf OUT "## ID\treads_count\n";
		my @mirIndex = sort { $a cmp $b } (keys %predictExp);
		foreach my $i (@mirIndex){
			if(defined($predictExp{$i}{score})){
				print OUT $predictExp{$i}{prefix}."_".$i."_".$predictExp{$i}{seq}."_".$predictExp{$i}{score}."\t".$predictExp{$i}{count}."\n";
			}else{
				print OUT $predictExp{$i}{prefix}."_".$i."_".$predictExp{$i}{seq}."\t".$predictExp{$i}{count}."\n";
			}
		}
	}else{
		printf OUT "\n## Mapped known miRNA Reaads:\t%d\n", 0;
	}

	close OUT;
}
