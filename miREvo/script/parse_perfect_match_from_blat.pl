use strict;
use warnings;
open(OUTFILE,">$ARGV[1]");
open(FILE,$ARGV[0]);
my $count=0;
while(my $line=<FILE>){
	next if $count-->0;
	next if $line=~ m/match/ || $line=~ m/------/ || $line=~ m/psLay/;
	my @ary=split(/\s+/,$line);
	next if !defined($ary[15])||$ary[10]!=$ary[16]-$ary[15] || $ary[10]!=$ary[0];
	$ary[15]++;
	print OUTFILE $ary[9]."\t".$ary[13]."\t".$ary[8]."\t".$ary[15]."\t".$ary[16]."\t".$ary[10]."\t".$ary[14]."\n";
}
