use strict;
use warnings;

open(GFF,$ARGV[0]);

while(my $line=<GFF>){
	chomp($line);
	my @ary=split(/\t/,$line);
	print "chr".$ary[0]."\t".$ary[3]."\t".$ary[4]."\t".$ary[8]."\n";
}
