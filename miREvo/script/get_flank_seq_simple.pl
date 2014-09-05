use strict;
use warnings;

my $genome=$ARGV[1];
my $flank=$ARGV[4]||0;
my ($path,$name)=($genome=~ m/(.*)\/(.*)/sg);
open(FILE,$ARGV[0]);
if(defined($ARGV[2])){
	open(GFF,">".$ARGV[2]);
}
if(defined($ARGV[3])){
	open(BED,">".$ARGV[3]);
}
while(my $line=<FILE>){
	my @ary=split(/\s+/,$line);
	my $loc_str;
	$ary[3]-=$flank;
	$ary[3]=1 if $ary[3]<0;
	$ary[4]+=$flank;
	if($ary[2] eq "+"){
		$loc_str=$ary[1].":".$ary[3]."..".$ary[4];
	}
	else{
		$loc_str=$ary[1].":".$ary[4]."..".$ary[3];
	}
	my @gff_info;
	if(defined($ARGV[3])){
		print BED $ary[1]."\t".$ary[3]."\t".$ary[4]."\t".$ary[0]."\n";
	}
	if(defined($ARGV[2])){
		$gff_info[0]=$ary[1];
		$gff_info[1]="maf";
		$gff_info[2]="mirna";
		$gff_info[3]=$ary[3];
		$gff_info[4]=$ary[4];
		$gff_info[5]=".";
		$gff_info[6]=$ary[2];
		$gff_info[7]=".";
		$gff_info[8]=$ary[0];
		print GFF join("\t",@gff_info)."\n";
	}
		
}
