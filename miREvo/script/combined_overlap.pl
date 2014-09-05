#!usr/bin/perl -w
use strict;
use warnings;

open(FILE,$ARGV[0]);
my $last_rec;
my $acc_count;
while(my $line=<FILE>){
	my @ary=split(/\t/,$line);
	my ($lib,$count,$mir,$score)=split(/_/,$ary[0]);
	$acc_count->{$ary[0]}++;
	#--------------------------------------------------
	# $ary[0]=$lib."_".$count."_".$score."_".$acc_count->{$ary[0]};
	#-------------------------------------------------- 
	$ary[0]=$lib."_".$count."_".$score;
	$line=join("\t",@ary);
	my $current_rec;
	$current_rec->{"line"}=$line;
	$current_rec->{"mir"}=$mir;
	$current_rec->{"chr"}=$ary[1];
	$current_rec->{"end"}=$ary[4];
	$current_rec->{"score"}=$score;
	$current_rec->{"len"}=$ary[4]-$ary[3];
	if(defined($last_rec) && $current_rec->{"chr"} eq $last_rec->{"chr"} && $last_rec->{"end"}-$ary[3] > 0.5*$last_rec->{"len"}){
		if ($last_rec->{"mir"} eq $current_rec->{"mir"}){
			if($current_rec->{"score"} !~ /\D/ ){
				if($last_rec->{"score"}<$current_rec->{"score"}){
					$last_rec=$current_rec;
				}
			}
		}
	}elsif(!defined($last_rec)){
		$last_rec=$current_rec;
	}else{
		print $last_rec->{"line"};
		$last_rec=$current_rec;
	}
}
print $last_rec->{"line"};

open(REP,">".$ARGV[1]);

foreach my $acc(keys(%$acc_count)){
	print REP $acc."\t".$acc_count->{$acc}."\n";
}
