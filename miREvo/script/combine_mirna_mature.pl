#!/usr/bin/perl -w

use strict;
use warnings;
use MyMod::Bio::Tools::NW;
use Getopt::Std;
use Bio::SeqIO;

my %opts;
my $usage = qq{
	Usage: $0 <lib.name>  <mature.fna>  <hairpin.fna> <other.mature>
};

getopts('', \%opts);

my ($libName,$mat,$hair, $other_file);
if (@ARGV == 3){
	($libName,$mat,$hair) = @ARGV;
}elsif(@ARGV == 4){
	($libName,$mat,$hair, $other_file) = @ARGV;
}else{
	die($usage) 
}

# for mature/seed alignment;
my %otherSeed;
if(defined($other_file)){
	my @sequences;
	my @ids;
	my $seqin=Bio::SeqIO->new(-file=>$other_file,-format=>"fasta");
	while(my $seq=$seqin->next_seq()){
		my ($newID) = ($seq->id =~ /(.*mir-?\d+)/i);
		my $new_seq = lc $seq->seq();
		$otherSeed{substr($new_seq,1,7)} = 1;
	}
}

my $payoff = { match      => 4,  # $match,
			mismatch   => -3, # $mismatch,
			gap_open   => -2, # $gap_open,
			gap_extend => -1 # $gap_extend
};

sub revcomp{
	my $tmpSeq = reverse $_[0];
	$tmpSeq =~ y/CGATUcgatu/GCTAAgctaa/;
	return $tmpSeq;
}

my(%matureSeq);
my $matin=Bio::SeqIO->new(-file=>$mat,-format=>"fasta");
while(my $seq=$matin->next_seq()){	
	my $id;
	next if($seq->display_id() =~ /\*$/);
	if($seq->display_id() =~ /\./){
		my @arr = split(/\./, $seq->display_id());
		$id = $arr[0];
	}elsif( $seq->display_id() =~ /\-[35]p/ ){
		$id = $seq->display_id();
		$id =~ s/\-[35]p//;
	}elsif( $seq->display_id() =~ /\*$/ ){
		$id = $seq->display_id();
		$id =~ s/\*$//;
	}else{
		$id = $seq->display_id();
	}
	my $tmpSeq = lc $seq->seq();
	$tmpSeq =~ tr/u/t/ ;
	#--------------------------------------------------
	# print $id,"\tmature\n";
	#-------------------------------------------------- 
	if(defined($matureSeq{lc($id)})){
		if(defined($otherSeed{substr($tmpSeq,1,7)})){
			$matureSeq{lc($id)} = $tmpSeq; 
		}
	}else{
		$matureSeq{lc($id)} = $tmpSeq; 
	}
}

my $hairin=Bio::SeqIO->new(-file=>$hair,-format=>"fasta");
my $index = 1;
while(my $seq=$hairin->next_seq()){	
	my $id;
	if($seq->display_id() =~ /\./){
		my @arr = split( /\./, $seq->display_id());
		$id = $arr[0];
	}else{
		$id = $seq->display_id();
	}
	if(defined($matureSeq{lc($id)})){
		my $matID;
		if ($seq->display_id() =~ /mir/i){
			$matID = substr($seq->display_id(),4);
		}else{
			my @tmpArr = split(/_/,$seq->display_id());
			if(@tmpArr == 4){
				$matID = $tmpArr[3];
			}else{
				$matID = substr($seq->display_id(),4);
			}
		}
		#--------------------------------------------------
		# print "matID\t$matID\t",$seq->id(),"\n";
		#-------------------------------------------------- 
		print ">$libName"."_".$index."_".$matureSeq{lc($id)}."_".$matID."\n";
		my $tmpSeq = $seq->seq();
		$tmpSeq =~ tr/Uu/Tt/ ;
		print $tmpSeq."\n";
		++$index;
	}else{
		print STDERR "Mature sequence missing: ", $seq->display_id(),"\n";
	}
}

