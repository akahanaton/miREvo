#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use File::Basename qw/basename dirname/;
use Data::Dumper;
use Bio::SeqIO;

my %opts;
my $usage = qq{
	Usage: $0 hairpin_file mature_file kn_bwt db_bwt
		or
	Usage: $0 hairpin_file mature_file kn_bwt db_bwt reads_file 
};

getopts('', \%opts);
my ($hp_file, $ma_file, $kn_bwt, $db_bwt, $reads );
if (@ARGV == 5){
	($hp_file, $ma_file, $kn_bwt, $db_bwt, $reads ) = @ARGV;
}elsif(@ARGV == 3){
	($hp_file, $ma_file, $kn_bwt) = @ARGV;
}else{
	die($usage) 
}

my %readsStat;
if(defined($reads)){
	open(READ, $reads)|| print "Warning: can not open $reads: $!\n";
	while(my $line = <READ>){
		if($line =~ /^>/){
			my(undef,undef,$count) = split(/_/,$line);
			$count =~ s/x// ;
			$readsStat{"Total"} += $count;
		}
	}
	close READ;
	printf "## total Reaads:\t%d\n", $readsStat{"Total"};
}

if(defined($db_bwt)){
	open(DB, $db_bwt)|| print "Warning: can not open $db_bwt: $!\n";
	while(my $line = <DB>){
		my @arr = split(/\s+/, $line);
		my(undef,undef,$count) = split(/_/,$arr[0]);
		$count =~ s/x// ;
		$readsStat{"DB"} += $count;
	}
	close DB;
	printf "## Filtered Reaads\t%d\n", $readsStat{"DB"};
}

if ( ! -e $kn_bwt){
	print "No mapping result for known miRNA, exit expresion analysis.\n";
	exit;
}

my %mature_seq;
my %hp_id;
my $seqin=Bio::SeqIO->new(-file=>$hp_file,-format=>"fasta");
while(my $seq=$seqin->next_seq()){
	$hp_id{lc($seq->display_id)} = 1;
}

my $matin=Bio::SeqIO->new(-file=>$ma_file,-format=>"fasta");
while(my $seq=$matin->next_seq()){
	my ($mature,$mature_tag);
	if ($seq->display_id() =~ /\*$/){
		my $new_seq = lc($seq->seq());
		$new_seq =~ tr/u/t/ ;
		$mature_seq{lc($seq->display_id)}{'star'} = $new_seq;
		next;
	}
	if(defined($hp_id{lc $seq->display_id})){
		my $new_seq = lc($seq->seq());
		$new_seq =~ tr/u/t/ ;
		$mature_seq{lc($seq->display_id)}{'mature'} = $new_seq;
		next;
	}
	my ($mir, $tag) = ($seq->display_id =~ m/(.*)(-[35]p)/i);
	if(defined($hp_id{lc $mir})){
		my $new_seq = lc($seq->seq());
		$new_seq =~ tr/uU/tT/ ;
		$mature_seq{lc($mir)}{lc $tag} = $new_seq;
	}
}

#read tags mapping information
my ($tags_info, %mapped_kn_reads, $total_mapped_copy);
# RL1_51_x5211    +       run1_45_gtggactgttgtcggccgcgcc_4.5      0       GTGGACTGTTGTCGGCCGCGCC  IIIIIIIIIIIIIIIIIIIIII  0
open(TAG,$kn_bwt);
while(my $line=<TAG>){
	chomp($line);
	my @ary=split(/\t/,$line);
	my ($lib,$ii,$copy)=split(/_/,$ary[0]);
	$copy =~ s/x//;
	# multiple hit
	$copy = $copy / ($ary[6] + 1);
	my $elem={"copy"=>$copy,"seq"=>$ary[4],"lib"=>$lib,"strand"=>$ary[1], "start" => $ary[3], "lng" => length $ary[4]};
	push(@{$tags_info->{lc($ary[2])}},$elem);
	$mapped_kn_reads{$ary[0]} += 1;
}

foreach my $read (keys %mapped_kn_reads){
	my(undef,undef,$count) = split(/_/,$read);
	$count =~ s/x// ;
	$total_mapped_copy += $count;
}
print "\n\n";
printf "## Mapped known miRNA Reads:\t%d\n", $total_mapped_copy;

my %mir_exp;
$seqin=Bio::SeqIO->new(-file=>$hp_file,-format=>"fasta");
while(my $seq=$seqin->next_seq()){	
	my $acc=lc($seq->display_id());
	foreach my $mat_tag( keys %{$mature_seq{$acc}}){
		my ( $mature_beg, $mature_end , $new_mat_tag) = get_mature_in_hairpin_seq($seq->seq(), $acc, $mat_tag, \%mature_seq);
		if (defined($tags_info->{$acc} ) && defined( $mature_beg )){
			my @tags=@{$tags_info->{$acc}};
			$mir_exp{$acc}{$new_mat_tag}=0;
			$mir_exp{$acc}{hp}=0;
			foreach my $elem(@tags){
				my $tagstart = $elem->{"start"};
				if ($tagstart >= ($mature_beg-3) && ($tagstart+$elem->{"lng"}) <= ($mature_end+3) && $elem->{"strand"} eq "+"){
					#--------------------------------------------------
					# print "$acc\t($mature_beg-3)\t($mature_lng+3)\t$tagstart\t",$elem->{"acc"},"\t",$elem->{"copy"},"\n";
					#-------------------------------------------------- 
					$mir_exp{$acc}{$new_mat_tag} = $mir_exp{$acc}{$new_mat_tag} + $elem->{"copy"};
					$mir_exp{$acc}{hp} = $mir_exp{$acc}{hp} + $elem->{"copy"};
				}else{
					$mir_exp{$acc}{hp} = $mir_exp{$acc}{hp} + $elem->{"copy"};
				}

			}
		}
	}
}


printf "## mirID\t5p_count\t3p_count\thairpin_count\n";

foreach my $acc( keys %mir_exp){
	my @tags = keys %{$mir_exp{$acc}};
	$mir_exp{$acc}{'-5p'} = "undef" if ! defined $mir_exp{$acc}{'-5p'};
	$mir_exp{$acc}{'-3p'} = "undef" if ! defined $mir_exp{$acc}{'-3p'};
	print "$acc", "\t", $mir_exp{$acc}{'-5p'}, "\t", $mir_exp{$acc}{'-3p'},"\t",$mir_exp{$acc}{hp},"\n";
}


sub rev_com_seq {
	my $seq = shift ;
	$seq =~ tr/atcgnATCGN/tagcnTAGCN/ ;
	$seq =reverse($seq);
	return($seq);
}

sub get_mature_in_hairpin_seq {
	my ($hairpin_seq, $hairpin_id, $tag, $mature_info) = @_;
	$hairpin_seq = lc($hairpin_seq);
	$hairpin_seq =~ tr/u/t/;
	my $hairpin_lng = length $hairpin_seq;
	my $id = lc($hairpin_id);
	my ($mature_beg, $mature_end) = (undef, undef);
	if(defined($mature_info->{$id}->{$tag})){
		my $mat_seq = $mature_info->{$id}->{$tag};
		my $mature_lng = length $mat_seq;
		my $tmp_pos;
		$tmp_pos = index($hairpin_seq, $mat_seq);
		if ($tmp_pos < 0){
			my $mat_seq_rec = rev_com_seq($mat_seq);
			$tmp_pos = index($hairpin_seq, $mat_seq);
			if ($tmp_pos < 0){
				print STDERR "Conn't defined mature position for $id$tag\t$mat_seq\t$hairpin_seq\n";
			}
		}
		$mature_beg = $tmp_pos;
		if( $mature_beg > ($hairpin_lng / 2)){
			$tag = '-3p';
		}else{
			$tag = '-5p';
		}
		$mature_end = $tmp_pos + $mature_lng ;
	}
	return ($mature_beg, $mature_end, $tag);
}
