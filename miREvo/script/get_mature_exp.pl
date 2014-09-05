#!/user/bin/perl -w
use strict;
use warnings;
use Bio::SeqIO;

if(scalar(@ARGV)<3){
	die "$0 hairpin.fa  mature.fa soap_result.out\n";
}



#read tags information
my $mapping_result=$ARGV[2]; #soap mapping results
my $tags_info;
# RL1_51_x5211    +       run1_45_gtggactgttgtcggccgcgcc_4.5      0       GTGGACTGTTGTCGGCCGCGCC  IIIIIIIIIIIIIIIIIIIIII  0
open(TAG,$mapping_result);
while(my $line=<TAG>){
	chomp($line);
	my @ary=split(/\t/,$line);
	my ($lib,$ii,$copy)=split(/_/,$ary[0]);
	$copy =~ s/x//;
	my $elem={"copy"=>$copy,"seq"=>$ary[4],"acc"=>$lib,"strand"=>$ary[1], "start" => $ary[3], "lng" => length $ary[4]};
	push(@{$tags_info->{uc($ary[2])}},$elem);
}

my $mature_file=$ARGV[1];
my %mature_seq;
my $matin=Bio::SeqIO->new(-file=>$mature_file,-format=>"fasta");
while(my $seq=$matin->next_seq()){
	my ($mature,$mature_tag);
	next if ($seq->display_id() =~ /\*$/);
	if($seq->display_id() =~ m/(.*mir\d+\w*)(.*)/i ){
		$mature = $1;
		$mature_tag = $2;
	}else{
		$mature = $seq->display_id();
		$mature_tag = '';
	}
	my $new_seq = lc($seq->seq());
	$new_seq =~ tr/uU/tT/ ;
	#--------------------------------------------------
	# print $mature_tag,"\t",$mature,"\n";
	#-------------------------------------------------- 
	if( $mature_tag eq '' ){
		$mature_seq{uc($mature)}{'mature'} = $new_seq;
	}else{
		$mature_seq{uc($mature)}{$mature_tag} = $new_seq;
	}
}

my $hp_file=$ARGV[0];
my $seqin=Bio::SeqIO->new(-file=>$hp_file,-format=>"fasta");
while(my $seq=$seqin->next_seq()){	
	my $acc=uc($seq->display_id());
	my @all_tags;
	my ($mature_tag_copy, $all_tag_copy) = (0,0);
	my ($mature_hairpin,$id_info, $mature_beg, $mature_lng ) = get_mature_in_hairpin_seq($seq->seq(), $acc, \%mature_seq);
	if (defined($tags_info->{$acc})){
		my @tags=@{$tags_info->{$acc}};
		my @tmp_all_tags;
		foreach my $elem(@tags){

			my $tagstart = $elem->{"start"};
			#--------------------------------------------------
			# print $tagstart,"\t",$mature_beg,"\t",$mature_lng,"\n";
			#-------------------------------------------------- 
			if ($tagstart >= ($mature_beg-3) && ($tagstart+$elem->{"lng"}) <= ($mature_beg + $mature_lng+3) && $elem->{"strand"} eq "+"){
				#--------------------------------------------------
				# print "$acc\t($mature_beg-3)\t($mature_lng+3)\t$tagstart\t",$elem->{"acc"},"\t",$elem->{"copy"},"\n";
				#-------------------------------------------------- 
				$mature_tag_copy = $mature_tag_copy + $elem->{"copy"};
				$all_tag_copy = $all_tag_copy + $elem->{"copy"};
			}else{
				$all_tag_copy = $all_tag_copy + $elem->{"copy"};
			}

		}
	}

	print "$acc\t$mature_tag_copy\t$all_tag_copy\n";
	
}

sub get_mature_in_hairpin_seq {
	my ($hairpin_seq, $hairpin_id, $mature_info) = @_;
	$hairpin_seq = lc($hairpin_seq);
	my $id = uc($hairpin_id);
	my (@mature_pos, @mature_lng, @mature_seq, @mature_id);
	#--------------------------------------------------
	# print Dumper($mature_info);
	#-------------------------------------------------- 
	if(defined($mature_info->{$id})){
		foreach my $tag (keys %{$mature_info->{$id}}){
			my $mat_seq = $mature_info->{$id}->{$tag};
			my $tmp_pos;
			$tmp_pos = index($hairpin_seq, $mat_seq);
			if ($tmp_pos < 0){
				my $mat_seq_rec = rev_com_seq($mat_seq);
				$tmp_pos = index($hairpin_seq, $mat_seq);
				if ($tmp_pos < 0){
					print STDERR "$id$tag\t$mat_seq\t$hairpin_seq\n"
				}else{
					push (@mature_pos, $tmp_pos);
					push (@mature_lng, length $mat_seq);
					push (@mature_seq, $mat_seq_rec);
					$tag =~ s/(mature)//;
					push (@mature_id, "$id$tag");
				}
			}else{
				push (@mature_pos, $tmp_pos);
				push (@mature_lng, length $mat_seq);
				push (@mature_seq, $mat_seq);
				$tag =~ s/(mature)//;
				push (@mature_id, "$id$tag");
			}
		}
	}

	my $mature_id_info = $hairpin_seq;
	$mature_id_info =~ tr/atcgnATCGN/ /;

	for (my $i =0; $i <= $#mature_seq; $i++){
		#--------------------------------------------------
		# print uc($mature_seq[$i]),"\n";
		#-------------------------------------------------- 
		substr($hairpin_seq, $mature_pos[$i], $mature_lng[$i]-1, uc($mature_seq[$i]));
		substr($mature_id_info, $mature_pos[$i], length($mature_id[$i])-1, lc($mature_id[$i]));
	}
	return ($hairpin_seq, $mature_id_info, $mature_pos[0], $mature_lng[0]);
}

sub rev_com_seq {
	my $seq = shift ;
	$seq =~ tr/atcgnATCGN/tagcnTAGCN/ ;
	$seq =reverse($seq);
	return($seq);
}
