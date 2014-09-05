use strict;
use warnings;
use MyMod::Bio::Tools::miRNA;
use MyMod::Bio::Tools::SeqAna;
use MyMod::Bio::Tools::NW;
use Bio::SeqIO;
use Bio::AlignIO;
use Data::Dumper;

my $seqana=MyMod::Bio::Tools::SeqAna->new();
my $mirna_fact=MyMod::Bio::Tools::miRNA->new();


# for mature/seed alignment;
my %otherMatures;
my %otherSeed;
if(defined($ARGV[2])){
	my @sequences;
	my @ids;
	my $seqin=Bio::SeqIO->new(-file=>$ARGV[2],-format=>"fasta");
	while(my $seq=$seqin->next_seq()){
		my $new_seq = lc $seq->seq();
		$new_seq =~ tr/u/t/;	
		$otherMatures{$seq->id} = $new_seq;
		my $curSeed = substr($new_seq, 1, 7);
		if(defined($otherSeed{$curSeed})){
			push (@{$otherSeed{$curSeed}}, $seq->id);
		}else{
			$otherSeed{$curSeed} = ();
			push (@{$otherSeed{$curSeed}}, $seq->id);
		}
	}
}
my @otherIDs = keys %otherMatures;

#--------------------------------------------------
# print Dumper(%otherSeed);
#-------------------------------------------------- 

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

#read tags information
#--------------------------------------------------
# dme_0_x3626256  +       dme_1_mir-1_droSec1     64      TGGAATGTAAAGAAGTATGGAG  IIIIIIIIIIIIIIIIIIIIII  9
#-------------------------------------------------- 
my $mapping_file=$ARGV[1]; #bwt mapping results
my $tags_info;
open(TAG,$mapping_file);
while(my $line=<TAG>){
	chomp($line);
	my @ary=split(/\t/,$line);
	my ($lib,$mir,$copy)=split(/_/,$ary[0]);
	$copy =~ s/x//;
	my $elem={"copy"=>$copy,"seq"=>$ary[4],"lib"=>$lib,"strand"=>$ary[1]};
	push(@{$tags_info->{$ary[2]}},$elem);
}


my $mirna_file=$ARGV[0];
my $seqin=Bio::SeqIO->new(-file=>$mirna_file,-format=>"fasta");
while(my $seq=$seqin->next_seq()){	
	my $acc=$seq->display_id();
	my @tmpArr = split(/_/, $acc);
	my $matureSeq = lc $tmpArr[2];
	my $curSeed= '';
	if (defined($ARGV[3]) && $ARGV[3] == 1){
		$curSeed = substr($matureSeq,1,7);
	}
	my $aln=Bio::SimpleAlign->new();
	my $tmplocseq=Bio::LocatableSeq->new(
		-display_id=>$seq->display_id(),
		-seq=>uc($seq->seq()),
		-start=>1,
		-end=>$seq->length()
	);
	$aln->add_seq($tmplocseq);
	my ($sst_aln)=$mirna_fact->get_sst_aln($aln,1);
	my @all_tags;
	next if !defined($tags_info->{$acc});
	my @tags=@{$tags_info->{$seq->display_id()}};
	my ($tmpseq)=$aln->each_seq_with_id($acc);
	my @tmp_all_tags;
	foreach my $elem(@tags){
		my $tmptag=Bio::LocatableSeq->new(
			-display_id=>"tag",
			-seq=>$elem->{"seq"},
			-start=>1,
			-end=>length($elem->{"seq"}),
			-strand=>1
		);
		my $tag_rc = reverse $elem->{"seq"};
		$tag_rc =~ tr/ACGTacgt/TGCAtgca/;
		my $tmptag_rc=Bio::LocatableSeq->new(
			-display_id=>"tag",
			-seq=>uc($tag_rc),
			-start=>1,
			-end=>length($tag_rc),
			-strand=>1
			);

		my ($aln_seq,$tagstart);
		($aln_seq,$tagstart)=$seqana->aln_identity_seq($tmpseq,$tmptag);
		if(! defined ($aln_seq)){
			($aln_seq,$tagstart)=$seqana->aln_identity_seq($tmpseq,$tmptag_rc);
		}
		#--------------------------------------------------
		# print $aln_seq,"\n";
		#-------------------------------------------------- 
		next if !defined($aln_seq);
		my $aln_seq_str=$aln_seq->seq();
		$aln_seq_str=~ tr/Tt/Uu/;
		$aln_seq_str=lc($aln_seq_str) if $elem->{"strand"} eq "-";
		#throw out negative strand
		#next if $elem->{"strand"} eq "-";
		my $tag_line=$aln_seq_str." ".$elem->{"lib"}." ".$elem->{"copy"}." ".$elem->{"strand"}."\n";
		push(@tmp_all_tags,{"tag"=>$tag_line,"pos"=>$tagstart});
	}

	my @sort_tags=sort{$a->{"pos"}<=>$b->{"pos"}}@tmp_all_tags;
	grep{push(@all_tags,$_)}@sort_tags;
	#--------------------------------------------------
	# next if scalar(@sort_tags)<1;
	#-------------------------------------------------- 
	print ">$acc\n";

	if(@otherIDs){
		my %matureAnno;
		my %seedAnno;
		if(defined($otherSeed{$curSeed})){
			my @IDs = @{$otherSeed{$curSeed}};
			for(my $i=0; $i <= $#IDs; $i++){
				$seedAnno{$IDs[$i]} = 1;
				my $nw = new Align::NW lc $matureSeq, lc $otherMatures{$IDs[$i]}, $payoff;
				$nw->score;
				$nw->align;
				#--------------------------------------------------
				# $nw->print_align;
				#-------------------------------------------------- 
				my $curIdentity = $nw->get_identity;
				#--------------------------------------------------
				# print $curIdentity,"\n";
				#-------------------------------------------------- 
				if($curIdentity >= 80){
					$matureAnno{$IDs[$i]} = 1;
				}
			}
		}
	
		my @seedAnnoArr = keys %seedAnno;
		if($#seedAnnoArr >= 0){
			my $annoStr = join(" ", @seedAnnoArr);
			print "\nexample miRBase miRNA with the same seed: $annoStr\n";
		}else{
			print "\nexample miRBase miRNA with the same seed: None\n";
		}
		my @matureAnnoArr = keys %matureAnno;
		if($#matureAnnoArr >= 0){
			my $annoStr = join(" ", @matureAnnoArr);
			print "\nexample miRBase miRNA with the similar mature (80% identity): $annoStr\n";
		}else{
			print "\nexample miRBase miRNA with the similar mature (80% identity): None\n";
		}
	}
	
	print "Sequences:\n";
	my ($disseq)=$aln->each_seq_with_id($acc);
	next if !defined($disseq);
	print $tmpseq->seq()."\n";
	
	print "\n";
	
	print "Structures:\n";
	$tmpseq=$sst_aln->{$acc};
	next if !defined($tmpseq);
	print $tmpseq."\n";
	
	print "\n";
	
	print "Reads:\n";
	if( $#all_tags >= 0 ){
		foreach my $tag_elem(@all_tags){
			print $tag_elem->{"tag"};
		}
	}else{
		print "No reads mapped to this miRNA,\n";
	}
	
	print "\n";
	
}
