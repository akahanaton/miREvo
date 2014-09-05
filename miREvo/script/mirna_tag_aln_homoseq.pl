#!/usr/bin/perl -w
use strict;
use warnings;
use MyMod::Bio::Tools::miRNA;
use MyMod::Bio::Tools::SeqAna;
use MyMod::Bio::Tools::NW;
use Bio::SeqIO;
use Bio::AlignIO;

my $seqana=MyMod::Bio::Tools::SeqAna->new();
my $mirna_fact=MyMod::Bio::Tools::miRNA->new();

# for mature/seed alignment;
my %otherMatures;
my %otherSeed;
if(defined($ARGV[4])){
	my @sequences;
	my @ids;
	my $seqin=Bio::SeqIO->new(-file=>$ARGV[4],-format=>"fasta");
	while(my $seq=$seqin->next_seq()){
		my ($newID) = ($seq->id =~ /(.*mir-?\d+)/i);
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

my $major_taxon = $ARGV[3];
my $seqin=Bio::SeqIO->new(-file=>$ARGV[2],-format=>"fasta");
my $seed;
my %mirna_seq;
my %mature_seq;
my $seed_info;
while(my $seq=$seqin->next_seq()){
	#--------------------------------------------------
	# my ($tmp_acc,$index,$mature,$grade)=split(/_/,$seq->display_id());
	#-------------------------------------------------- 
	my @arr = split(/_/,$seq->display_id());
	my $acc = $arr[0].'_'.$arr[1]."_".$arr[3];
	$seed->{$acc} = substr $arr[2],1,7;
	$mirna_seq{$acc} = $seq->seq();
	$mature_seq{$acc} = lc $arr[2];
}



#read tags information
my $mapping_result=$ARGV[1]; #bwt mapping results
my $tags_info;
open(TAG,$mapping_result);
my $total_libs;
while(my $line=<TAG>){
	chomp($line);
	my @ary=split(/\t/,$line);
	my ($lib,$mir,$copy)=split(/_/,$ary[0]);
	$copy =~ s/x// ;
	my $elem={"copy"=>$copy,"seq"=>$ary[4],"lib"=>$lib,"strand"=>$ary[1]};
	my ($lib2, $index, $seq, $score)=split(/_/, $ary[2]);
	my $acc = $lib2."_".$index."_".$score;
	push(@{$tags_info->{$acc}},$elem);
	$total_libs->{$lib}=1;
}

my @all_libs=sort(keys(%$total_libs));
my $align_file=$ARGV[0];
my $alnin=Bio::AlignIO->new(-file=>$align_file,-format=>"clustalw", -displayname_flat=>1);

while(my $aln=$alnin->next_aln()){	
	my $acc;
	my @seqs=$aln->each_seq();
	my @arr = split(/_/,$seqs[0]->display_id());
	$acc = $arr[0]."_".$arr[1]."_".$arr[2];
	my ($sst_aln)=$mirna_fact->get_sst_aln($aln,1);

	my @aln_taxons;
	foreach my $seq(@seqs){
		my @tmpArr = split(/_/, $seq->display_id());
		push(@aln_taxons, $tmpArr[-1]);
	}
	my %seen = ();
	@aln_taxons = grep { ! $seen{$_}++ } @aln_taxons;

	my @all_tags;
	my $hit_acc=$acc."_".$major_taxon;
	next if !defined($tags_info->{$acc});
	my @tags=@{$tags_info->{$acc}};
	#combined tags
	my $all_copy;
	foreach my $tmptag (@tags){
		$all_copy->{$tmptag->{"seq"}."_".$tmptag->{"strand"}}->{$tmptag->{"lib"}}=$tmptag->{"copy"};
	}

	my @new_tags;
	foreach my $tmpseqid(keys %{$all_copy}){
		my ($tmptag,$tmpstrand)=split(/_/,$tmpseqid);
		my @lib_copys;
		foreach my $lib(@all_libs){
			my $tmpcopy=0;
			$tmpcopy=$all_copy->{$tmpseqid}->{$lib} if defined($all_copy->{$tmpseqid}->{$lib});
			push(@lib_copys,$tmpcopy);
		}
		my $new_copy=join("\t",@lib_copys);
		push(@new_tags,{"lib"=>"","copy"=>$new_copy,"seq"=>$tmptag,"strand"=>$tmpstrand});
	}

	my $tmpseq;
	foreach my $seq($aln->each_seq()){
		if($seq->display_id =~ /$major_taxon/){
			$tmpseq=$seq;
		}
	}
	next if !defined($tmpseq);
	$tmpseq->seq(uc($tmpseq->seq()));
	my @tmp_all_tags;
	foreach my $elem(@new_tags){
		my $tmptag=Bio::LocatableSeq->new(
			-display_id=>"tag",
			-seq=>uc($elem->{"seq"}),
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

		next if !defined($aln_seq);
		my $aln_seq_str=$aln_seq->seq();
		$aln_seq_str=~ s/T/U/sg;
		$aln_seq_str=lc($aln_seq_str) if $elem->{"strand"} eq "-";
		my $tag_line=$aln_seq_str." ".$elem->{"lib"}.$elem->{"copy"}." ".$elem->{"strand"}."\n";
		push(@tmp_all_tags,{"tag"=>$tag_line,"pos"=>$tagstart});
	}
	my @sort_tags=sort{$a->{"pos"}<=>$b->{"pos"}}@tmp_all_tags;
	grep{push(@all_tags,$_)}@sort_tags;

	if(defined($seed->{$acc})){
		my $tmpseed=Bio::LocatableSeq->new(
			-display_id=>"seed",
			-seq=>uc($seed->{$acc}),
			-start=>1,
			-end=>length($seed->{$acc}),
			-strand=>1
		);
		my $seed_rc = reverse $seed->{$acc};
		$seed_rc =~ tr/ACGTacgt/TGCAtgca/;
		my $tmpseed_rc=Bio::LocatableSeq->new(
			-display_id=>"seed",
			-seq=>uc($seed_rc),
			-start=>1,
			-end=>length($seed_rc),
			-strand=>1
		);
		my ($aln_seq,$tagstart);
		($aln_seq,$tagstart)=$seqana->aln_identity_seq($tmpseq,$tmpseed);
		if (! defined($aln_seq)){
			($aln_seq,$tagstart)=$seqana->aln_identity_seq($tmpseq,$tmpseed_rc);
		}
		next if !defined($aln_seq);
		my $aln_seq_str=$aln_seq->seq();
		$aln_seq_str=~ s/T/U/sg;
		$seed_info->{$acc}= $aln_seq_str;
	}
	
	print ">$acc\n";
	print "$mirna_seq{$acc}\n";
	
	if(@otherIDs){
		my %matureAnno;
		my %seedAnno;
		my $curSeed = lc substr($mature_seq{$acc} ,1 ,7);
		if(defined($otherSeed{$curSeed})){
			my @IDs = @{$otherSeed{$curSeed}};
			for(my $i=0; $i <= $#IDs; $i++){
				my $nw = new Align::NW lc $mature_seq{$acc}, lc $otherMatures{$IDs[$i]}, $payoff;
				$seedAnno{$IDs[$i]} = 1;
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
			print "\nexample miRBase miRNA with the same seed:$annoStr\n\n";
		}else{
			print "\nexample miRBase miRNA with the same seed:None\n\n";
		}
		my @matureAnnoArr = keys %matureAnno;
		if($#matureAnnoArr >= 0){
			my $annoStr = join(" ", @matureAnnoArr);
			print "\nexample miRBase miRNA with the similar mature (80% identity):$annoStr\n\n";
		}else{
			print "\nexample miRBase miRNA with the similar mature (80% identity):None\n\n";
		}
	}
	
	print "Sequences:\n\n";
	foreach my $other_taxon(@aln_taxons){
		my ($tmpseq)=$aln->each_seq_with_id($acc."_".$other_taxon);
		next if !defined($tmpseq);
		print $tmpseq->seq()."\t\t".$other_taxon."\n";
	}
	print "\n\n";
	
	print "Seed:\n\n";
	if(defined($seed_info->{$acc})){
		print $seed_info->{$acc}
	}
	print "\n\n";

	print "Structures:\n\n";
	foreach my $other_taxon(@aln_taxons){
		#--------------------------------------------------
		# my ($tmpseq)=$sst_aln->each_seq_with_id($acc."_".$other_taxon."_"."sst");
		# next if !defined($tmpseq);
		# print $tmpseq->seq()."\t\t".$other_taxon."\n";
		#-------------------------------------------------- 
		my ($tmpseq)=$sst_aln->{$acc."_".$other_taxon};
		next if !defined($tmpseq);
		print $tmpseq."\t\t".$other_taxon."\n";
	}
	print "\n\n";
	
	print "Reads:\n\n";
	if( $#all_tags >= 0){
		foreach my $tag_elem(@all_tags){
			print $tag_elem->{"tag"};
		}
	}else{
		print "No reads mapped to this miRNA.\n";
	}
	print "\n\n";
}
