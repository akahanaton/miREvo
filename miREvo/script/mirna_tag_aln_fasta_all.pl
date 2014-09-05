#!/user/bin/perl -w
use strict;
use warnings;
use MyMod::Bio::Tools::miRNA;
use MyMod::Bio::Tools::SeqAna;
use Bio::SeqIO;
use Data::Dumper;

my $seqana=MyMod::Bio::Tools::SeqAna->new();
my $mirna_fact=MyMod::Bio::Tools::miRNA->new();


my $anno_info;
if(defined($ARGV[2])){
	open(ANNO,$ARGV[2]);
	while(my $line=<ANNO>){
		chomp($line);
		my @ary=split(/\t/,$line);
		my $acc=$ary[0];
		$anno_info->{$acc}=$ary[1];
	}
}


#read tags information

my $mapping_result=$ARGV[1]; #bowtie mapping results
# dme_0_x3626256  +       dme-mir-1       55      TGGAATGTAAAGAAGTATGGAG  IIIIIIIIIIIIIIIIIIIIII  0
my $tags_info;
open(TAG,$mapping_result);
while(my $line=<TAG>){
	chomp($line);
	my @ary=split(/\t/,$line);
	my ($lib,$mir,$copy)=split(/_/,$ary[0]);
	$copy =~ s/x//;
	my $elem={"copy"=>$copy,"seq"=>$ary[4],"acc"=>$lib,"strand"=>$ary[1]};
	#if($ary[8]=~ m/(.*)_(.*)/sg){
	#	$ary[8]=$1;
	#}
	push(@{$tags_info->{$ary[2]}},$elem);
}


my $align_file=$ARGV[0];
my $seqin=Bio::SeqIO->new(-file=>$align_file,-format=>"fasta");

while(my $seq=$seqin->next_seq()){	
	my $acc=$seq->display_id();
	my $aln=Bio::SimpleAlign->new();
	my $tmplocseq=Bio::LocatableSeq->new(
		-display_id=>$seq->display_id(),
		-seq=>$seq->seq(),
		-start=>1,
		-end=>$seq->length()
	);
	$aln->add_seq($tmplocseq);
	my ($sst_aln)=$mirna_fact->get_sst_aln($aln,1);
	my @all_tags;
	my ($tmpseq) =$aln->each_seq_with_id($acc);
	if (defined($tags_info->{$acc})){
		my @tags=@{$tags_info->{$seq->display_id()}};
		my @tmp_all_tags;
		foreach my $elem(@tags){
			my $tmptag=Bio::LocatableSeq->new(
				-display_id=>"tag",
				-seq=>$elem->{"seq"},
				-start=>1,
				-end=>length($elem->{"seq"}),
				-strand=>1
			);
			my ($aln_seq,$tagstart)=$seqana->aln_identity_seq($tmpseq,$tmptag);
			next if !defined($aln_seq);
			my $aln_seq_str=$aln_seq->seq();
			$aln_seq_str=~ s/T/U/sg;
			$aln_seq_str=lc($aln_seq_str) if $elem->{"strand"} eq "-";
			#throw out negative strand
			#next if $elem->{"strand"} eq "-";
			my $tag_line=$aln_seq_str." ".$elem->{"acc"}." ".$elem->{"copy"}." ".$elem->{"strand"}."\n";
			push(@tmp_all_tags,{"tag"=>$tag_line,"pos"=>$tagstart});
		}
		my @sort_tags=sort{$a->{"pos"}<=>$b->{"pos"}} @tmp_all_tags;
		grep{push(@all_tags,$_)}@sort_tags;
		next if scalar(@sort_tags)<1;
	}

	print ">$acc\n";
	
	if(defined($anno_info->{$acc})){
		#my @rep_taxons=split(/__/,$rep_info->{$acc});
		$anno_info->{$acc}=~ s/__/; /sg;
		print "Another annotations: ".$anno_info->{$acc}."\n\n";
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
		print "No reads mapped to this miRNA.\n";
	}
	
	print "\n";
	
}
