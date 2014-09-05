use strict;
use warnings;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;

#extract pre-miRNA,calculate the RNAforester score and so on

my $alnin=Bio::AlignIO->new(-file=>$ARGV[0],-format=>"clustalw");
my $alnout=Bio::AlignIO->new(-file=>">$ARGV[1]",-format=>"clustalw");
$alnout->line_length(500);
my $seqout=Bio::SeqIO->new(-file=>">$ARGV[2]",-format=>"fasta");

my $skip_flag=0;	
while(my $aln=$alnin->next_aln()){	
	my $final_aln=$aln;
	my $out_aln=Bio::SimpleAlign->new();
	foreach my $seq($final_aln->each_seq()){
		my ($acc,$taxon)=($seq->display_id()=~ m/(.*)_(.*)\.(.*)/sg);
		#--------------------------------------------------
		# next if($taxon_used->{$taxon}==0);
		#-------------------------------------------------- 
		#--------------------------------------------------
		# $seq->display_id($acc."_".$taxon);	
		#-------------------------------------------------- 
		my $seqstr=$seq->seq();
		$seqstr=~ s/-//sg;
		my $tmpseq=Bio::PrimarySeq->new(
			-display_id=>$seq->display_id(),
			-seq=>$seqstr
		);
		$seqout->write_seq($tmpseq);
		#--------------------------------------------------
		# $seq->display_id($acc."_".$taxon);
		#-------------------------------------------------- 
		$out_aln->add_seq($seq);
	}
	$out_aln=$out_aln->remove_gaps("-",1);
	$alnout->write_aln($out_aln);
}

