use strict;
use warnings;
use Bio::SeqIO;
use Bio::AlignIO;
use MyMod::Bio::Tools::SeqAna;
#--------------------------------------------------
# use Data::Dumper;
#-------------------------------------------------- 

my $seqana=MyMod::Bio::Tools::SeqAna->new();


my $aln_dir = $ARGV[0];
my $gff_file = $ARGV[1];
my $acc = $ARGV[2];


open (GFF, $gff_file)  || die("** Fail to open $gff_file \n");
while (my $line = <GFF>){
chomp $line;
my @tmpArr = split(/\s+/, $line);
my $maf_file=$aln_dir."/".$tmpArr[8].".maf";
my $start = $tmpArr[3];
my $end = $tmpArr[4];
my $strand = $tmpArr[6];
my $mirAcc = $tmpArr[8];
next if ( ! -e $maf_file );
my $alnin=Bio::AlignIO->new(-file=>$maf_file,-format=>"maf", -verbose=>-1); #the multiple alignment with maf format
my $alnout=Bio::AlignIO->new(-file=>">".$maf_file."-aln",-format=>"clustalw", -verbose=>-1);

my $aln_ary;
my $done=0;

my $major_taxon=$ARGV[2]||"dm3";

while(my $aln=$alnin->next_aln()){
	push(@{$aln_ary},$aln);
}

next if (!defined($aln_ary));

my @seqs=$aln_ary->[0]->each_seq();
my $major_id;
foreach my $seq(@seqs){
	$major_id=$seq->display_id() if $seq->display_id()=~ m/$major_taxon/sg;
}
my ($major_seq)=$aln_ary->[0]->each_seq_with_id($major_id);
my ($taxon,$chr_num)=split(/\./,$major_id);
my $aln_start=$major_seq->start();
my ($final_aln,$pos_str)=$seqana->concat_aln_ary_mulaln_all($aln_ary);

#--------------------------------------------------
# print Dumper($final_aln) if $maf_file =~ /2b-1/;;
#-------------------------------------------------- 
		
my $tmpaln=Bio::SimpleAlign->new();
my %aln_seq;
foreach my $seq($final_aln->each_seq){
	my ($taxon,$chr)=($seq->display_id()=~ m/(.*)\.(.*)/sg);
	if(!defined($aln_seq{$taxon})){
		$aln_seq{$taxon}{"Seq"}= $seq->seq();
		$seq->display_id($mirAcc."_".$seq->display_id());
		$aln_seq{$taxon}{"AlnSeq"}= $seq;
	}else{
		my $cur_seq = $seq->seq();
		$cur_seq =~ s/-//g;
		my $old_seq = $aln_seq{$taxon}{"Seq"};
		$old_seq =~ s/-//g;
		if( $cur_seq ne $old_seq && length($cur_seq) >= length($old_seq)){
				$aln_seq{$taxon}{"Seq"}= $seq->seq();
				$seq->display_id($mirAcc."_".$seq->display_id());
				$aln_seq{$taxon}{"AlnSeq"}= $seq;
		}
	}
}
#--------------------------------------------------
# print Dumper(%aln_seq) if $maf_file =~ /2b-1/;;
#-------------------------------------------------- 
foreach my $taxon (keys %aln_seq){
	$tmpaln->add_seq($aln_seq{$taxon}{"AlnSeq"});
}

#--------------------------------------------------
# print Dumper($tmpaln) if $maf_file =~ /2b-1/;;
#-------------------------------------------------- 

$final_aln=$tmpaln;
if($strand eq '-'){
	$final_aln=$seqana->revcom_aln($final_aln);
}
my ($final_major_seq)=$final_aln->each_seq_with_id($mirAcc."_".$major_id);
if(!defined($final_major_seq)){
	print $mirAcc."not exist\n";
	next;
}
$alnout->write_aln($final_aln);
my $seqstr=$final_major_seq->seq();
$seqstr=~ s/-//sg;
$final_major_seq->seq($seqstr);
$final_major_seq->display_id($mirAcc."".$major_id);
#--------------------------------------------------
# if(abs($final_major_seq->start()-$start)>10 || abs($final_major_seq->end()-$end)>10){
# 	print STDERR $mirAcc." not perfect"."\n";		
# }
#-------------------------------------------------- 

#$seqout->write_seq($final_major_seq);
$aln_ary=undef;
$done=1;
}

#--------------------------------------------------
# foreach my $acc(keys(%$pos_str)){
# 	my @ary=@{$pos_str->{$acc}};
# 	my $indel_info;
# 	for(my $i=0;$i<scalar(@ary);$i++){
# 		my ($start1,$end1,$strand1)=split(/_/,$ary[$i]->{"seq1"});
# 		my ($start2,$end2,$strand2)=split(/_/,$ary[$i]->{"seq2"});
# 		my $indel;
# 		my $start;
# 		my $end;
# 		if($strand1 eq "+"){
# 			$start=$start2;
# 			$end=$end1;
# 		}
# 		else{
# 			$start=$start1;
# 			$end=$end2;
# 		}
# 		$indel=$start-$end-1;
# 		if($indel>0){
# 			push(@$indel_info,$indel."_".$start."_".$end);
# 		}
# 	}
#   print $acc."\t".join("\t",@$indel_info)."\n" if defined($indel_info);
# }
#-------------------------------------------------- 
