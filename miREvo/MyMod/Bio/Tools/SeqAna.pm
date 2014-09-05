#
#===============================================================================
#
#         FILE:  SeqAna.pm
#
#  DESCRIPTION:  Seqence analysis modules
#
#        FILES:  MyMod/Bio/Tools/SeqAna.pm
#         BUGS:  
#        NOTES:  
#       AUTHOR:  Shen Yang
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  2006年03月23日 14时29分42秒 CST
#     REVISION:  
#===============================================================================

package MyMod::Bio::Tools::SeqAna;
use vars qw(@ISA $Default_Source);
use strict;
use warnings;

use Bio::Root::Root; 
use Bio::Root::IO;
use MyMod::MyRoot;

@ISA= qw(Bio::Root::Root Bio::Root::IO MyMod::MyRoot);
sub new {
    my ($caller,@args)=@_;
    my $class=ref($caller)||$caller;
    my $self=$class->SUPER::new(@args);
    return $self;
}

sub concat_aln{
	my ($self,$aln1,$aln2,$noorder,$same_acc)=@_;
	my $new_aln=Bio::SimpleAlign->new(); 
	my @seqary1=$aln1->each_seq;
	my @seqary2=$aln2->each_seq;
	my $all_accs;
	grep{$all_accs->{$_->display_id()}++}(@seqary1,@seqary2); 
	my @accs;
	foreach my $acc(keys(%$all_accs)){
		push(@accs,$acc) if !defined($same_acc)||$all_accs->{$acc}==2;
		print $acc." only one sequence\n" if $all_accs->{$acc}==1;
	} 
	my $add_seq_id;
	my $max_add_len=0;
	foreach my $acc(@accs){ 
		my ($seq1)=$aln1->each_seq_with_id($acc);
		my ($seq2)=$aln2->each_seq_with_id($acc); 
		next if !defined($seq1)||!defined($seq2); 
		#next if !defined($seq1)&&!defined($seq2); 
		if(defined($noorder)||$seq1->end()+1 ==$seq2->start()||$seq1->start()-1==$seq2->end()){
		#if(defined($seq1)&&defined($seq2)){
			if(defined($noorder)||abs($seq1->end()-$seq2->start())<500||abs($seq1->start()-$seq2->end())<500){
			}
			else{
				next;
			}
			if($seq1->end()+1 ==$seq2->start()||$seq1->start()-1==$seq2->end()){
			}
			else{
				my $add_len=abs($seq1->end()-$seq2->start())+1<abs($seq1->start()-$seq2->end())+1?abs($seq1->end()-$seq2->start())+1:abs($seq1->start()-$seq2->end())+1;
				$add_seq_id->{$acc}=$add_len;
				$max_add_len=$add_len if $max_add_len<$add_len;
			}
		}
		
			
 		my $seqstr1=defined($seq1)?$seq1->seq():"-"x$aln1->length();
		my $seqstr2=defined($seq2)?$seq2->seq():"-"x$aln2->length(); 
		my $seq_str=$seqstr1.$seqstr2;
		my $seq1_start=defined($seq1)?$seq1->start():$seq2->start();
		my $seq2_start=defined($seq2)?$seq2->start():$seq1->start(); 
		my $seq1_end=defined($seq1)?$seq1->end():$seq2->end();
		my $seq2_end=defined($seq2)?$seq2->end():$seq1->end();
		my $min_start=$seq1_start<$seq2_start?$seq1_start:$seq2_start;
		my $max_end=$seq1_end>$seq2_end?$seq1_end:$seq2_end;
		my $new_seq=Bio::LocatableSeq->new(
			-start=>$min_start,
			-end=>$max_end,
			-strand=>$seq1->strand(),
			-display_id=>$acc,
			-verbose=>-1,
			-seq=>$seq_str
		);
		$new_aln->add_seq($new_seq);
	}
	#foreach my $acc(@accs){
	#	my $seq=$new_aln->each_seq_with_id($acc);
	#	if(defined($add_seq_id->{$acc})){
	#		my $seqstr="N"x$add_seq_id->{$acc}."-"x$max_add_len-$add_seq_id->{$acc};
	#		
	return $new_aln;
}

sub concat_aln_for_mulaln{
	my ($self,$aln1,$aln2,$no_insert)=@_;
	my $new_aln=Bio::SimpleAlign->new(); 
	my @seqary1=$aln1->each_seq;
	my @seqary2=$aln2->each_seq;
	my $all_accs;
	grep{$all_accs->{$_->display_id()}++}(@seqary1,@seqary2); 
	my @accs=keys(%$all_accs);
	my $add_seq_id;
	my $max_add_len=0;
	foreach my $acc(@accs){ 
		my ($seq1)=$aln1->each_seq_with_id($acc);
		my ($seq2)=$aln2->each_seq_with_id($acc); 
		next if defined($no_insert) && (!defined($seq1)||!defined($seq2)); 
		if(!defined($seq1)&&defined($seq2)){
			my $seq_str="-"x$aln1->length().$seq2->seq();
			my $new_seq=Bio::LocatableSeq->new(
				-start=>$seq2->start(),
				-end=>$seq2->end(),
				-strand=>$seq2->strand(),
				-display_id=>$acc,
				-verbose=>-1,
				-seq=>$seq_str
			);
			$new_aln->add_seq($new_seq);
		} 
		elsif(!defined($seq2)&&defined($seq1)){
			my $seq_str=$seq1->seq()."-"x$aln2->length();
			my $new_seq=Bio::LocatableSeq->new(
				-start=>$seq1->start(),
				-end=>$seq1->end(),
				-strand=>$seq1->strand(),
				-display_id=>$acc,
				-verbose=>-1,
				-seq=>$seq_str
			);
			$new_aln->add_seq($new_seq);
		}
		else{
			if($seq1->strand()==$seq2->strand()&&($seq1->end()+1 ==$seq2->start()||$seq1->start()-1==$seq2->end())){
				my $seqstr1=$seq1->seq();
				my $seqstr2=$seq2->seq(); 
				my $seq_str=$seqstr1.$seqstr2;
				my $seq1_start=$seq1->start();
				my $seq2_start=$seq2->start(); 
				my $seq1_end=$seq1->end();
				my $seq2_end=$seq2->end();
				my $min_start=$seq1_start<$seq2_start?$seq1_start:$seq2_start;
				my $max_end=$seq1_end>$seq2_end?$seq1_end:$seq2_end;
				my $new_seq=Bio::LocatableSeq->new(
					-start=>$min_start,
					-end=>$max_end,
					-strand=>$seq1->strand(),
					-display_id=>$acc,
					-verbose=>-1,
					-seq=>$seq_str
				);
				$new_aln->add_seq($new_seq);
			} 
			else{
				next;
			}
		}	
	}#
	return $new_aln;
}

sub concat_aln_for_mulaln_all{
	my ($self,$aln1,$aln2,$pos_str)=@_;
	my $new_aln=Bio::SimpleAlign->new(); 
	my @seqary1=$aln1->each_seq;
	my @seqary2=$aln2->each_seq;
	my $all_accs;
	grep{$all_accs->{$_->display_id()}++}(@seqary1,@seqary2); 
	my @accs=keys(%$all_accs);
	my $add_seq_id;
	my $max_add_len=0;
	foreach my $acc(@accs){ 
		my ($seq1)=$aln1->each_seq_with_id($acc);
		my ($seq2)=$aln2->each_seq_with_id($acc); 
		if(!defined($seq1)&&defined($seq2)){
			my $seq_str="-"x$aln1->length().$seq2->seq();
			my $new_seq=Bio::LocatableSeq->new(
				-start=>$seq2->start(),
				-end=>$seq2->end(),
				-strand=>$seq2->strand(),
				-display_id=>$acc,
				-verbose=>-1,
				-seq=>$seq_str
			);
			$new_aln->add_seq($new_seq);
		} 
		elsif(!defined($seq2)&&defined($seq1)){
			my $seq_str=$seq1->seq()."-"x$aln2->length();
			my $new_seq=Bio::LocatableSeq->new(
				-start=>$seq1->start(),
				-end=>$seq1->end(),
				-strand=>$seq1->strand(),
				-display_id=>$acc,
				-verbose=>-1,
				-seq=>$seq_str
			);
			$new_aln->add_seq($new_seq);
		}
		else{
			my $pos_str1=$seq1->start()."_".$seq1->end()."_".$seq1->strand();
			my $pos_str2=$seq2->start()."_".$seq2->end()."_".$seq2->strand();
			push(@{$pos_str->{$acc}},{"seq1"=>$pos_str1,"seq2"=>$pos_str2});
			my $seqstr1=$seq1->seq();
			my $seqstr2=$seq2->seq(); 
			my $seq_str=$seqstr1.$seqstr2;
			my $seq1_start=$seq1->start();
			my $seq2_start=$seq2->start(); 
			my $seq1_end=$seq1->end();
			my $seq2_end=$seq2->end();
			my $min_start=$seq1_start<$seq2_start?$seq1_start:$seq2_start;
			my $max_end=$seq1_end>$seq2_end?$seq1_end:$seq2_end;
			my $new_seq=Bio::LocatableSeq->new(
				-start=>$min_start,
				-end=>$max_end,
				-strand=>$seq1->strand(),
				-display_id=>$acc,
				-verbose=>-1,
				-seq=>$seq_str
			);
			$new_aln->add_seq($new_seq);
		} 
	}
	return ($new_aln,$pos_str);
}

sub concat_aln_ary_mulaln{
	my ($self,$aln_ary,$noorder,$same_acc)=@_;
	my $all_aln=$aln_ary->[0];
	for(my $i=1;$i<scalar(@$aln_ary);$i++){
		$all_aln=$self->concat_aln_for_mulaln($all_aln,$aln_ary->[$i],$noorder,$same_acc);
	} 
	return $all_aln;
} 

sub concat_aln_ary_mulaln_all{
	my ($self,$aln_ary)=@_;
	my $all_aln=$aln_ary->[0];
	my $pos_str;
	for(my $i=1;$i<scalar(@$aln_ary);$i++){
		($all_aln,$pos_str)=$self->concat_aln_for_mulaln_all($all_aln,$aln_ary->[$i],$pos_str);
	} 
	return ($all_aln,$pos_str);
} 

sub concat_aln_ary{
	my ($self,$aln_ary,$noorder,$same_acc)=@_;
	my $all_aln=$aln_ary->[0];
	for(my $i=1;$i<scalar(@$aln_ary);$i++){
		$all_aln=$self->concat_aln($all_aln,$aln_ary->[$i],$noorder,$same_acc);
	} 
	return $all_aln;
} 

sub revcom_aln{
	use Bio::SimpleAlign;
	my ($self,$aln)=@_;
	my $new_aln=Bio::SimpleAlign->new(-verbose=>-1);
	my @seqs=$aln->each_seq();
	foreach my $seq(@seqs){
		$seq=$seq->revcom();
		$new_aln->add_seq($seq);
	}
	return $new_aln;
}

sub aln_identity_seq{
	my ($self,$seq1,$seq2)=@_;
	my $seq1str=$seq1->seq();
	my $seq2str=$seq2->seq();
	$seq1str=~ s/-//sg;
	$seq2str=~ s/-//sg;
	my $loc=index($seq1str,$seq2str);
	return undef if $loc==-1;
	my $new_seq2str="-"x$loc.$seq2str."-"x(length($seq1str)-$loc+length($seq2str));
	my @origin_chars=split(//,$seq1->seq());
	my @seq1_chars=split(//,$seq1str);
	my @seq2_chars=split(//,$new_seq2str);
	my @new_seq2_chars;
	foreach my $char(@origin_chars){
		if($char eq "-"){
			push(@new_seq2_chars,$char);
		}
		else{
			push(@new_seq2_chars,shift(@seq2_chars));
		}
	}
	my $new_seq2=Bio::LocatableSeq->new(
		-display_id=>$seq2->display_id(),
		-verbose=>-1,
		-seq=>join("",@new_seq2_chars),
		-start=>1,
		-end=>$seq2->length()
	);
	return ($new_seq2,$loc);
}

sub revcom{
	my ($self,$sequence) = @_;

	my $res = reverse($sequence);

	$res =~ tr/ACGT/TGCA/;

	return $res;
}

1;

