package MyMod::Bio::Tools::miRNA;
use vars qw(@ISA $Default_Source);
use strict;
use warnings;

use Bio::Root::Root; 
use Bio::Root::IO;
use MyMod::MyRoot;
use Bio::SimpleAlign;
use Bio::Symbol::Symbol;
use Bio::Symbol::Alphabet;
use MyMod::Bio::Tools::RNA;

@ISA= qw(Bio::Root::Root Bio::Root::IO MyMod::MyRoot);
sub new {
    my ($caller,@args)=@_;
    my $class=ref($caller)||$caller;
    my $self=$class->SUPER::new(@args);
    return $self;
}

sub get_sst_aln{
	my ($self,$aln,$ifenergy)=@_;
	$aln->map_chars('\.','-');
	my @seqs=$aln->each_seq();
	my $sst_aln=Bio::SimpleAlign->new();
	my $sst_hash;
	my @ids;
	foreach my $seq(@seqs){
		my $origin_seq=$seq->seq();
		push(@ids,$seq->display_id());
		push(@ids,$seq->display_id()."_sst");
		my $seqstr=$seq->seq();
		$seqstr=~ s/[-\.]//sg;
		my ($sst,$mfe)=RNA::fold($seqstr);;
		$mfe=int($mfe*100)/100;
		my @chars=split(//,$origin_seq);
		my @sst_chars=split(//,$sst);
		my @new_sst_chars;
		foreach my $char(@chars){
			if($char eq "-" || $char eq "."){
				push(@new_sst_chars,$char);
			}
			else{
				push(@new_sst_chars,shift(@sst_chars));
			}
		}
		my $new_sst_str=join('',@new_sst_chars);
		$new_sst_str=$new_sst_str."($mfe)" if defined($ifenergy);
		#--------------------------------------------------
		# my $sst_seq=Bio::LocatableSeq->new(
		# 	-display_id=>$seq->display_id()."_sst",
		# 	-seq=>$new_sst_str,
		# 	-start=>$seq->start(),
		# 	-end=>$seq->end(),
		# 	-strand=>$seq->strand(),
		# 	-alphabet=>$alphabet
		# );
		# $sst_aln->add_seq($sst_seq);
		#-------------------------------------------------- 
		$sst_hash->{$seq->display_id()} = $new_sst_str;
	}
	#--------------------------------------------------
	# return ($sst_aln,@ids);
	#-------------------------------------------------- 
	return ($sst_hash,@ids);
}
		
1;
