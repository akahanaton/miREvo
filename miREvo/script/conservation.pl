# /usr/bin/perl -w
#===============================================================================
#         FILE:  dnaDistanceK2P.pl
#
#  DESCRIPTION:  Calculate DNA Distance using Kimura 2-Parameters Model
#        NOTES:  ---
#       AUTHOR:  Lv Yang
#      VERSION:  1.0
#      CREATED:  01/22/2010
#     REVISION:  11/06/2011
#===============================================================================

my $version = 1.00;

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SimpleAlign;
use Data::Dumper;
use MyMod::Bio::Tools::NW;

#necessary arguments
my %opts;
GetOptions(\%opts, "i=s", "o=s", "s=s", "m=s","f=s","p=s","h");
&usage if (!(defined $opts{i} and defined $opts{o} and defined $opts{s} and defined $opts{m} and defined $opts{p} and defined $opts{f} or defined $opts{h} ));

my $payoff = { match   => 3,  # $match,
			mismatch   => -3, # $mismatch,
			gap_open   => 1, # $gap_open,
			gap_extend => 1 # $gap_extend
};

my $aln_input   =   $opts{i};
my $ma_file = $opts{m};
my $hp_file = $opts{f};
my $major_taxon =   $opts{s};
my $prj_name =   $opts{p};
my $output  =   $opts{o};

my (%hp_id, $index);
my $seqin=Bio::SeqIO->new(-file=>$hp_file,-format=>"fasta");
while(my $seq=$seqin->next_seq()){
	$hp_id{lc($seq->display_id)} = ++$index;
}

my (%mature_seq, %index_info);
my $matin=Bio::SeqIO->new(-file=>$ma_file,-format=>"fasta");
my $mir_index = 1;
while(my $seq=$matin->next_seq()){
	my $new_seq = lc($seq->seq());
	$new_seq =~ tr/uU/tT/ ;
	my($mir, $new_id, $mir_fam, $tag);
	my @arr= split(/_/,$seq->display_id);
	if($arr[1] =~ m/^\d+$/){  # mirID alread index
		$new_id = join("_",$arr[0], $arr[1], $arr[3]);
		my $tmp_id = $seq->display_id;
		if(defined($hp_id{lc $tmp_id})){
			$mature_seq{lc($new_id)}{'mature'} = $new_seq;
			next;
		}else{
			$tmp_id=~ s/\*//;
			if(defined($hp_id{lc $tmp_id})){
				$mature_seq{lc($new_id)}{'star'} = $new_seq;
				next;
			}
		}
	}else{
		$mir = substr($seq->display_id,4);
		$mir_fam = $seq->display_id;
		$mir_fam =~ s/\-[35]p//;
		$mir_fam =~ s/\*$//;
	}
	if (defined($hp_id{lc $mir_fam})){
		$mir_index = $hp_id{lc $mir_fam};
	}else{
		print STDERR "Warrning: can't find corresponding hairpin sequecing for  ".$seq->display_id."\n";
	}
	if($mir =~ /\-5p/){
		$mir =~ s/\-5p//;
		$new_id = $prj_name."_".$mir_index."_".$mir;
		$tag = "5p";
	}elsif($mir =~ /\-3p/){
		$mir =~ s/\-3p//;
		$new_id = $prj_name."_".$mir_index."_".$mir;
		$tag = "3p";
	}elsif($mir =~ /\*$/){
		$mir =~ s/\-3p//;
		$new_id = $prj_name."_".$mir_index."_".$mir;
		$tag = "star";
	}else{
		$new_id = $prj_name."_".$mir_index."_".$mir;
		$tag = "mature";
	}
	$mature_seq{lc($new_id)}{lc $tag} = $new_seq;
}

#create file handle
open OUT, ">", "$output" or die "Cannot create file $output: $!.\n";

my %aln_seq = ();
my @aln_taxons;
my (%kmir_hp, %kmir_mt, $kmir_3p, %identity_hp, %identity_mt, %identity_sd);
my (%identity_info, %kmir_info);
my $alnin=Bio::AlignIO->new(-file=>$aln_input,-format=>"clustalw",-verbose => -1, -displayname_flat=>1 );
while(my $aln=$alnin->next_aln(-verbose => -1)){
	my @seqs=$aln->each_seq(-verbose => -1);
	# get reference sequence
	# dme_100_mir-969_droAna3/1727754-1727880
	my($ref_seq_hp,$refAcc);
	foreach my $seq(@seqs){
		if( $seq->display_id() =~ /$major_taxon/){
			$ref_seq_hp = $seq;
			$refAcc = $seq->display_id();
			last;
		}
	}
	my @tmpArr = split(/_/, $seqs[0]->display_id());
	my $acc = $tmpArr[0]."_".$tmpArr[1]."_".$tmpArr[2];
	next if $acc =~ s/rep\d+//;

	# hairpin conservation
	my $tmp_aln_hp = Bio::SimpleAlign->new();
	$tmp_aln_hp ->add_seq($ref_seq_hp);
	foreach my $seq(@seqs){
		@tmpArr = split(/_/, $seq->display_id());
		my ($curTaxon) = ($tmpArr[3] =~ /(^\w+\d)?/);
		if ($curTaxon ne $major_taxon) {
			push(@aln_taxons, $curTaxon);

			# kmir_info
			my($div, $p, $q) = &div_k2p($ref_seq_hp->seq(), $seq->seq());
			if(defined($div) && $div < 5){
				$kmir_info{$acc}{hairpin}{$curTaxon} = $div;
			} else {
				$kmir_info{$acc}{hairpin}{$curTaxon} = 'undef';
			}

			$tmp_aln_hp->add_seq($seq);
			my $identity = $tmp_aln_hp->percentage_identity();
			$tmp_aln_hp->remove_seq($seq);
			$identity_info{$acc}{hairpin}{$curTaxon} = $identity;
		}
	}
	 
	foreach my $tag (keys %{$mature_seq{$acc}}){
		my $cur_mature_seq = $mature_seq{$acc}{$tag};
		next if !defined($ref_seq_hp);  
		my ($mature_beg, $mature_end, $new_tag) = get_residue_number($ref_seq_hp, $cur_mature_seq);
		next if (! defined $mature_beg);
		my $column_start = $ref_seq_hp->column_from_residue_number($mature_beg);
		my $column_end = $ref_seq_hp->column_from_residue_number($mature_end);
		my $mature_aln = $aln->slice($column_start, $column_end, -verbose=>-1);
		my ($ref_seq_mt)=$mature_aln->each_seq_with_id($refAcc, -verbose => -1);
		#--------------------------------------------------
		# my ($ref_seq_mt)=$mature_aln->get_seq_by_id($refAcc, -verbose => -1);
		#-------------------------------------------------- 

		# mature conservation
		my $tmp_aln_mt = Bio::SimpleAlign->new();
		$tmp_aln_mt->add_seq($ref_seq_mt);
		my @seqs=$mature_aln->each_seq(-verbose => -1);
		foreach my $seq(@seqs){
			@tmpArr = split(/_/, $seq->display_id());
			my ($curTaxon) = ($tmpArr[3] =~ /(^\w+\d)?/);
			if ($curTaxon ne $major_taxon) {

				#--------------------------------------------------
				# print "$acc\t$new_tag\t$curTaxon\n";
				#-------------------------------------------------- 

				push(@aln_taxons, $curTaxon);

				# kmir_info
				my($div, $p, $q) = &div_k2p($ref_seq_mt->seq(), $seq->seq());
				if(defined($div) && $div < 5){
					$kmir_info{$acc}{$new_tag}{$curTaxon} = $div;
				} else {
					$kmir_info{$acc}{$new_tag}{$curTaxon} = 'undef';
				}

				$tmp_aln_mt->add_seq($seq);
				my $identity = $tmp_aln_mt->percentage_identity();
				#--------------------------------------------------
				# print Dumper($tmp_aln_mt);
				#-------------------------------------------------- 
				$tmp_aln_mt->remove_seq($seq);
				$identity_info{$acc}{$new_tag}{$curTaxon} = $identity;
			}
		}

		$column_start = $ref_seq_mt->column_from_residue_number($mature_beg+1);
		$column_end = $ref_seq_mt->column_from_residue_number($mature_beg+7);
		my $seed_aln = $mature_aln->slice($column_start, $column_end,-verbose=>-1);
		my ($ref_seq_sd)=$seed_aln->each_seq_with_id($refAcc, -verbose=>-1);
		my $tmp_aln_sd = Bio::SimpleAlign->new();
		$tmp_aln_sd->add_seq($ref_seq_sd);
		my @seqs=$seed_aln->each_seq(-verbose => -1);
		foreach my $seq(@seqs){
			@tmpArr = split(/_/, $seq->display_id());
			my ($curTaxon) = ($tmpArr[3] =~ /(^\w+\d)?/);
			if ($curTaxon ne $major_taxon) {

				#--------------------------------------------------
				# print "$acc\t$new_tag\t$curTaxon\n";
				#-------------------------------------------------- 

				push(@aln_taxons, $curTaxon);

				# kmir_info
				my($div, $p, $q) = &div_k2p($ref_seq_sd->seq(), $seq->seq());
				if(defined($div) && $div < 5){
					$kmir_info{$acc}{$new_tag." seed"}{$curTaxon} = $div;
				} else {
					$kmir_info{$acc}{$new_tag." seed"}{$curTaxon} = 'undef';
				}

				$tmp_aln_sd->add_seq($seq);
				my $identity = $tmp_aln_sd->percentage_identity();
				$tmp_aln_sd->remove_seq($seq);
				$identity_info{$acc}{$new_tag." seed"}{$curTaxon} = $identity;
			}
		}
	}

	#--------------------------------------------------
	# print join(" ", $acc,$mature_seq, $mature_beg, $mature_end,"\n");
	#-------------------------------------------------- 
}

my %seen = ();
@aln_taxons = grep { ! $seen{$_}++ } @aln_taxons;

#--------------------------------------------------
# print Dumper(%identity_info);
# print Dumper(%kmir_hp);
#-------------------------------------------------- 


print OUT "MiRNA structue\t";
foreach my $taxon ( @aln_taxons ){
	print OUT "$taxon\t";
}
print OUT "\n";

foreach my $acc (keys %kmir_info){
	print OUT "$acc\n";
	foreach my $tag ("hairpin","5p","5p seed","3p","3p seed"){
		print OUT "Kmir of $tag: \t";
		foreach my $taxon ( @aln_taxons){
			if (defined($kmir_info{$acc}{$tag}{$taxon})){
				printf OUT "%2.2f\t", $kmir_info{$acc}{$tag}{$taxon};
			}else{
				print OUT "undef\t";
			}
		}
		print OUT "\n";
	}
	foreach my $tag ("hairpin","5p","5p seed","3p","3p seed"){
		print OUT "Identity of $tag: \t";
		foreach my $taxon ( @aln_taxons){
			if (defined($identity_info{$acc}{$tag}{$taxon})){
				printf OUT "%2.2f\t", $identity_info{$acc}{$tag}{$taxon};
			}else{
				print OUT "undef\t";
			}
		}
		print OUT "\n";
	}
	print OUT "\n";
}


#close file handle
close OUT;
	
#========================================SUB==============================================
sub div_k2p {
    my($seq1, $seq2) = @_;
    my ($transition, $transversion, $n, $k, $p, $q) = (0, 0, 0, 0, 0, 0);
    for my $i (0..length($seq1)-1) {
        my $base1 = substr($seq1, $i, 1);
        my $base2 = substr($seq2, $i, 1);
        next if($base1 eq "-" or $base2 eq "-" or $base1 eq "N" or $base2 eq "N");
        $n++;
        $base1 =~ tr/atcgu/ATCGT/;
        $base1 =~ tr/atcgu/ATCGT/;
        next if($base1 eq $base2); 
        if($base1 eq "A" and $base2 eq "G") {
            $transition++;
        } 
        elsif($base1 eq "T" and $base2 eq "C") {
            $transition++;
        } 
        elsif($base1 eq "C" and $base2 eq "T") {
            $transition++;
        } 
        elsif($base1 eq "G" and $base2 eq "A") {
            $transition++;
        }
        else {
            $transversion++;
        }
    }

    #1) Two sequences are not available
    if(!$n) {
        return (0, 0, 0);
    }

    #2) K2P distance
    $p = $transition / $n;
    $q = $transversion / $n;
    if(1-2*$p-$q > 0 and 1-2*$q >0 ) {
        $k = (-1/2)*log((1-2*$p-$q)*sqrt(1-2*$q));
        return (abs($k), $p, $q);
    } else {

    #3) Two sequences are far diverged and cannot calculate K2P distance
        return (5, -1, -1);
    }
}

#========================================SUB==============================================

sub usage{
    print <<"USAGE";
Usage:
	$0 -i <input file> -s <reference species name> -o <output file>
options:
	-i  input file (clustalw alighment file)
	-m  fasta file with query miRNAs, required
	-b  fasta file with MATURE sequence of miRBase annotated miRNAs in related speceis, optional
	-s  reference species name in your MAF file (e.g. dm3)
	-o  output file
	-h  help
USAGE
exit(1);
}

sub rev_com_seq {
	my $seq = shift ;
	$seq =~ tr/atcgnATCGN/tagcnTAGCN/ ;
	$seq =reverse($seq);
	return($seq);
}

sub get_residue_number{
	my ($raw_hairpin_seq, $mat_seq) = @_;
	my $tag;
	my $hairpin_seq = lc($raw_hairpin_seq->seq());
	$hairpin_seq =~ s/\-//g;
	my $hairpin_lng = length $hairpin_seq;
	my ($mature_beg, $mature_end);
	my $tmp_pos;
	$tmp_pos = index($hairpin_seq, $mat_seq);
	if ($tmp_pos < 0){
		my $mat_seq_rec = rev_com_seq($mat_seq);
		$tmp_pos = index($hairpin_seq, $mat_seq);
		if ($tmp_pos < 0){
			print STDERR "here\t$mat_seq\t$hairpin_seq\n";
			return (undef, undef, undef);
		}
	}
	if( $tmp_pos> ($hairpin_lng / 2)){
		$tag = '3p';
	}else{
		$tag = '5p';
	}

	$mature_beg = $tmp_pos + $raw_hairpin_seq->start;
	$mature_end = $mature_beg + length $mat_seq;
	return ($mature_beg, $mature_end, $tag);
}

#========================================END==============================================
			#--------------------------------------------------
			# my $nw = new Align::NW $ref_seq_hp->seq, $seq->seq(), $payoff;
			# $nw->score;
			# $nw->align;
			# my $identity = $nw->get_identity;
			#-------------------------------------------------- 
