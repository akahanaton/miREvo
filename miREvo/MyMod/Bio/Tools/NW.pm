package Align::NW;

use 5.005;
use strict;
use integer;
use vars qw($VERSION @ISA);

require Exporter;

@ISA = qw(Exporter);
$VERSION = '1.01';

$| = 1;

sub new
{
    my($package, $a, $b, $payoff) = @_;

    $a = ['', split //, $a];
    $b = ['', split //, $b];

    my $rows = @$a;
    my $cols = @$b;

    my $dp;
    for (my $row=0; $row<$rows; $row++)
    {
	for (my $col=0; $col<$cols; $col++)
	{
	    my $cell = { row   => $row,
			 col   => $col,
		         score => 0   };

	    $dp->[$row][$col] = $cell;
	}
    }

    my $nw = { a       => $a,
	       b       => $b,
	       rows    => $rows,
	       cols    => $cols,
	       dp      => $dp,
	       payoff  => $payoff};

    bless $nw, $package
}


sub score
{
    my $nw      = shift;
    my $dp      = $nw->{dp};
    my $a       = $nw->{a};
    my $b       = $nw->{b};

    my $rows = @$a;
    my $cols = @$b;

    my $payoff     = $nw->{payoff};
    my $match      = $payoff->{match};
    my $mismatch   = $payoff->{mismatch};
    my $gap_open   = $payoff->{gap_open};
    my $gap_extend = $payoff->{gap_extend};


    for (my $row=1; $row<$rows; $row++)
    {
	my $a1 = $a->[$row];
	for (my $col=1; $col<$cols; $col++)
	{
	    my $cell = $dp->[$row][$col];

	    my $b1 = $b->[$col];
	    my $compare = $a1 eq $b1 ? $match : $mismatch;
	    my $prev = $dp->[$row-1][$col-1];
	    $cell->{score} = $prev->{score} + $compare;
	    $cell->{prev}  = $prev;

	    for (my $r=0; $r<$row; $r++)
	    {
		my $prev  = $dp->[$r][$col];
		my $score = $prev->{score} + $gap_open + $gap_extend * ($row-$r);
		$score < $cell->{score} and next;
		$cell->{score} = $score;
		$cell->{prev}  = $prev;
	    }

	    for (my $c=0; $c<$col; $c++)
	    {
		my $prev  = $dp->[$row][$c];
		my $score = $prev->{score} + $gap_open + $gap_extend * ($col-$c);
		$score < $cell->{score} and next;
		$cell->{score} = $score;
		$cell->{prev}  = $prev;
	    }
	}
    }

}


sub dump_score
{
    my $nw = shift;
    my $a  = $nw->{a};
    my $b  = $nw->{b};
    my $dp = $nw->{dp};

    my @b = join('   ', @$b);
    print " @b\n";

    my $rows = @$a;
    my $cols = @$b;

    for (my $row=1; $row<$rows; $row++)
    {
	my $a1 = $a ->[$row];
	my $r1 = $dp->[$row];
	my @s1 = map { sprintf "%3d", $_->{score} } @$r1;
	shift @s1;
	print "$a1 @s1\n";
    }

    print "\n";
}


sub align
{
    my $nw = shift;

    $nw->{align} = { a => [],
		     s => [],
		     b => [] };

    my $cell = $nw->_max_cell;
               $nw->{score} = $cell->{score};

	       $nw->_align_tail($cell);
       $cell = $nw->_align_body($cell);
    	       $nw->_align_head($cell);

	       $nw->_join_align;
}


sub get_score
{
    my $nw = shift;
    $nw->{score}
}


sub _max(&@)
{
    my $less = shift;
    my $max  = shift;

    &$less($max, $_) and $max = $_ for (@_);

    $max
}


sub _max_cell
{
    my $nw  = shift;
    my $dp  = $nw->{dp};

    my $rows = $nw->{rows};
    my $cols = $nw->{cols};

    my @right  = map { $dp->[$_][-1] } 1..$rows-1;
    my @bottom = map { $dp->[-1][$_] } 1..$cols-1;

    _max { $_[0]->{score} < $_[1]->{score} } @right, @bottom
}


sub _align_tail
{
    my($nw, $cell) = @_;

    my $a = $nw->{a};
    my $b = $nw->{b};

    my $row = $cell->{row};
    my $col = $cell->{col};

    my @a = @$a[$row+1..$#$a];
    my @s = ();
    my @b = @$b[$col+1..$#$b];
    
    $nw->_unshift_align(\@a, \@s, \@b);
}


sub _align_body
{
    my($nw, $cell) = @_;

    my $dp = $nw->{dp};
    my $a  = $nw->{a};
    my $b  = $nw->{b};
 
    my(@a, @s, @b);

    for (;;)
    {	
	my $row = $cell->{row};
	my $col = $cell->{col};
	$row and $col or last;

	my $prev = $cell->{prev};

	if ($prev->{row} < $row and $prev->{col} < $col)
	{
	    my $a1 = $a->[$row];
	    my $b1 = $b->[$col];
	    unshift @a, $a1;
	    unshift @s, $a1 eq $b1 ? '|' : ' ';
	    unshift @b, $b1;
	}
	elsif ($prev->{row} < $row)
	{
	    my $gap = $row - $prev->{row};
	    unshift @a, @$a[$row-$gap+1..$row];
	    unshift @s, ' ' x $gap;
	    unshift @b, '.' x $gap;
	}
	else
	{
	    my $gap = $col - $prev->{col};
	    unshift @a, '.' x $gap;
	    unshift @s, ' ' x $gap;
	    unshift @b, @$b[$col-$gap+1..$col];
	}

	$cell = $prev;
    }

    $nw->_unshift_align(\@a, \@s, \@b);
    $cell
}


sub _align_head
{
    my($nw, $cell) = @_;

    my $a  = $nw->{a};
    my $b  = $nw->{b};

    my $row = $cell->{row};
    my $col = $cell->{col};
    my $max = _max { $_[0] < $_[1] } $row, $col;

    my @a = (' ' x $col, @$a[1..$row]);
    my @s = (' ' x $max		     );
    my @b = (' ' x $row, @$b[1..$col]);

    $nw->_unshift_align(\@a, \@s, \@b);
}


sub _unshift_align
{
    my($nw, $a, $s, $b) = @_;
    my $align = $nw->{align};

    unshift @{$align->{a}}, @$a;
    unshift @{$align->{s}}, @$s;
    unshift @{$align->{b}}, @$b;
}


sub _join_align
{
    my $nw = shift;
    my $align = $nw->{align};

    for my $key (keys %$align)
    {
	my $x = $nw->{align}{$key};
	$nw->{align}{$key} = join('', @$x);
    }
}    

sub seed_search{
	my $nw = shift;
	my $a = join("",@{$nw->{a}});
	my $b = join("",@{$nw->{b}});
	$a =~ s/\s+//g;
	$b =~ s/\s+//g;
	my $a_len = length $a;
	my $b_len = length $b;
	my $seed_size = 7;


	my $cur_seed_a = substr( $a, 1, $seed_size);
	my $cur_seed_b = substr( $b, 1, $seed_size);
	return ($cur_seed_a eq $cur_seed_b);
}

sub get_align
{
    my $nw = shift;
    $nw->{align}
}

sub get_match_line{
    my $nw = shift;
    return $nw->{align}->{s};
}

sub get_identity{
	# the first sequence is  the reference
    my $nw = shift;
    my $ref_seq = $nw->{align}->{a};
    my $match_str =$nw->{align}->{s};
	$ref_seq =~ s/[\s+\.\-]//g;
	$match_str =~ s/\s+//g;
	#--------------------------------------------------
	# print "$ref_seq\t$match_str\n";
	#-------------------------------------------------- 
	my $ref_len =length $ref_seq;
	my $match_len =length $match_str;
	if($ref_len>0){
		return (100*$match_len/$ref_len);
	}else{
		return 0;
	}
}

sub print_align
{
    my $nw = shift;
    my $align = $nw->{align};
    my $a = $align->{a};
    my $s = $align->{s};
    my $b = $align->{b};
    my $lineLen = 60;

    $a =~ tr[ -~][^]c;
    $b =~ tr[ -~][^]c;
	my $i;
    for ($i=0; $i<(length($a)-$lineLen); $i+=$lineLen)
    {
		print substr($a, $i, $lineLen), "\n";
		print substr($s, $i, $lineLen), "\n";
		print substr($b, $i, $lineLen), "\n";
		print "\n";
    }
	print substr($a, $i), "\n";
	print substr($s, $i), "\n" if($i < length($s));
	print substr($b, $i), "\n" if($i < length($b));
	print "\n";
}

1

__END__

=head1 NAME

Align::NW - Needleman-Wunsch algorithm for optimal global sequence alignment

=head1 SYNOPSIS

    use Align::NW;
  
    $payoff = { match      => $match,
		mismatch   => $mismatch,
		gap_open   => $gap_open,
		gap_extend => $gap_extend };
  
    $nw = new Align::NW $a, $b, $payoff;
    $nw->score;
    $nw->align;
  
    $score = $nw->get_score;
    $align = $nw->get_align;

    $nw->print_align;
    $nw->dump_score;

=head1 DESCRIPTION

C<Align::NW> finds the optimal global alignment of the sequences
C<$a> and C<$b>, subject to the C<$payoff> matrix.

=head2 Algorithm

C<Align::NW> uses the Needleman-Wunsch dynamic programming algorithm.
This algorithm runs in O(a*b*(a+b)), where a and b are the 
lengths of the two sequences to be aligned. 

=head2 Alignments

An alignment of two sequences is represented by three lines.
The first line shows the first sequence,
and the third line shows the second sequence.

The second line has a row of symbols.
The symbol is a vertical bar where ever characters in the two sequences match,
and a space where ever they do not.

Dots may be inserted in either sequence to represent gaps.

For example, the two sequences

    abcdefghajklm
    abbdhijk

could be aligned like this

    abcdefghajklm
    || |   | || 
    abbd...hijk

As shown, there are 6 matches, 2 mismatches, and one gap of length 3.

C<Align::NW> retuns an alignment as a hash

    $align = { a => $a,
	       s => $s,
	       b => $b };

I<$a> and I<$b> are the two sequences. 
I<$s> is the line of symbols.

=head2 The Payoff Matrix

The alignment is scored according to a payoff matrix

    $payoff = { match      => $match,
		mismatch   => $mismatch,
		gap_open   => $gap_open,
		gap_extend => $gap_extend };

The entries in the matrix are the number of points added to the score 

=over

=item *

for each match

=item *

for each mismatch

=item *

when a gap is opened in either sequence

=item *

for each position that a gap is extended (including the first)

=back

For correct operation, match must be positive, 
and the other entries must be negative.

=head2 Example

Given the payoff matrix

   $payoff = { match      =>  4,
	       mismatch   => -3,
	       gap_open   => -2,
	       gap_extend => -1 };

The sequences

    abcdefghajklm
    abbdhijk

are aligned and scored like this

                a b c d e f g h a j k l m
                | |   |       |   | | 
                a b b d . . . h i j k

    match       4 4   4       4   4 4  
    mismatch       -3          -3
    gap_open           -2
    gap_extend         -1-1-1

for a total score of 24-6-2-3 = 15. 
The algorithm guarantees that no other alignment of these two sequences
has a higher score under this payoff matrix.

=head1 METHODS

=over 4

=item I<$nw> = C<new> C<Align::NW> I<$a>, I<$b>, I<$payoff>, I<%options>

Creates and returns a new C<Align::NW> object.
I<$a> and I<$b> are the sequences to be aligned.
I<$payoff> is the payoff matrix, described above.
Additional options maybe passed in the I<%options> hash;
see L</OPTIONS> for details.

=item I<$nw>->C<score>

Fills in the score matrix for I<$nw>.
This is the O(a*b*(a+b)) operation.

=item I<$nw>->C<align>

Backtracks through the score matrix and generates an alignment for 
the two sequences.
C<score> must be called before C<align>.

=item I<$score> = I<$nw>->C<get_score>

Returns the score of the alignment.
C<score> must be called before C<get_score>.

=item I<$align> = I<$nw>->C<get_align>

Returns the alignment of the two sequences, as described above in 
L</Alignments>.
C<align> must be called before C<get_align>.

=item I<$nw>->C<print_align>

Pretty prints the alignment to STDOUT.
C<align> must be called before C<print_align>.

=item I<$nw>->C<dump_score>

Dumps the score matrix to STDOUT.
This is useful mainly for debugging.
The matrix is I<not> pretty printed;
line wrapping makes large matrices difficult to read.

=back

=head1 OPTIONS

Options may be passed to C<new> in the C<%options> hash.
The following options are defined.

=over 4

=item B<-v>

Verbose output.
Prints some dots to STDERR.
Useful for monitoring the progress of large alignments.

=back

=head1 SEE ALSO

=over 4

=item *

Needleman, S.B. and Wunsch, C.D. 1970. "A general method applicable to
the search for similarities in the amino acid sequences of two
proteins" I<Journal of Molecular Biology>. 48: 443-453.

=item *

Smith, T.F. and Waterman, M.S. 1981. "Identification of common
molecular subsequences" I<Journal of Molecular Biology>. 147: 195-197

=back

There are usually some some tutorials on Needleman-Wunsch and
Smith-Waterman alignment floating around on the web. I used to provide
links to some, but they kept going 404. If you Google around a bit you
can probably find a current one.

=head1 ACKNOWLEDGMENTS

=over 4

=item *

Andreas Doms <ad11@inf.tu-dresden.de>

=back

=head1 AUTHOR

Steven McDougall <swmcd@world.std.com>

=head1 COPYRIGHT

Copyright 1999-2003 by Steven McDougall. This module is free
software; you can redistribute it and/or modify it under the same
terms as Perl itself.

=cut
