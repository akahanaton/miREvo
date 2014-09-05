use strict;
use warnings; 

open(IN,$ARGV[0]);
open(OUT,">$ARGV[1]");

while(my $line = <IN>)
{
	if($line =~/^>/){
		chomp $line;
		my ($id,undef) = split(/\s+/, $line);
		print OUT $id,"\n";
	}else{
		$line =~ tr/uU/tT/;	
		print OUT $line;
	}
}
close IN;
close OUT;
