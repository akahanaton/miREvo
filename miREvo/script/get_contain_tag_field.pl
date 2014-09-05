#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use Data::Dumper;

my %opts;
my $usage = qq{
	Usage: $0 <tab.sep.file> <id_file> <field> <Y/N>
	field count begin with 1, 
	Y or N
};

getopts('', \%opts);
die($usage) if (@ARGV != 4);
my ($tab_file, $id_file,$field, $YN) = @ARGV;

my (@tags, @lib1_info);


my %tagInfo;
open(F1, $id_file) || die("** Fail to open '$id_file'.\n");
while( my $line = <F1>)
{
	chomp $line;
	$tagInfo{$line} = 1;
}
close F1;

open(COMP, $tab_file) || die("** Fail to open $tab_file.\n");
while(my $line = <COMP>)
{
	chomp $line;
	my @arr = split(/\s+/, $line);
	if($YN eq 'Y'){
		if (defined($tagInfo{$arr[$field - 1]})){
			print $line,"\n";
		}
	}else{
		if (!defined($tagInfo{$arr[$field - 1]})){
			print $line,"\n";
		}
	}
}
close COMP;

