#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Std;
use File::Copy;
use File::Path;
use File::Basename qw/basename dirname/;
use Cwd qw(abs_path);



my $usage=
"$0 file_command_line file_structure rounds_controls prj_dir

-a   Output progress to screen
";

my $file_command_line=shift or die $usage;
my $file_structure=shift or die $usage;
my $rounds=shift or die $usage;
my $prj_dir=shift or die $usage;

#options
my %options=();
getopts("a",\%options);


my $dir="$prj_dir/dir_perform_controls";
my $path = abs_path($0);
$path = dirname($path);



my $command_line=parse_file_command_line($file_command_line);

print $command_line,"\n";

if($options{a}){print STDERR "total number of rounds controls=$rounds\n";}

perform_controls();


sub perform_controls{

    mkdir $dir;
    
    my $round=1;

    while($round<=$rounds){

	if($options{a}){print STDERR "$round\r";}

		system("$path/permute_structure.pl $file_structure > $dir/precursors_permuted.str 2> /dev/null");
		
		my $ret=`$command_line 2> /dev/null`;

		print "permutation $round\n\n";

		print "$ret";
		
		$round++;
    }

	rmtree($dir);

    if($options{a}){print STDERR "controls performed\n\n";}
}


sub parse_file_command_line{

    my ($file) = @_;

    open (FILE, "<$file") or die "can not open $file\n";
    while (my $line=<FILE>){

	if($line=~/(\S+)/){
	    
	    chomp $line;

	    $line=~s/$file_structure/$dir\/precursors_permuted.str/;
	    
	    $line=~s/>.+//;

	    return $line;

	}
    }
    die "$file is empty\n";
}




