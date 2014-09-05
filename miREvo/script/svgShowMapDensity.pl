#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Std;
use File::Basename qw/basename dirname/;
use Data::Dumper;
use Bio::SeqIO;

my %opts;
my $usage = qq{
	Usage: $0 <lib.name> <predict.mirna>  <bwt_result> <pred/homo>
};

getopts('', \%opts);
die($usage) if (@ARGV != 4);
my ($libName,$fas,$bwt, $runTag) = @ARGV;


my(@svgFiles,@mirnaTags,%mirnasCov,%mapReadsNum, %mirnasLen);

my $seqin=Bio::SeqIO->new(-file=>$fas,-format=>"fasta");
while(my $seq=$seqin->next_seq()){
		my $seqID = $seq->display_id();

		my $inSVG = $libName."/image/".$seqID."_ss.svg";
		if ( $runTag eq 'homo' )
		{
			$inSVG =~ s/image/homo_image/
		}
		push(@svgFiles,$inSVG);
		my $mirnaTag;
		if ($seqID =~ m/([a-z]+_\d+)/){
			$mirnaTag = $1;
			push(@mirnaTags,$1);
		}else{
			print "unmatch\n";
		}
		my $seqLen = length($seq->seq());
		print $mirnaTag,"\n";
		$mirnasLen{$mirnaTag} =  $seqLen;
		@{$mirnasCov{$mirnaTag}} = (0) x $seqLen;
}


open(BWT, $bwt) || die("** Fail to open '$bwt'.\n");
while(my $line = <BWT>)
{
	my @tmpArr = split(/\s+/,$line);
	$tmpArr[2] =~ m/([a-z]+_\d+)/; #dme_1_cagcgtcagtcaagcggaagg_450.5
	my $mirnaTag = $1;
	print "$mirnaTag\n";
	my $readLen = length $tmpArr[4];
	next if $tmpArr[3] > $mirnasLen{$mirnaTag} - $readLen;
	if(defined($mirnasCov{$mirnaTag})){
		for(my $i =0; $i<$readLen; ++$i)
		{
			@{$mirnasCov{$mirnaTag}}[$tmpArr[3]+$i] += 1;
		}
	}
	if(defined($mapReadsNum{$mirnaTag})){
		$mapReadsNum{$mirnaTag} += 1;
	}else{
		$mapReadsNum{$mirnaTag} = 1;
	}
}
close BWT;

#--------------------------------------------------
# print Dumper($mirnasCov{'dmm_191'});
# print Dumper($mapReadsNum{'dmm_191'});
#-------------------------------------------------- 

# normolize each pos by total mapped reads num, one (mirna) by one
foreach my $mirna (keys %mirnasCov)
{
	if ( defined($mapReadsNum{$mirna})){
		for( 0..$#{$mirnasCov{$mirna}}){
			@{$mirnasCov{$mirna}}[$_] = 100 * @{$mirnasCov{$mirna}}[$_] / $mapReadsNum{$mirna};
		}
	}
}


for (my $i = 0; $i <= $#svgFiles; ++$i) 
{
	my $inSVG = $svgFiles[$i];
	open(INSVG, $inSVG) || die("** Fail to open '$inSVG'.\n");
	$_ = $inSVG;
	s/_ss/_ss_new/;
	my $outSVG= $_;
	open OUTSVG, "> $outSVG";

	my @baseColors;
	my @basePos = ();
	$/="<g";
	my ($polylineInfo,$pairsInfo,$seqInfo,$headInfo);

	while(my $line = <INSVG>)
	{
		if($line =~ /id="outline" points=/) {
			my @tmpArr = split(/\)\">\n/,$line);
			$headInfo = $headInfo.$tmpArr[0].")\">";
			$polylineInfo = $tmpArr[1];
		}elsif($line =~ /id="pairs"/) {
			$pairsInfo = $line;
		}elsif($line =~ /id="seq"/) {
			$seqInfo = $line;
			my @tmpArr = split(/\n/,$line);
			foreach my $tag (@tmpArr){
				if($tag =~m/<text\s+x=\"(.*)\"\s+y=\"(.*)\">([AGCTUagctu])<\/text>/){
					my $tmpBase = {"x"=>$1, "y"=>$2,"base"=>$3};
					push(@basePos, $tmpBase);
				}
			}
		}else{
			$headInfo = $line;
		}
	}
	close INSVG;

	my @headInfoArr = split( /\n/ , $headInfo );
	my $recordNum = $#headInfoArr-2;
	for (0..$recordNum){
		print OUTSVG $headInfoArr[$_],"\n";
	}

	# 0,255,0 green
	# 255 ,255,0 yellow
	# 255,0,0 red

	my $step = (255+255)/50;
	my ($red, $green, $blue) = (0,255,0);

	#print legend
	printf OUTSVG "<g id=\"legend_posen\" transform=\"translate(40 350)\" font-family=\"Arial,Helvetica\">\n";
	printf OUTSVG "<text x=\"0\" y=\"30\" font-size=\"12px\" fill=\"dimgray\">0</text>\n";
	my $xPos = 0;
	for (0..49){
		if($red < 255){
			$red += $step;
			$red = 255 if($red>255);
		}else{
			$red = 255;
			$green -= $step;
			$green = 0 if($green < 0);
		}
		$xPos += 2;
		printf OUTSVG "<rect style=\"stroke:rgb(%d,%d,0); stroke-width:0; fill:rgb(%d,%d,0)\" height=\"10\" x=\"%d\" y=\"30\" width=\"2\" />\n",
			$red,$green,$red,$green,$xPos;
	}
	printf OUTSVG "<text x=\"84\" y=\"30\" font-size=\"12px\" fill=\"dimgray\">100</text>\n";
	printf OUTSVG "</g>\n\n";

	print OUTSVG $headInfoArr[-1],"\n";
	#print circles
	my $totalBaseNum = $#basePos;
	for (0..$totalBaseNum){ 
		my $colorTag = @{$mirnasCov{$mirnaTags[$i]}}[$_]/2 * $step;
		if ($colorTag > 255 ){
			$red = 255;
			$green = 510 - $colorTag;
		}else{
			$red = $colorTag;
			$green = 255;
		}
		printf OUTSVG "<circle cx=\"%f\" cy=\"%f\" r=\"8\" style=\"fill:rgb(%d,%d,0)\" stroke=\"black\" stroke-width=\"0\"/>\n", 
					$basePos[$_]->{"x"},$basePos[$_]->{"y"},$red, $green;
	}

	#print polyline
	print OUTSVG $polylineInfo,$pairsInfo,$seqInfo;
	close OUTSVG;
	# convert svg to png
	$_ = $outSVG;
	s/svg/png/;
	my $pngfile= $_;
	system("convert $outSVG $pngfile");
} # end FAS
