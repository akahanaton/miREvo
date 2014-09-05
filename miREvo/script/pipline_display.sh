#!/bin/bash

shopt -s -o nounset


function USAGE {
	echo ""
	echo "Usage: miREvo display -o prefix -i reads.fasta -S miRNA.fasta [options]"
	echo "    Options for miRNA query"
	echo "        -o  <str>   abbreviation for project name, required"
	echo "        -i  <str>   fasta file with sequencing reads, required"
	echo "        -M  <str>   fasta file with MATURE sequence of known miRNAs, optional"
	echo "                    the miRNAs' ID must be in miRBase format, for example '>dme-miR-1-5p'"
	echo "        -H  <str>   fasta file with PRECURSOR sequence of known miRNAs, optional"
	echo "                    the miRNAs' ID must be in miRBase format, for example '>dme-miR-1'"
	echo "        -t  <int>   temperature cut-off for RNAfold when calculating secondary structures of RNAs, default=22"
	echo ""
	echo "    Options for Bowtie:"
	echo "        -v  <int>   maximum number of mismatches allowed on a read, <=2. default=0bp"
	echo "        -p  <int>   number of processors to use, default=1"
	echo "        -n  <int>   max mismatches in seed (can be 0-3, default: -n 0) (mismatch mapping takes longer)"
	echo "        -k  <int>   a read is allowed to map up to this number of positions in the genome, default=5; for plant, 15 is recommaned"
}

if [ $# -eq 0 ]; then
	USAGE;
	exit 192;
fi

declare -rx SCRIPT=${0##*/}
declare -r OPTSTRING="i:o:M:H:t:p:v:n:k:z:h"
declare SWITCH
declare mature
declare hairpin 
declare READS=' '
declare prj
declare PROJNAME 
declare -i TEMP=22
declare -i MIS=0
declare -i CPU=1
declare -i SEED_MIS=0
declare -i MAXREP=5

program_dir="$MIREVO/script"

while getopts "$OPTSTRING" SWITCH ; do
	case $SWITCH in
		h)  USAGE;
			exit 192;
		;;
		i) READS="$OPTARG"
		;;
		o) prj="$OPTARG"
		;;
		M) mature="$OPTARG"
		;;
		H) hairpin="$OPTARG"
		;;
		t) TEMP="$OPTARG"
		;;
		p) CPU="$OPTARG"
		;;
		v) MIS="$OPTARG"
		;;
		n) SEED_MIS="$OPTARG"
		;;
		k) MAXREP="$OPTARG"
		;;
		#--------------------------------------------------
		# z) DROSHA="$OPTARG"
		# ;;
		#-------------------------------------------------- 
		\?) exit 192
		;;
		*) printf "miREvo predict: $LINENO: %s\n" "script error: unhandled argument"
		exit 192
		;;
	esac
done

PROJNAME=`basename $prj`

CMDLOG=$prj/display.cmd

if [ ! -e $prj ]; then
	mkdir $prj
fi

if [  -e $CMDLOG ]; then
	rm $CMDLOG
fi

echo "Display begin ...." > $CMDLOG

MIR=$prj/display.mirna.fas
perl $program_dir/combine_mirna_mature.pl $PROJNAME $mature $hairpin > $MIR
perl $program_dir/clean_seq.pl $MIR $prj/display.mirna.fas.cleaned

num_of_mirna=`grep '>' $MIR | wc -l`
if [ $num_of_mirna -eq 0 ] ; then
	echo "No miRNA in query file $MIR, exit now."
	exit 192;
fi


if [ ! -e $prj/image ] ; then
	echo "mkdir $prj/image" >> $CMDLOG
	mkdir $prj/image
fi

echo "cd $prj/image" >> $CMDLOG
cd $prj/image
echo "cat ../display.mirna.fas.cleaned | RNAfold -noPS -T $TEMP -d0 | grep -A2 '>' | grep -v '\-\-' | awk '{print }' > display.mirna.fas.ss" >> ../display.cmd
cat ../display.mirna.fas.cleaned | RNAfold -noPS -T $TEMP -d0 | grep -A2 '>' | grep -v '\-\-' | awk '{print }' > display.mirna.fas.ss
echo "cat display.mirna.fas.ss | RNAplot -o svg" >> ../display.cmd
cat display.mirna.fas.ss | RNAplot -o svg > image.id
cd - >/dev/null

if [[ $READS == ' ' ]] ; then
	echo "No resequencing reads are provided, create a fake mapping result"
	echo "touch $prj/display.mirna.bwt" >> $CMDLOG
	touch $prj/display.mirna.bwt
else
	bowtie-build -f $MIR.cleaned $prj/display.mirna > $prj/bowtie-build.log
	echo "bowtie -p $CPU -v $MIS -f -n 0 -e 80 -l 18 -a -m 15 --best --strata $prj/display.mirna $READS > $prj/display.mirna.bwt" >> $CMDLOG
	bowtie -p $CPU -v $MIS -f -n 0 -e 80 -l 18 -a -m 15 --best --strata $prj/display.mirna $READS > $prj/display.mirna.bwt
fi

echo "perl $program_dir/svgShowMapDensity.pl $PROJNAME $MIR.cleaned $prj/display.mirna.bwt pred" >> $CMDLOG
perl $program_dir/svgShowMapDensity.pl $PROJNAME $MIR.cleaned $prj/display.mirna.bwt pred

echo "perl $program_dir/mirna_tag_aln_fasta_all.pl $MIR.cleaned $prj/display.mirna.bwt > $prj/display.map" >> $CMDLOG
perl $program_dir/mirna_tag_aln_fasta_all.pl $MIR.cleaned $prj/display.mirna.bwt > $prj/display.map

echo "perl $program_dir/analysis_filter.pl  $hairpin $mature /$prj/display.mirna.bwt $prj/display.statistic" >> $CMDLOG
perl $program_dir/analysis_filter.pl  $hairpin $mature /$prj/display.mirna.bwt $prj/display.statistic 


echo ""
echo ""
echo "Display successfully done."
