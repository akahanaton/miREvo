#!/bin/bash
#--------------------------------------------------
# shopt -s -o nounset
#-------------------------------------------------- 

function USAGE {
	echo ""
	echo "Usage: miREvo filter -o prefix -i reads.fas -d database.fasta -H known.miRNA -M mature.fa [options]"
	echo "    Options for filter reads"
	echo "        -i  <str>   sequence reads file in FASTA format, uniqe merged, required"
	echo "        -d  <str>   The prefix of the bowtie index for Repeat Database, constructed from a Fasta file contain known tRNAs, sRNA, etc., recommmanded"
	echo "                    For instance, if the reference is 'database.fasta', then the prefix is 'database' and building-command is:"
	echo "                    'bowtie-build -f database.fasta database'"
	echo "        -H  <str>   Botwtie index for known miRNAs' hairpin sequences. These should be the known haripin sequences for the species being analyzed."
	echo "                    the miRNAs' ID must be in miRBase format, for example '>dme-miR-1'"
	echo "        -M  <str>   fasta file with miRNAs mature. These should be the known mature sequences for the species being analyzed."
	echo "                    the miRNAs' ID must be in miRBase format, for example '>dme-miR-1-5p'"
	echo "        -o  <str>   abbreviation for project name, 3 letter code for the sequencing library or the species of interest,required"
	echo ""
	echo "    Options for Bowtie:"
	echo "        -v  <int>   maximum number of mismatches allowed on a read when mapping to the Repeat Database, <=2. default=2bp"
	echo "        -p  <int>   number of processors to use, default=1"
}

if [ $# -eq 0 ]; then
	USAGE;
	exit 192;
fi

declare -rx SCRIPT=${0##*/}
declare -r OPTSTRING="i:d:H:M:o:v:p:n:k:h"
declare SWITCH
declare DATABASE 
declare hairpin=''
declare MATURE=''
declare INPUT
declare prj
declare PROJNAME 
declare MAXREP=5
declare -i MIS=2
declare -i CPU=1
declare -i Ns=5
declare -i SEED=8

program_dir="$MIREVO/script"

while getopts "$OPTSTRING" SWITCH ; do
	case $SWITCH in
		h) USAGE;
		   exit 0
		;;
		i) INPUT="$OPTARG"
		;;
		d) DATABASE="$OPTARG"
		;;
		H) hairpin="$OPTARG"
		;;
		M) MATURE="$OPTARG"
		;;
		o) prj="$OPTARG"
		;;
		p) CPU="$OPTARG"
		;;
		v) MIS="$OPTARG"
		;;
		k) MAXREP="$OPTARG"
		;;
		\?) exit 192
		;;
		*) printf "miREvo filter: $LINENO: %s\n" "script error: unhandled argument"
		exit 192
		;;
	esac
done

if [ ! -e $prj ]; then
	mkdir $prj
fi

PROJNAME=`basename $prj`

DB_BWT=$prj/filter.db.bwt
DB_NOMAP=$prj/filter.db.nomap.fas
DB_MAP=$prj/filter.db.map.fas

KN_BWT=$prj/filter.kn.bwt
KN_NOMAP=$prj/filter.kn.nomap.fas
KN_MAP=$prj/filter.kn.map.fas
CMDLOG=$prj/filter.cmd

declare db_name
declare hairpin_seq

if [[  -e $hairpin.fa  ]]; then
	hairpin_seq=$hairpin.fa
elif [[ -e $hairpin.fas  ]]; then
	hairpin_seq=$hairpin.fas
elif [[ -e $hairpin.fasta  ]]; then
	hairpin_seq=$hairpin.fasta
else
	echo "Can't locate the fasta sequence for $hairpin";
	echo "Please provide an avaliabe $hairpin sequence, such as";
	echo "$hairpin.fa, $hairpin.fas, $hairpin.fasta, etc."
	echo "and build the Bowtie index using a command like:"
	echo "bowtie-build -f $hairpin.fa $hairpin"
	echo "exit now."
	exit 192;
fi

echo ""
if [ -e $DATABASE.1.ebwt ] ; then
	db_name=`basename $DATABASE`
	echo "Filtering reads mapped to $db_name out"
	echo "bowtie -f -v $MIS -a --best --strata --suppress 5,6,7 $DATABASE $INPUT --al $DB_MAP --un $DB_NOMAP -p $CPU > $DB_BWT" > $CMDLOG
	bowtie -f -v $MIS -a --best --strata --suppress 5,6,7 $DATABASE $INPUT --al $DB_MAP --un $DB_NOMAP -p $CPU > $DB_BWT

else
	echo "Warnning: miREvo cann't locate the Bowtie index files for repeat database $DATABASE"
	echo "skip filtering reads in the $DATABASE"
	echo "cp $INPUT $DB_NOMAP" > $CMDLOG
	cp $INPUT $DB_NOMAP
fi



echo ""
echo "Filtering reads mapped to known miRNAs out"

if [ -e $hairpin.1.ebwt ] ; then
	echo "bowtie -f -v 1 -a --best --strata $hairpin $DB_NOMAP --al $KN_MAP --un $KN_NOMAP -p $CPU > $KN_BWT" >> $CMDLOG
	bowtie -f -v 1 -a --best --strata $hairpin $DB_NOMAP --al $KN_MAP --un $KN_NOMAP -p $CPU > $KN_BWT

	echo "perl $program_dir/analysis_filter.pl  $hairpin_seq $MATURE $KN_BWT $DB_BWT $INPUT > $prj/filter.statistic" >> $CMDLOG
	perl $program_dir/analysis_filter.pl  $hairpin_seq $MATURE $KN_BWT $DB_BWT $INPUT > $prj/filter.statistic

	echo "perl $program_dir/mirna_tag_aln_fasta.pl $hairpin_seq $KN_BWT > $prj/filter.kn.map" >> $CMDLOG
	perl $program_dir/mirna_tag_aln_fasta.pl $hairpin_seq $KN_BWT > $prj/filter.kn.map
else
	echo "Warnning: miREvo cann't locate the Bowtie index files for known miRNAs $hairpin_seq"
	echo "skip filtering reads in the $hairpin_seq"
	echo "ln -s `basename $DB_NOMAP` $KN_NOMAP" >> $CMDLOG
	if [ -e $KN_NOMAP ] ;  then
		rm $KN_NOMAP
		ln -s `basename $DB_NOMAP` $KN_NOMAP
	else
		ln -s `basename $DB_NOMAP` $KN_NOMAP
	fi
fi

if [[ ! -e $prj/filter.db.bwt && ! -e $prj/filter.kn.bwt ]] ; then
	echo "Warrning: "
	echo "Failed to locate Botie Index files, nothing done for filter, exit."
	echo "If you wish to continue, run predict command nextly."
	exit 192;
fi

echo ""
echo "Filter successfully done."
