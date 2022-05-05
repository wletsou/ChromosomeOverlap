#! /bin/bash
set -e

# Convert a vcf file into UCSC bed format at select SNPs

# sh /home/wletsou/scripts/vcf2bed.sh eur.indiv.chr11.69000000-70000000.haplotypes.vcf "rs79373485,rs1122316,rs76809977,rs559664,rs657315,rs79442425,rs111929748,rs2046494,rs4980661,rs637185"

# sh /home/wletsou/scripts/vcf2bed.sh ukbb_bca_cases.indiv.chr11.68700000-69700000.haplotypes.vcf "rs61881030,rs149949063,rs1123665,rs3829241,rs74674490,rs72930631,rs12577213,rs73520512,rs3892895,rs74897859,rs17308715,rs67039008,rs17149775,rs117236867,rs11228530,rs17309046,rs72932540,rs78033785,rs74342245,rs4495899,rs117694794,rs10792029,rs10896445,rs78208050,rs61881109,rs75590802,rs7937094,rs10896449,rs111762835,rs4930672,rs116951063,rs118004919,rs35637432,rs34737133,rs11825015,rs78656034,rs74551015,rs117821797,rs12275849,rs56103266,rs117131950,rs7103126,rs116926312,rs11539762,rs142581206,rs12274095,rs73512137,rs11601693,rs10896461,rs78919175,rs7124547,rs11603814,rs72932105,rs12789955,rs115223540,rs11600497,rs72932198,rs12802601,rs61882193,rs12790064,rs11263638,rs76855139,rs74471298,rs10160464,rs4980785,rs7105934,rs6606643,rs79373485,rs1122316,rs76809977,rs559664,rs657315,rs79442425,rs111929748,rs637185,rs4980661,rs2046494"

VCF=$1 # vcf file
SNPS=$2 # comma-separated list of rsids
STR=$3 # string to prepend field 1 with, e.g., "chr" if vcf does not supply "chr" to chromosome numbers
OUTPUT=$4 # optional output .bed
DIRECTORY=$5 # current working directory, specify replicate folder if REPLICATE_NO > 0
HOME_DIR=$6 # location of scripts to be run

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY

if [ -z $SNPS ]
then
  SNPS=($(cut -f3 $VCF)) # array of rsids
  SNPS=$(echo ${SNPS[*]} | sed 's/ /,/g') # print array on one line and replace spaces with commas
fi

if [ -z $OUTPUT ]
then
  OUTPUT=${VCF%.*}.bed
fi

awk 'BEGIN{OFS="\t"; n=split("'$SNPS'",array,","); for (i=1;i<=n;i++) {snps[array[i]]}; printf "track name=\"Haplotype SNPs\" description=\"Haplotype SNPs\" visibility=2 itemRgb=\"On\"\n" } ($3 in snps){print "'$STR'"$1,$2-1,$2,$3,0,".",($2-1),$2,"255,0,0"}' $VCF > $OUTPUT
