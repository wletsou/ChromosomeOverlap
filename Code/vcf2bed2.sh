#! /bin/bash
set -e

# sh /home/wletsou/scripts/vcf2bed2.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.haplotypes.hg38.vcf "rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs74471298_A=0,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "" "h2.full-length.bed"

VCF=$1 # vcf file
SNPS=$2 # comma-separated list of rsids,or haplotype rsid_allele=[0/1],...
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
echo $SNPS
# chops allele and captures its value, if it exists, in "b"; otherwise returns rsid (a) and 0 (b)
awk 'BEGIN{OFS="\t"; n=split("'$SNPS'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)_[A-Z]*=(.*)","\\1","g",array[i]); b=gensub("(.*)_[A-Z]*=(.*)","\\2","g",array[i]); snps[a]=b+0}; printf "track name=\"Haplotype SNPs\" description=\"Haplotype SNPs\" visibility=2 itemRgb=\"On\"\n" } ($3 in snps){print "'$STR'"$1,$2-1,$2,$3,0,".",($2-1),$2,(snps[$3]==0?"255,0,0":"0,0,255")}' $VCF > $OUTPUT
