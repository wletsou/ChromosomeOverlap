#! /bin/bash
set -e

VCF=$1 # vcf file
SNPS=$2 # comma-separated list of rsids,or haplotype rsid_allele=[0/1],...
STR=$3 # string to prepend field 1 with, e.g., "chr" if vcf does not supply "chr" to chromosome numbers
OUTPUT=$4 # optional output .bed
HEADER=$5 # 1 to include a header line in bed file
DIRECTORY=$6 # current working directory, specify replicate folder if REPLICATE_NO > 0
HOME_DIR=$7 # location of scripts to be run

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

if [ -z $HEADER ]
then
  HEADER=1
fi

echo $SNPS
# chops allele and captures its value, if it exists, in "b"; otherwise returns rsid (a) and 0 (b)
awk 'BEGIN{OFS="\t"; n=split("'$SNPS'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)_[A-Z]*=(.*)","\\1","g",array[i]); b=gensub("(.*)_[A-Z]*=(.*)","\\2","g",array[i]); snps[a]=b+0}; if ('$HEADER'==1) {printf "track name=\"Haplotype SNPs\" description=\"Haplotype SNPs\" visibility=2 itemRgb=\"On\"\n"} } ($3 in snps){print "'$STR'"$1,$2-1,$2,$3,0,".",($2-1),$2,(snps[$3]==0?"0,0,255":"255,0,0")}' $VCF > $OUTPUT
