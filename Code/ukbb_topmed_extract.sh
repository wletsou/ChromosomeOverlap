#! /bin/bash
set -e

# bsub -P SJLIFE -J ukbb_topmed_extract.chr16.52499188-52699188 -oo ukbb_topmed_extract.chr16.52499188-52699188.out -eo ukbb_topmed_extract.chr16.52499188-52699188.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_topmed_extract.sh ukbb_snp_list.chr16.52499188-52699188.txt /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr16.51560000.53560000.qced.anno.info 16 52499188,52699188 /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/ukb.bca.topmedr2.chr16.51560000.53560000.qced.hg38.vcf.gz"

SNPS=$1 # optional comma-separated list of snps or file with one snp per line
INFO=$2 # for converting between hg19 and hg38
CHR=$3 # single chromosome number
BP_RANGE=$4 # comma-separated list from_bp,to_bp of region on chromosome CHR, use hg19
VCF_FILE=$5 # vcf file to be extracted from
DIRECTORY=$6 # PWD by default
HOME_DIR=$7 # location of program files

module load bcftools/1.10.2
module load tabix/0.2.6

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY

BP_RANGE=($(echo $BP_RANGE | perl -pne 's/([0-9]+)[,]*/$1 /g'))

if [ ! -z $SNPS ]
then
  if [ -f $SNPS ]
  then
    SNPS=$(cat $SNPS | awk '($1 ~ /rs[0-9]*/){printf "%s%s",(found>0?",":""),$1; found+=1} END{printf "\n"}') # snps supplied as a file, translated to comma-separated list, only keep rs snps
  fi
  region_str=$(awk 'BEGIN{n=split("'${SNPS}'",array,","); for (i=1;i<=n;i++) {snps[array[i]]} } ($8 in snps){printf "%schr%s:%s",(found>0?",":""),$2,$3; found+=1} END{printf "\n"}' $INFO) # get hg38 chromosome & coordinates of each snp
fi
if [ -z $region_str ]
then
  # find hg38 limits corresponding to hg19 limits, chromosome field does begin with "chr"
  range=($(awk 'NR>1{b=gensub("chr'$CHR':([0-9]*):.*","\\1","g",$7); if (b>='${BP_RANGE[0]}' && b<='${BP_RANGE[1]}') {print $3} }' $INFO | sort -gk1 | awk 'NR==1{print $0} {row=$0} END{print row}')) # prints first and last line of sorted list of within-limits hg38 positions (3rd field), based on h19 value in 7th field
  # range=($(awk 'NR>1{b=gensub("'$CHR'[.]([0-9]*)[.].*","\\1","g",$7); if (b+0>='${BP_RANGE[0]}'+0 && b+0<='${BP_RANGE[1]}'+0) {print $3} }' $INFO | sort -gk1 | awk 'NR==1{print $0} {row=$0} END{print row}')) # prints first and last line of sorted list of within-limits hg38 positions (3rd field), based on h19 value in 7th field
  region_str=chr${CHR}:${range[0]}\-${range[1]}
fi
if [ -z $SNPS ]
then
  str=".imputed" # string indicating all imputed SNPs in region
fi
echo $str
echo $region_str
printf "\n"

if [ -z $VCF_FILE ]
then
  VCF_FILE="/research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/drive4william/ukb.bca.topmedr2.chr16.51560000.53560000.qced.hg38.vcf.gz"
fi

echo Compress and index vcf file:
echo bcftools view \-r $region_str \-m2 \-M2 \-v snps \--output-type z \--output-file ukbb.topmed.hg19_chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${str}.hg38.vcf.gz $VCF_FILE
bcftools view -r $region_str -m2 -M2 -v snps --output-type z --output-file ukbb.topmed.hg19_chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${str}.hg38.vcf.gz $VCF_FILE
printf "\n"

echo tabix \-p vcf ukbb.topmed.hg19_chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${str}.hg38.vcf.gz
tabix -p vcf ukbb.topmed.hg19_chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${str}.hg38.vcf.gz
printf "\n"
