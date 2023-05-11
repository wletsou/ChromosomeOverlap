#! /bin/bash
set -e

# bsub -P SJLIFE -J ukbb_topmed_merge.chr11.129361171-129561171 -oo ukbb_topmed_merge.chr11.129361171-129561171.out -eo ukbb_topmed_merge.chr11.129361171-129561171.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_topmed_merge.sh ukbb_snp_list.txt /research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/chr11/chr11.qced.anno.info 11 129361171,129561171 /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/chr11"

SNPS=$1 # optional comma-separated list of snps or file with one snp per line
INFO=$2
CHR=$3 # single chromosome number
BP_RANGE=$4 # comma-separated list from_bp,to_bp of region on chromosome CHR
VCF=$5 # location of vcf files to be merged
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
  range=($(awk 'NR>1{b=gensub("chr'$CHR'[:]([0-9]*)[:].*","\\1","g",$7); if (b+0>='${BP_RANGE[0]}'+0 && b+0<='${BP_RANGE[1]}'+0) {print $3} }' $INFO | sort -gk1 | awk 'NR==1{print $0} {row=$0} END{print row}')) # prints first and last line of sorted list of within-limits hg38 positions (3rd field), based on h19 value in 7th field
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

if [ -z $VCF ]
then
  VCF=" /research_jude/rgs01_jude/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/chr${CHR}"
fi

test -f vcf_list.txt && rm vcf_list.txt
touch vcf_list.txt
for subdir in $VCF/b*
do
  echo bcftools view \-m2 \-M2 \-v snps \-r $region_str ${subdir}/chr${CHR}.dose.vcf.gz \-O z \-o ${subdir##*/}.vcf.gz
  bcftools view -m2 -M2 -v snps -r $region_str ${subdir}/chr${CHR}.dose.vcf.gz -O z -o ${subdir##*/}.vcf.gz
  echo tabix \-p vcf ${subdir##*/}.vcf.gz
  tabix -p vcf ${subdir##*/}.vcf.gz
  test -f ${subdir##*/}.vcf.gz && echo echo ${subdir##*/}.vcf.gz \>\> vcf_list.txt && printf "\n"
  test -f ${subdir##*/}.vcf.gz && echo ${subdir##*/}.vcf.gz >> vcf_list.txt
done

# echo bcftools merge \-l vcf_list.txt \-O z \-m snps -o $VCF/ukbb.topmed.hg19_chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${str}.hg38.vcf.gz
# bcftools merge -l vcf_list.txt -O z -m snps -o $VCF/ukbb.topmed.hg19_chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${str}.hg38.vcf.gz
# echo tabix \-p vcf $VCF/ukbb.topmed.hg19_chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${str}.hg38.vcf.gz
# tabix -p vcf $VCF/ukbb.topmed.hg19_chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${str}.hg38.vcf.gz

echo bcftools merge \-l vcf_list.txt \-O z \-m snps -o ukbb.topmed.hg19_chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${str}.hg38.vcf.gz
bcftools merge -l vcf_list.txt -O z -m snps -o ukbb.topmed.hg19_chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${str}.hg38.vcf.gz
echo tabix \-p vcf ukbb.topmed.hg19_chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${str}.hg38.vcf.gz
tabix -p vcf ukbb.topmed.hg19_chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${str}.hg38.vcf.gz

for file in b*.vcf*
do
  test -f $file && echo rm $file
  test -f $file && rm $file
done

test -f vcf_list.txt && echo rm vcf_list.txt
test -f vcf_list.txt && rm vcf_list.txt
