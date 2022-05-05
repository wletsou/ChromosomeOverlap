#! /bin/bash
set -e

# bsub -P SJLIFE -J dbgap_haplotype_extract4 -oo dbgap_haplotype_extract4.out -eo dbgap_haplotype_extract4.err -R "rusage[mem=10000]" -q large_mem "sh /home/wletsou/scripts/dbgap_haplotype_extract4.sh dbgap28544_cases.indiv,dbgap28544_controls.indiv "/research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/chr11/chr11.qced.anno.info" CCND1_TAD_snps.txt 11 68850000,69231641 /research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/vcfs/dbgap28544.topmedr2.cleaned.hg38.chr11.vcf.gz"

POPULATION=$1 # commas-separated list of population indiv files
INFO=$2 # SNP info file for converting hg38 to hg19
SNPS=$3 # comma-separated list of snps or file with one snp per line
CHR=$4 # single chromosome number, appended with "chr" if necessary to match VCF_FILE
BP_RANGE=$5 # comma-separated list from_bp,to_bp of region on chromosome CHR
VCF_FILE=$6 # location of vcf files
DIRECTORY=$7 # PWD by default
HOME_DIR=$8 # location of program files

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

BP_RANGE=($(echo $BP_RANGE | perl -pne 's/([0-9]+)[,]*/$1 /g'))

POPULATION=($(echo $POPULATION | perl -pne 's/[,]/ /g'))

if [ -z $seed ]
then
  seed=20200116
fi

test ! -z $INFO && (test -f $INFO || (>&2 echo "SNP info file not found."; exit 1) ) || (>&2 echo "No info file supplied."; exit 1)

if [ -z $VCF_FILE ]
then
  VCF=$(echo "/research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/vcfs/")
  cd $VCF
  array=($(ls | grep ${CHR[0]}'[_.].*vcf.gz$')) # format of vcf.gz files with chr CHR in their name; use ukbb_gds2vcf.R to conver gds files to vcf if vcf does not exist
  declare -p array
  file_str=${array[0]} # file name of the the .vcf.gz file
  cd $DIRECTORY
  VCF_FILE=${VCF}/${file_str}
else
  VCF=${VCF_FILE##*/} # location of VCF_FILE found by trimming before last forward slash /
fi
echo Phased haplotypes files is $VCF_FILE

if [ ! -f ${VCF_FILE}.csi ] && [ ! -f ${VCF_FILE}.tbi ] # create index if it does not exist
then
  echo Generate .tbi index file
  echo tabix \-p vcf $VCF_FILE
  bcftools tabix -p vcf $VCF_FILE
  printf "\n"
fi

# specify region to extract: selected snps or entire range
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
  region_str=chr${CHR}:${BP_RANGE[0]}\-${BP_RANGE[1]}
fi
echo $region_str
printf "\n"

samples_file=$(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').samples
echo Get unique individuals in supplied .indiv file in order:
echo awk \'{seen[NR]=\$1\;} END{n = length\(seen\)\; for \(i=1\;i\<=n\;i++\) {print seen[i]} }\' ${POPULATION[*]} \> $samples_file # unique samples in .indiv file for population i
awk '{seen[NR]=$1;} END{n = length(seen); for (i=1;i<=n;i++) {print seen[i]} }' ${POPULATION[*]} > $samples_file # unique samples in .indiv file for population i, printed in same order as .indiv file
printf "\n"

transpose_file=haplotype_estimates_transpose.$(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
echo Get doubled header \(each subject id twice\):
echo bcftools view  \-r $region_str \--samples-file $samples_file \--force-samples \--output-type v \--header-only $VCF_FILE \| awk \'{printf \"sjlife\\t\"\; for \(i=10\;i\<=NF\;i++\) {printf \"%s\\t%s%s\",\$i,\$i,\(i\<NF?\"\\t\":\"\\n\"\)} }\' \| tail \-1 \> $transpose_file
bcftools view -r $region_str --samples-file $samples_file --force-samples -m2 -M2 -v snps --output-type v --header-only $VCF_FILE | awk '{printf "sjlife\t"; for (i=10;i<=NF;i++) {printf "%s\t%s%s",$i,$i,(i<NF?"\t":"\n")} }' | tail -1 > $transpose_file
printf "\n"

new_vcf_file=$(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.haplotypes.hg38.vcf
echo Extract individuals from .vcf:
echo bcftools view \-r $region_str \--samples-file $samples_file \--force-samples \--output-type v \--no-header \--output-file $new_vcf_file $VCF_FILE
bcftools view -r $region_str --samples-file $samples_file --force-samples -m2 -M2 -v snps --output-type v --no-header --output-file $new_vcf_file $VCF_FILE
printf "\n"

echo Compress and index vcf file:
echo bcftools view \-r $region_str \--samples-file $samples_file \--force-samples \--output-type z \--output-file $(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.hg38.vcf.gz $VCF_FILE
bcftools view -r $region_str --samples-file $samples_file --force-samples -m2 -M2 -v snps --output-type z --output-file $(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.hg38.vcf.gz $VCF_FILE
printf "\n"

echo tabix \-p vcf $(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.hg38.vcf.gz
tabix -p vcf $(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.hg38.vcf.gz
printf "\n"


echo Translate to hg19: # column 6 of $INFO is chr:hg38_pos:ref:alt, column 8 is rsid, column 16 is hg19 position
echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{b=gensub\(\"chr.*:\([0-9]\*\):.*\",\"\\\\1\",\"g\",\$6\)\; rsid[b]=\$8\; next} \(\$2 in rsid\){\$3=rsid[\$2]\; print \$0}\' $INFO $new_vcf_file \> ${new_vcf_file}.tmp
awk 'BEGIN{OFS="\t"} NR==FNR{b=gensub("chr.*:([0-9]*):.*","\\1","g",$6); rsid[b]=$8; next} ($2 in rsid){$3=rsid[$2]; print $0}' $INFO $new_vcf_file > ${new_vcf_file}.tmp
printf "\n"

test -f ${new_vcf_file}.tmp && echo mv ${new_vcf_file}.tmp $new_vcf_file && printf "\n"
test -f ${new_vcf_file}.tmp && mv ${new_vcf_file}.tmp $new_vcf_file

echo Extract haplotypes as rsid x chromosome: # use most likely allele call
echo awk \'BEGIN{OFS=\"\\t\"} {printf \"%s_%s\\t\",\$3,\$5\; for \(i=10\;i\<=NF\;i++\) {p=gensub\(/\([0-9.]+\)[\|/]\([0-9.]+\).*/,\"\\\\1\",\"g\",\$i\)\; m=gensub\(/\([0-9.]+\)[\|/]\([0-9.]+\).*/,\"\\\\2\",\"g\",\$i\)\; printf \"%s\\t%s%s\",p,m,\(i\<NF?\"\\t\":\"\\n\"\)} }\' $new_vcf_file \>\> $transpose_file
awk 'BEGIN{OFS="\t"} {printf "%s_%s\t",$3,$5; for (i=10;i<=NF;i++) {p=gensub(/[0-9]+[|/][0-9]+:[0-9.]+:([0-9.]+),[0-9.]+:.*/,"\\1","g",$i); m=gensub(/[0-9]+[|/][0-9]+:[0-9.]+:[0-9.]+,([0-9.]+):.*/,"\\1","g",$i); printf "%s\t%s%s",p,m,(i<NF?"\t":"\n")} }' $new_vcf_file >> $transpose_file
printf "\n"

echo Check for alleles not present in some population:
echo awk \'NR==1{print \$0} NR\>1{for \(i=2\;i\<=NF\;i++\) {if \(\$i !\~ /[0-9]/\) {next} }\; print \$0 }\' $transpose_file \> ${transpose_file%.txt}.tmp
awk 'NR==1{print $0} NR>1{for (i=2;i<=NF;i++) {if ($i !~ /[0-9]/) {next} }; print $0 }' $transpose_file > ${transpose_file%.txt}.tmp
printf "\n"

test -f ${transpose_file%.txt}.tmp && echo mv ${transpose_file%.txt}.tmp $transpose_file && printf "\n"
test -f ${transpose_file%.txt}.tmp && mv ${transpose_file%.txt}.tmp $transpose_file

echo Remove unshared alleles from vcf file:
echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{b=gensub\(\"\(.*\)_.*\",\"\\\\1\",\"g\",\$1\)\; snp[b]\; next} \(\$3 in snp\){print \$0}\' $transpose_file $new_vcf_file \> ${new_vcf_file}.tmp
awk 'BEGIN{OFS="\t"} NR==FNR{b=gensub("(.*)_.*","\\1","g",$1); snp[b]; next} ($3 in snp){print $0}' $transpose_file $new_vcf_file > ${new_vcf_file}.tmp

test -f ${new_vcf_file}.tmp && echo mv ${new_vcf_file}.tmp $new_vcf_file && printf "\n"
test -f ${new_vcf_file}.tmp && mv ${new_vcf_file}.tmp $new_vcf_file

haplotypes_file=haplotype_estimates.$(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
echo Transpose haplotypes into chromosome x rsid:
echo awk \'BEGIN{OFS=\"\\t\"}\; {for\(j=1\;j\<=NF\;j++\) {a[NR,j]=\$j\; n_rows=NR\; n_cols=\(n_cols\<NF?NF:n_cols\)} } END{for \(j=1\;j\<=n_cols\;j++\) {for \(i=1\;i\<=n_rows\;i++\) {printf \"%s%s\",a[i,j],\(i==n_rows?\"\\n\":\"\\t\"\)} } }\' $transpose_file \> $haplotypes_file
awk 'BEGIN{OFS="\t"}; {for(j=1;j<=NF;j++) {a[NR,j]=$j; n_rows=NR; n_cols=(n_cols<NF?NF:n_cols)} } END{for (j=1;j<=n_cols;j++) {for (i=1;i<=n_rows;i++) {printf "%s%s",a[i,j],(i==n_rows?"\n":"\t")} } }' $transpose_file > $haplotypes_file
printf "\n"

for ((i=0;i<${#POPULATION[@]};i++))
do
  echo Extract ${POPULATION[i]%.indiv*} subjects as chromosome x rsid:
  echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{id[\$1]\; next} \(\$1 in id \|\| FNR==1\){print \$0}\' ${POPULATION[i]} $haplotypes_file \> haplotype_estimates.${POPULATION[i]%.indiv*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
  awk 'BEGIN{OFS="\t"} NR==FNR{id[$1]; next} ($1 in id || FNR==1){print $0}' ${POPULATION[i]} $haplotypes_file > haplotype_estimates.${POPULATION[i]%.indiv*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
  printf "\n"

  echo Transpose ${POPULATION[i]%.indiv*} haplotypes into rsid x chromosome:
  echo awk \'BEGIN{OFS=\"\\t\"}\; {for\(j=1\;j\<=NF\;j++\) {a[NR,j]=\$j\; n_rows=NR\; n_cols=\(n_cols\<NF?NF:n_cols\)} } END{for \(j=1\;j\<=n_cols\;j++\) {for \(i=1\;i\<=n_rows\;i++\) {printf \"%s%s\",a[i,j],\(i==n_rows?\"\\n\":\"\\t\"\)} } }\' haplotype_estimates.${POPULATION[i]%.indiv*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt \> haplotype_estimates_transpose.${POPULATION[i]%.indiv*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
  awk 'BEGIN{OFS="\t"}; {for(j=1;j<=NF;j++) {a[NR,j]=$j; n_rows=NR; n_cols=(n_cols<NF?NF:n_cols)} } END{for (j=1;j<=n_cols;j++) {for (i=1;i<=n_rows;i++) {printf "%s%s",a[i,j],(i==n_rows?"\n":"\t")} } }' haplotype_estimates.${POPULATION[i]%.indiv*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt > haplotype_estimates_transpose.${POPULATION[i]%.indiv*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
  printf "\n"
done
