#! /bin/bash
set -e

# bsub -P SJLIFE -J 1000gp_haplotype_extract -oo 1000gp_haplotype_extract.out -eo 1000gp_haplotype_extract.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/1000gp_haplotype_extract.sh eur.txt 11 68881416,69455872 /scratch_space/wletsou/sjlife/GWAS/1000gp_chr11.0 /home/wletsou/scripts"

# bsub -P SJLIFE -J 1000gp_haplotype_extract -oo 1000gp_haplotype_extract.out -eo 1000gp_haplotype_extract.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/1000gp_haplotype_extract.sh eur.txt 11 69000000,70000000 /scratch_space/wletsou/sjlife/GWAS/1000gp_chr11.1 /home/wletsou/scripts"

# bsub -P SJLIFE -J 1000gp_haplotype_extract -oo 1000gp_haplotype_extract.out -eo 1000gp_haplotype_extract.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/1000gp_haplotype_extract.sh eur.txt 11 68850000,69500000 /scratch_space/wletsou/sjlife/GWAS/1000gp_chr11.2 /home/wletsou/scripts"

# bsub -P SJLIFE -J 1000gp_haplotype_extract -oo 1000gp_haplotype_extract.out -eo 1000gp_haplotype_extract.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/1000gp_haplotype_extract.sh eur.txt 11 68850000,69500000 /scratch_space/wletsou/sjlife/GWAS/1000gp_chr11.3 /home/wletsou/scripts"

# bsub -P SJLIFE -J 1000gp_haplotype_extract -oo 1000gp_haplotype_extract.out -eo 1000gp_haplotype_extract.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/1000gp_haplotype_extract.sh eur.txt 11 129361171,129561171 /scratch_space/wletsou/sjlife/GWAS/1000gp_chr11.4 /home/wletsou/scripts"

POPULATION=$1 # txt file of population ids
CHR=$2 # single chromosome number
BP_RANGE=$3 # comma-separated list from_bp,to_bp of regions on chromosome CHR
DIRECTORY=$4 # current working directory, specify replicate folder if REPLICATE_NO > 0
HOME_DIR=$5 # location of scripts to be run

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY
echo Directory is $DIRECTORY
printf "\n"

POPULATION=($(echo $POPULATION | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))
n_pop=${#POPULATION[@]}

BP_RANGE=($(echo $BP_RANGE | perl -pne 's/([0-9]+)[,]*/$1 /g'))

queue="large_mem" # specify queue for large-memory tasks

VCF=$(echo "/research/rgs01/project_space/yasuigrp/EpiGenetics/common/1000genomes/GRCh37") # location of 100gp vcf files

echo Create ${POPULATION[0]%.*}.indiv file:
echo awk \'BEGIN{OFS=\"\\t\"} {print \$1,\$1}\' ${POPULATION[0]} \> ${POPULATION[0]%.*}.indiv
awk 'BEGIN{OFS="\t"} {print $1,$1}' ${POPULATION[0]} > ${POPULATION[0]%.*}.indiv # doubled file of subject ids
printf "\n"

echo Extract haplotypes:
echo bsub \-P SJLIFE \-J 1000gp_hybrid_haplotype.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo ${DIRECTORY}/1000gp_hybrid_haplotype.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY}/1000gp_hybrid_haplotype.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-R \"rusage[mem=10000]\" \-q $queue \-K \"sh ${HOME_DIR}/ukbb_hybrid_haplotype.sh ${POPULATION[0]%.*}.indiv ${CHR} ${3} ${DIRECTORY} ${HOME_DIR} ${VCF}\"
bsub -P SJLIFE -J 1000gp_hybrid_haplotype.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo ${DIRECTORY}/1000gp_hybrid_haplotype.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY}/1000gp_hybrid_haplotype.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -R "rusage[mem=10000]" -q $queue -K "sh ${HOME_DIR}/ukbb_hybrid_haplotype.sh ${POPULATION[0]%.*}.indiv ${CHR} ${3} ${DIRECTORY} ${HOME_DIR} ${VCF}" # use original input BP_RANGE=$3 to resample haplotypes in larger population
printf "\n"

echo Wait until haplotypes extracted:
echo bsub \-P SJLIFE \-J sleep_0.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.overlap \-oo ${DIRECTORY}/sleep_0.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.overlap.out \-eo ${DIRECTORY}/sleep_0.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.overlap.err \-w \"done\(1000gp_hybrid_haplotype.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}\)\" \-R \"rusage[mem=32]\" \-K \"sleep 10\"
bsub -P SJLIFE -J sleep_0.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.overlap -oo ${DIRECTORY}/sleep_0.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.overlap.out -eo ${DIRECTORY}/sleep_0.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.overlap.err -w "done(1000gp_hybrid_haplotype.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]})" -R "rusage[mem=32]" -K "sleep 10"
printf "\n"

echo Keep only rs SNPs:
echo awk \'NR==1{printf \"%s\",\$1\; for \(i=2\;i\<NF\;i++\) {if \(\$i \~ /rs[0-9]*_[A-Z]\$/\) {included[i]=1\; printf \"\\t%s\",\$i} } if \(\$NF \~ /rs[0-9]*_[A-Z]\$/\) {included[NF]=1\; printf \"\\t%s\\n\",\$NF} else {printf \"\\n\"} } NR\>1{printf \"%s\",\$1\; for \(i=2\;i\<NF\;i++\) {if \(included[i]==1\) {printf \"\\t%s\",\$i} } if \(included[NF]==1\) {printf \"\\t%s\\n\",\$NF} else {printf \"\\n\"} }\' haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt \> haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.temp.txt
awk 'NR==1{printf "%s",$1; for (i=2;i<NF;i++) {if ($i ~ /rs[0-9]*_[A-Z]$/) {included[i]=1; printf "\t%s",$i} } if ($NF ~ /rs[0-9]*_[A-Z]$/) {included[NF]=1; printf "\t%s\n",$NF} else {printf "\n"} } NR>1{printf "%s",$1; for (i=2;i<NF;i++) {if (included[i]==1) {printf "\t%s",$i} } if (included[NF]==1) {printf "\t%s\n",$NF} else {printf "\n"} }' haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt > haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.temp.txt
printf "\n"

echo cat haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.temp.txt \> haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
cat haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.temp.txt > haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
echo test \-f haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.temp.txt \&\& rm haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.temp.txt
test -f haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.temp.txt && rm haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.temp.txt
printf "\n"

echo Included SNPs:
head -1 haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
awk 'NR==1{printf "%s SNP%s\n",NF-1,(NF-1==1?"":"s")}' haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
printf "\n"

# echo Prune LD SNPs from haplotypes files: # use ${REF_POPULATION[0]} as to determine which SNPs are in LD
# echo sh ${HOME_DIR}/1000G_LD_prune.sh haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt ${POPULATION[0]} $CHR ${BP_RANGE[0]},${BP_RANGE[1]} $DIRECTORY $HOME_DIR
# sh ${HOME_DIR}/1000G_LD_prune.sh haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt ${POPULATION[0]} $CHR ${BP_RANGE[0]},${BP_RANGE[1]} $DIRECTORY $HOME_DIR

echo Included SNPs:
head -1 haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
awk 'NR==1{printf "%s SNP%s\n",NF-1,(NF-1==1?"":"s")}' haplotype_estimates.${POPULATION[0]%.*}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
printf "\n"

prune_array=($(for file in *prune.out
do
  test -f $file && echo $file
done)) # look for file of pruned snps
if ((${#prune_array[@]}>0))
then
  declare -p prune_array
  echo Make pruned vcf file:
  for vcf_file in *.vcf
  do
    echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{rsid[\$1]\; next} \(!\(\$3 in rsid\) \&\& \(\$3 \~ /rs[0-9]*\$/\)\){print \$0}\' ${prune_array[0]} ${vcf_file} \> ${vcf_file%*.vcf}.prune.vcf
    awk 'BEGIN{OFS="\t"} NR==FNR{rsid[$1]; next} (!($3 in rsid) && ($3 ~ /rs[0-9]*$/)){print $0}' ${prune_array[0]} ${vcf_file} > ${vcf_file%*.vcf}.prune.vcf # print only lines not containing pruned SNPs
  done
  printf "\n"
fi
