#! /bin/bash
set -e

# for generating a DRIVE vcf file from the region of interest

# cut -f13,14 dbgap28544.pheno.txt | awk 'BEGIN{OFS="\t"} (NR > 1 && $2 == 1){print $1,$1}' | sort -gk1 > dbgap28544_cases.indiv
# cut -f13,14 dbgap28544.pheno.txt | awk 'BEGIN{OFS="\t"} (NR > 1 && $2 == 0){print $1,$1}' | sort -gk1 > dbgap28544_controls.indiv

# bsub -P SJLIFE -J dbgap_haplotype_merge.chr11.69231642-69431642 -oo dbgap_haplotype_merge.chr11.69231642-69431642.out -eo dbgap_haplotype_merge.chr11.69231642-69431642.err -R "rusage[mem=1000]" "sh dbgap_haplotype_merge.sh dbgap28544_cases,dbgap28544_controls "dbgap28544.pheno.txt" "chr11.qced.anno.info" ukbb_snp_list.chr11.69231642-69431642.txt 11 69231642,69431642.txt dbgap28544.topmedr2.cleaned.hg38.chr11.vcf.gz"

POPULATION=$1 # commas-separated list of population names or indiv files
PHENO=$2 # phenotype file for determining case/control status of POPULATION
INFO=$3 # SNP info file for converting hg38 to hg19
SNPS=$4 # comma-separated list of snps or file with one snp per line
CHR=$5 # single chromosome number, appended with "chr" if necessary to match VCF_FILE
BP_RANGE=$6 # comma-separated list from_bp,to_bp of region on chromosome CHR
VCF_FILE=$7 # vcf file
DIRECTORY=$8 # PWD by default
HOME_DIR=$9 # location of program files

module load bcftools/1.10.2
module load tabix/0.2.6

if [ -z $HOME_DIR ];
then
  unset HOME_DIR
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

test ! -z $PHENO && (test -f $PHENO || (>&2 echo "Phenotype file not found."; exit 1) ) || (>&2 echo "No phenotype file supplied."; exit 1)
test ! -z $INFO && (test -f $INFO || (>&2 echo "SNP info file not found"; exit 1) ) || (>&2 echo "no info file supplied."; exit 1)

echo Extract cases and controls:
for ((i=0;i<${#POPULATION[@]};i++))
do
  if [ ! -f ${POPULATION[i]}.indiv ]
  then
    cut -f1,2 $PHENO | awk 'BEGIN{OFS="\t"} (NR > 1 && $2 == '$((${#POPULATION[@]}-i-1))'){print $1,$1}' | sort -gk1 > ${POPULATION[i]}.indiv # list of ids in columns 1 and 2, sorted in increasing order, based on affected status (1 = case POPULATION[0], 0 = control POPULATION[1])
  else
    echo ${POPULATION[i]}.indiv already exists.
  fi
done
printf "\n"

test -z $VCF_FILE && (>&2 echo "VCF file not supplied"; exit 1)
test -f $VCF_FILE || (>&2 echo "VCF file $VCF_FILE not found"; exit 1)

VCF=${VCF_FILE##*/} # location of VCF_FILE found by trimming before last forward slash /

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

samples_file=$(echo ${POPULATION[*]%.*} | sed 's/ /+/g').samples
echo Get unique individuals in supplied .indiv file in order:
echo awk \'{seen[NR]=\$1\;} END{n = length\(seen\)\; for \(i=1\;i\<=n\;i++\) {print seen[i]} }\' ${POPULATION[*]/%/.indiv} \> $samples_file # unique samples in .indiv file for population i
awk '{seen[NR]=$1;} END{n = length(seen); for (i=1;i<=n;i++) {print seen[i]} }' ${POPULATION[*]/%/.indiv} > $samples_file # unique samples in .indiv file for population i, printed in same order as .indiv file
printf "\n"

echo Compress and index vcf file:
echo bcftools view \-r $region_str \--samples-file $samples_file \--force-samples \--output-type z \--output-file $(echo ${POPULATION[*]%.*} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.hg38.vcf.gz $VCF_FILE
bcftools view -r $region_str --samples-file $samples_file --force-samples -m2 -M2 -v snps --output-type z --output-file $(echo ${POPULATION[*]%.*} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.hg38.vcf.gz $VCF_FILE
printf "\n"

echo tabix \-p vcf $(echo ${POPULATION[*]%.*} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.hg38.vcf.gz
tabix -p vcf $(echo ${POPULATION[*]%.*} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.hg38.vcf.gz
printf "\n"
