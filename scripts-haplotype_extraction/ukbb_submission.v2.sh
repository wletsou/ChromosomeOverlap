#! /bin/bash
set -e

DISCOVERY=$1 # comma-separated list of case,control discovery populations with .indiv files
REPLICATION=$2 # comma-separated list of case,control replication populations with .indiv files
VALIDATION=$3 # comma-separated list of case,control validation populations with .indiv files
CHR=$4 # single chromosome number
BP_RANGE=$5 # comma-separated list from_bp,to_bp of regions on chromosome CHR
HAPLOTYPE=$6 # comma-separated list of rsid_allele,... to condition on
VALUE=$7 # 0 or 1: do we condition on non-carriers or carriers (default)?
SNPS=$8 # SNPs to be included, supplied as a file with one rsid per line
VCF_FILES=$9 # comma-separated list of full path VCF file to extract haplotypes
INFO=${10} # snp.info.txt, used for converting hg38 to hg19
N_SAMPLES=${11} # how many rows of the haplotypes file to select for initial overlaps
COMPLEMENT=${12} # whether to take complement of sampled set
SEED=${13} # random seed for drawing rows from haplotypes file
DIRECTORY=${14} # current working directory
HOME_DIR=${15} # location of program files

module load R/3.6.1
module load python/3.7.0
module load python/conda/3.8.1
module load conda3/5.1.0

if [ -z $HOME_DIR ];
then
  HOME_DIR="/home/wletsou/scripts" # in case HOME_DIR is on your PATH
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
echo Directory is $DIRECTORY
cd $DIRECTORY
printf "\n"

DISCOVERY=($(echo $DISCOVERY | sed -E 's/([^,]+)[,]*/\1 /g'))

for ((i=0;i<${#DISCOVERY[@]};i++))
do
  test -f ${DISCOVERY[i]}.indiv || (>&2 echo "Discovery $i file not found."; exit 1)
done
(( ${#DISCOVERY[@]}==2 )) || (>&2 echo "Must supply discovery cases and controls."; exit 1)

REPLICATION=($(echo $REPLICATION | sed -E 's/([^,]+)[,]*/\1 /g'))

for ((i=0;i<${#REPLICATION[@]};i++))
do
  test -f ${REPLICATION[i]}.indiv || (>&2 echo "Replication $i file not found."; exit 1)
done

VALIDATION=($(echo $VALIDATION | sed -E 's/([^,]+)[,]*/\1 /g'))

for ((i=0;i<${#VALIDATION[@]};i++))
do
  test -f ${VALIDATION[i]}.indiv || (>&2 echo "Validation $i file not found."; exit 1)
done

POPULATION=(${DISCOVERY[*]} ${REPLICATION[*]} ${VALIDATION[*]}) # up to six populations

BP_RANGE=($(echo $BP_RANGE | sed -E 's/([0-9]+)[,]*/\1 /g'))

if [ -z $SEED ]
then
  SEED=20200116
fi

start_time=$(date +%s)

queue="standard" # specify queue for large-memory tasks
mem_req=50000 # memory requested for large-memory tasks (in mb)

echo Extract haplotypes:
echo bsub \-P SJLIFE \-J ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo ${DIRECTORY/%/\/}ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY/%/\/}ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-R \"rusage[mem=$mem_req]\" \-q $queue \-K \"sh ${HOME_DIR/%/\/}ukbb_haplotype_extract3.sh $(echo ${POPULATION[*]/%/.indiv} | sed 's/ /,/g') \\\"${SNPS}\\\" \\\"${INFO}\\\" ${CHR} ${BP_RANGE[0]},${BP_RANGE[1]} \\\"${VCF_FILES}\\\" ${DIRECTORY} ${HOME_DIR}\"
job_id=$(bsub -P SJLIFE -J ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo ${DIRECTORY/%/\/}ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY/%/\/}ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -R "rusage[mem=$mem_req]" -q $queue -K "sh ${HOME_DIR/%/\/}ukbb_haplotype_extract3.sh $(echo ${POPULATION[*]/%/.indiv} | sed 's/ /,/g') \"${SNPS}\" \"${INFO}\" ${CHR} ${BP_RANGE[0]},${BP_RANGE[1]} \"${VCF_FILES}\" ${DIRECTORY} ${HOME_DIR}")
echo $job_id
job_id=$(echo $job_id | awk 'b=gensub(/.*<([0-9]*)>.*/,"\\1","g",$0) {print b}') #extract job_id (number) from output
printf "\n"

echo Wait until haplotypes extracted:
echo bsub \-P SJLIFE \-J ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo ${DIRECTORY/%/\/}sleep.ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY/%/\/}sleep.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-w \"done\(${job_id}\)\" \-R \"rusage[mem=32]\" \-K \"sleep 10\"
bsub -P SJLIFE -J sleep.ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo ${DIRECTORY/%/\/}sleep.ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY/%/\/}sleep.ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -w "done(${job_id})" -R "rusage[mem=32]" -K "sleep 10"
printf "\n"

if [ ! -z $REPLICATION ] || [ ! -z $VALIDATION ]
then
  echo Make combined cases+controls files:
  echo awk \'BEGIN{OFS=\"\\t\"} {print \$0}\' $(echo ${DISCOVERY[*]/%/.indiv}) \> $(echo ${DISCOVERY[*]} | sed -E 's/[ ]+/+/g').indiv
  awk 'BEGIN{OFS="\t"} {print $0}' $(echo ${DISCOVERY[*]/%/.indiv}) > $(echo ${DISCOVERY[*]} | sed -E 's/[ ]+/+/g').indiv
  echo awk \'NR==FNR{id[\$1]\; next} \(NR!=FNR \&\& FNR==1\){print \$0} \(NR!=FNR \&\& FNR\>1 \&\& \$1 in id\){print \$0}\' $(echo ${DISCOVERY[*]} | sed -E 's/[ ]+/+/g').indiv haplotype_estimates.$(echo ${POPULATION[*]}| sed -E 's/[ ]+/+/g').chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt \> haplotype_estimates.$(echo ${DISCOVERY[*]} | sed -E 's/[ ]+/+/g').chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
  awk 'NR==FNR{id[$1]; next} (NR!=FNR && FNR==1){print $0} (NR!=FNR && FNR>1 && $1 in id){print $0}' $(echo ${DISCOVERY[*]} | sed -E 's/[ ]+/+/g').indiv haplotype_estimates.$(echo ${POPULATION[*]}| sed -E 's/[ ]+/+/g').chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt > haplotype_estimates.$(echo ${DISCOVERY[*]} | sed -E 's/[ ]+/+/g').chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt && printf "\n"
fi

if [ ! -z $REPLICATION ]
then
  echo awk \'BEGIN{OFS=\"\\t\"} {print \$0}\' $(echo ${REPLICATION[*]/%/.indiv}) \> $(echo ${REPLICATION[*]} | sed -E 's/[ ]+/+/g').indiv
  awk 'BEGIN{OFS="\t"} {print $0}' $(echo ${REPLICATION[*]/%/.indiv}) > $(echo ${REPLICATION[*]} | sed -E 's/[ ]+/+/g').indiv
  echo awk \'NR==FNR{id[\$1]\; next} \(NR!=FNR \&\& FNR==1\){print \$0} \(NR!=FNR \&\& FNR\>1 \&\& \$1 in id\){print \$0}\' $(echo ${REPLICATION[*]} | sed -E 's/[ ]+/+/g').indiv haplotype_estimates.$(echo ${POPULATION[*]}| sed -E 's/[ ]+/+/g').chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt \> haplotype_estimates.$(echo ${REPLICATION[*]} | sed -E 's/[ ]+/+/g').chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
  awk 'NR==FNR{id[$1]; next} (NR!=FNR && FNR==1){print $0} (NR!=FNR && FNR>1 && $1 in id){print $0}' $(echo ${REPLICATION[*]} | sed -E 's/[ ]+/+/g').indiv haplotype_estimates.$(echo ${POPULATION[*]}| sed -E 's/[ ]+/+/g').chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt > haplotype_estimates.$(echo ${REPLICATION[*]} | sed -E 's/[ ]+/+/g').chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt && printf "\n"
fi

if [ ! -z $VALIDATION ]
then
  echo awk \'BEGIN{OFS=\"\\t\"} {print \$0}\' $(echo ${VALIDATION[*]/%/.indiv}) \> $(echo ${VALIDATION[*]} | sed -E 's/[ ]+/+/g').indiv
  awk 'BEGIN{OFS="\t"} {print $0}' $(echo ${VALIDATION[*]/%/.indiv}) > $(echo ${VALIDATION[*]} | sed -E 's/[ ]+/+/g').indiv
  echo awk \'NR==FNR{id[\$1]\; next} \(NR!=FNR \&\& FNR==1\){print \$0} \(NR!=FNR \&\& FNR\>1 \&\& \$1 in id\){print \$0}\' $(echo ${VALIDATION[*]} | sed -E 's/[ ]+/+/g').indiv haplotype_estimates.$(echo ${POPULATION[*]}| sed -E 's/[ ]+/+/g').chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt \> haplotype_estimates.$(echo ${VALIDATION[*]} | sed -E 's/[ ]+/+/g').chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
  awk 'NR==FNR{id[$1]; next} (NR!=FNR && FNR==1){print $0} (NR!=FNR && FNR>1 && $1 in id){print $0}' $(echo ${VALIDATION[*]} | sed -E 's/[ ]+/+/g').indiv haplotype_estimates.$(echo ${POPULATION[*]}| sed -E 's/[ ]+/+/g').chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt > haplotype_estimates.$(echo ${VALIDATION[*]} | sed -E 's/[ ]+/+/g').chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt && printf "\n"
fi

# find haplotype_estimates files
test -f haplotype_estimates.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt && CASES_HAPLOTYPES=haplotype_estimates.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
test -f haplotype_estimates.${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt && CONTROLS_HAPLOTYPES=haplotype_estimates.${DISCOVERY[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
echo $CASES_HAPLOTYPES
test -f $CASES_HAPLOTYPES || (>&2 echo "haplotype_estimates.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt does not exist in ${DIRECTORY}"; exit 1)
echo $CONTROLS_HAPLOTYPES
test -f $CONTROLS_HAPLOTYPES || (>&2 echo "haplotype_estimates.${DISCOVERY[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt does not exist in ${DIRECTORY}"; exit 1)
printf "\n"

echo Included SNPs:
head -1 ${CASES_HAPLOTYPES}
awk 'NR==1{printf "%s SNP%s\n",NF-1,(NF-1==1?"":"s")}' ${CASES_HAPLOTYPES}
printf "\n"

echo Choose only rs SNPs:
echo awk \'NR==1{printf \"%s\",\$1\; for \(i=2\;i\<NF\;i++\) {if \(\$i \~ /rs[0-9]*_[A-Z]\$/\) {included[i]=1\; printf \"\\t%s\",\$i} } if \(\$NF \~ /rs[0-9]*_[A-Z]\$/\) {included[NF]=1\; printf \"\\t%s\\n\",\$NF} else {printf \"\\n\"} } NR\>1{printf \"%s\",\$1\; for \(i=2\;i\<NF\;i++\) {if \(included[i]==1\) {printf \"\\t%s\",\$i} } if \(included[NF]==1\) {printf \"\\t%s\\n\",\$NF} else {printf \"\\n\"} }\' ${CONTROLS_HAPLOTYPES} \> ${CONTROLS_HAPLOTYPES%.*}.tmp
awk 'NR==1{printf "%s",$1; for (i=2;i<NF;i++) {if ($i ~ /rs[0-9]*_[A-Z]$/) {included[i]=1; printf "\t%s",$i} } if ($NF ~ /rs[0-9]*_[A-Z]$/) {included[NF]=1; printf "\t%s\n",$NF} else {printf "\n"} } NR>1{printf "%s",$1; for (i=2;i<NF;i++) {if (included[i]==1) {printf "\t%s",$i} } if (included[NF]==1) {printf "\t%s\n",$NF} else {printf "\n"} }' ${CONTROLS_HAPLOTYPES} > ${CONTROLS_HAPLOTYPES%.*}.tmp
printf "\n"

test -f ${CONTROLS_HAPLOTYPES%.*}.tmp && echo mv ${CONTROLS_HAPLOTYPES%.*}.tmp ${CONTROLS_HAPLOTYPES} && printf "\n"
test -f ${CONTROLS_HAPLOTYPES%.*}.tmp && mv ${CONTROLS_HAPLOTYPES%.*}.tmp ${CONTROLS_HAPLOTYPES}

if [ ! -z $HAPLOTYPE ]
then
  if [ -z $VALUE ]
  then
    VALUE=1 # defualt is to condition on carrying the the haplotype
  fi
  echo Subset haplotype_estimates files based on haplotype $HAPLOTYPE
  echo sh ${HOME_DIR/%/\/}conditional_haplotype_subset.sh ${CASES_HAPLOTYPES},${CONTROLS_HAPLOTYPES} ${HAPLOTYPE} ${VALUE} ""
  sh ${HOME_DIR/%/\/}conditional_haplotype_subset.sh ${CASES_HAPLOTYPES},${CONTROLS_HAPLOTYPES} ${HAPLOTYPE} ${VALUE} ""
  printf "\n"

  test -f ${CASES_HAPLOTYPES%.*txt}.subset.txt || (>&2 echo "${CASES_HAPLOTYPES%.*txt}.subset.txt does not exist in ${DIRECTORY}"; exit 1)
  test -f ${CONTROLS_HAPLOTYPES%.*txt}.subset.txt || (>&2 echo "${CONTROLS_HAPLOTYPES%.*txt}.subset.txt does not exist in ${DIRECTORY}"; exit 1)
  printf "\n"

  echo Replace HAPLOTYPES variables
  test -f ${CASES_HAPLOTYPES%.*txt}.subset.txt && echo CASES_HAPLOTYPES=${CASES_HAPLOTYPES%.*txt}.subset.txt
  test -f ${CASES_HAPLOTYPES%.*txt}.subset.txt && CASES_HAPLOTYPES=${CASES_HAPLOTYPES%.*txt}.subset.txt
  echo $CASES_HAPLOTYPES
  test -f ${CONTROLS_HAPLOTYPES%.*txt}.subset.txt && echo CONTROLS_HAPLOTYPES=${CONTROLS_HAPLOTYPES%.*txt}.subset.txt
  test -f ${CONTROLS_HAPLOTYPES%.*txt}.subset.txt && CONTROLS_HAPLOTYPES=${CONTROLS_HAPLOTYPES%.*txt}.subset.txt
  echo $CONTROLS_HAPLOTYPES
  printf "\n"
fi

echo Transpose cases haplotypes file, no header:
echo bsub \-P SJLIFE \-J haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo ${DIRECTORY/%/\/}haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY/%/\/}haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-R \"rusage[mem=$mem_req]\" \-q $queue \-K \"sh ${HOME_DIR/%/\/}ukbb_haplotype_transpose.sh ${CASES_HAPLOTYPES} 0 ${DIRECTORY}\"
job_id=$(bsub -P SJLIFE -J haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo ${DIRECTORY/%/\/}haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY/%/\/}haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -R "rusage[mem=$mem_req]" -q $queue -K "sh ${HOME_DIR/%/\/}ukbb_haplotype_transpose.sh ${CASES_HAPLOTYPES} 0 ${DIRECTORY}")
echo $job_id
job_id=$(echo $job_id | awk 'b=gensub(/.*<([0-9]*)>.*/,"\\1","g",$0) {print b}') #extract job_id (number) from output

printf "\n"

echo Wait until cases haplotypes transposed:
echo bsub \-P SJLIFE \-J sleep.haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-w \"done\(${job_id}\)\" \-R \"rusage[mem=32]\" \-oo ${DIRECTORY/%/\/}sleep.haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY/%/\/}sleep.haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-K \"sleep 10\"
bsub -P SJLIFE -J sleep.haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -w "done(${job_id})" -R "rusage[mem=32]" -oo ${DIRECTORY/%/\/}sleep.haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY/%/\/}sleep.haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -K "sleep 10"
printf "\n"

echo Transpose controls haplotypes file, no header:
echo bsub \-P SJLIFE \-J haplotypes_transpose.${DISCOVERY[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo ${DIRECTORY/%/\/}haplotypes_transpose.${DISCOVERY[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY/%/\/}haplotypes_transpose.${DISCOVERY[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-R \"rusage[mem=$mem_req]\" \-q $queue \-K \"sh ${HOME_DIR/%/\/}ukbb_haplotype_transpose.sh ${CONTROLS_HAPLOTYPES} 0 ${DIRECTORY}\"
job_id=$(bsub -P SJLIFE -J haplotypes_trasnpose.${DISCOVERY[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo ${DIRECTORY/%/\/}haplotypes_transpose.${DISCOVERY[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY/%/\/}haplotypes_transpose.${DISCOVERY[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -R "rusage[mem=$mem_req]" -q $queue -K "sh ${HOME_DIR/%/\/}ukbb_haplotype_transpose.sh ${CONTROLS_HAPLOTYPES} 0 ${DIRECTORY}")
echo $job_id
job_id=$(echo $job_id | awk 'b=gensub(/.*<([0-9]*)>.*/,"\\1","g",$0) {print b}') #extract job_id (number) from output
printf "\n"

echo Wait until controls haplotypes transposed:
echo bsub \-P SJLIFE \-J sleep.haplotypes_transpose.${DISCOVERY[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-w \"done\(${job_id}\)\" \-R \"rusage[mem=32]\" \-oo ${DIRECTORY/%/\/}sleep.haplotypes_transpose.${DISCOVERY[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY/%/\/}sleep.haplotypes_transpose.${DISCOVERY[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-K \"sleep 10\"
bsub -P SJLIFE -J sleep.haplotypes_transpose.${DISCOVERY[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -w "done(${job_id})" -R "rusage[mem=32]" -oo ${DIRECTORY/%/\/}sleep.haplotypes_transpose.${DISCOVERY[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY/%/\/}sleep.haplotypes_transpose.${DISCOVERY[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -K "sleep 10"
printf "\n"

echo Examine haplotypes files:
echo head ${CASES_HAPLOTYPES}
head ${CASES_HAPLOTYPES}
echo $(cat ${CASES_HAPLOTYPES} | wc -l) line$( (( $(cat ${CASES_HAPLOTYPES} | wc -l) > 1 )) && echo "s" || echo "")
printf "\n"

echo head ${CONTROLS_HAPLOTYPES}
head ${CONTROLS_HAPLOTYPES}
echo $(cat ${CONTROLS_HAPLOTYPES} | wc -l) line$( (( $(cat ${CONTROLS_HAPLOTYPES} | wc -l) > 1 )) && echo "s" || echo "")
printf "\n"

if [ ! -z $N_SAMPLES ]
then
  echo python ${HOME_DIR/%/\/}ukbb_row_sample.py \-n $N_SAMPLES \-s $SEED \-f ${CASES_HAPLOTYPES}
  python ${HOME_DIR/%/\/}ukbb_row_sample.py -n $N_SAMPLES -s $SEED -f ${CASES_HAPLOTYPES} # sample rows using random seed
  ROWS=".row_sample_${N_SAMPLES}" # for renaming row-sampled file
  printf "\n"

  if [ ! -z $COMPLMENT ]
  then
    COMPLEMENT=0 # whether to take complement of sampled rows; default is not to
  fi
  if (($COMPLEMENT==1))
  then
    echo Get complement set:
    echo awk \'\(NR==FNR \&\& NR\>1\){seen[\$0]=\$0\; next} \(FNR!=NR \&\& FNR==1\){print \$0} \(\$0 in seen==0\){print \$0}\' ${CASES_HAPLOTYPES%.*}${ROWS}.txt ${CASES_HAPLOTYPES} \> ${CASES_HAPLOTYPES%.*}${ROWS}.tmp
    awk '(NR==FNR && NR>1){seen[$0]=$0; next} (FNR!=NR && FNR==1){print $0} ($0 in seen==0){print $0}' ${CASES_HAPLOTYPES%.*}${ROWS}.txt ${CASES_HAPLOTYPES} > ${CASES_HAPLOTYPES%.*}${ROWS}.tmp
    test -f ${CASES_HAPLOTYPES%.*}${ROWS}.tmp && echo mv ${CASES_HAPLOTYPES%.*}${ROWS}.tmp ${CASES_HAPLOTYPES%.*}${ROWS}.txt && printf "\n"
    test -f ${CASES_HAPLOTYPES%.*}${ROWS}.tmp && mv ${CASES_HAPLOTYPES%.*}${ROWS}.tmp ${CASES_HAPLOTYPES%.*}${ROWS}.txt
  fi

  echo Transpose row-sampled haplotypes file, no header:
  echo bsub \-P SJLIFE \-J haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${ROWS} \-oo ${DIRECTORY/%/\/}haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${ROWS}.out \-eo ${DIRECTORY/%/\/}haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${ROWS}.err \-R \"rusage[mem=$mem_req]\" \-q $queue \-K \"sh ${HOME_DIR/%/\/}ukbb_haplotype_transpose.sh ${CASES_HAPLOTYPES%.*}${ROWS}.txt 0 ${DIRECTORY}\"
  bsub -P SJLIFE -J haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${ROWS} -oo ${DIRECTORY/%/\/}haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${ROWS}.out -eo ${DIRECTORY/%/\/}haplotypes_transpose.${DISCOVERY[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${ROWS}.err -R "rusage[mem=$mem_req]" -q $queue -K "sh ${HOME_DIR/%/\/}ukbb_haplotype_transpose.sh ${CASES_HAPLOTYPES%.*}${ROWS}.txt 0 ${DIRECTORY}"
  printf "\n"

  echo Examine sampled haplotypes file:
  echo head ${CASES_HAPLOTYPES%.*}${ROWS}.txt
  head ${CASES_HAPLOTYPES%.*}${ROWS}.txt
  echo $(cat ${CASES_HAPLOTYPES%.*}${ROWS}.txt | wc -l) line$( (( $(cat ${CASES_HAPLOTYPES%.*}${ROWS}.txt | wc -l) > 1 )) && echo "s" || echo "")
  printf "\n"
fi
