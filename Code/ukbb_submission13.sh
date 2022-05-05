#! /bin/bash
set -e

POPULATION=$1 # comma-separated list of case,control populations with .indiv files
CHR=$2 # single chromosome number
BP_RANGE=$3 # comma-separated list from_bp,to_bp of regions on chromosome CHR
HAPLOTYPE=$4 # comma-separated list of rsid_allele,... to condition on
SNPS=$5 # SNPs to be included, supplied as a file with one rsid per line
VCF_FILES=$6 # comma-separated list of full path VCF file to extract haplotypes
INFO=$7 # snp.info.txt, used for converting hg38 to hg19
ALPHA=$8 # p value cutoff for selecting patterns to carry forward
N_SAMPLES=$9 # how many rows of the haplotypes file to select for initial overlaps
COMPLEMENT=${10} # whether to take complement of sampled set
SEED=${11} # random seed for drawing rows from haplotypes file
DIRECTORY=${12} # current working directory
HOME_DIR=${13} # location of program files

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
echo Directory is $DIRECTORY
cd $DIRECTORY
printf "\n"

module load python/3.7.0
module load python/conda/3.8.1
module load conda3/5.1.0

POPULATION=($(echo $POPULATION | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))

for ((i=0;i<${#POPULATION[@]};i++))
do
  test -f ${POPULATION[i]}.indiv || (>&2 echo "Population $i file not found."; exit 1)
done

BP_RANGE=($(echo $BP_RANGE | perl -pne 's/([0-9]+)[,]*/$1 /g'))

if [ -z $SEED ]
then
  SEED=20200116
fi

start_time=$(date +%s)

queue="large_mem" # specify queue for large-memory tasks
mem_req=20000 # memory requested for large-memory tasks (in mb)

echo Extract haplotypes:
echo bsub \-P SJLIFE \-J ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo ${DIRECTORY}/ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY}/ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-R \"rusage[mem=$mem_req]\" \-q $queue \-K \"sh ${HOME_DIR}/ukbb_haplotype_extract3.sh $(echo ${POPULATION[*]/%/.indiv} | sed 's/ /,/g') \\\"${SNPS}\\\" \\\"${INFO}\\\" ${CHR} ${BP_RANGE[0]},${BP_RANGE[1]} \\\"${VCF_FILES}\\\" ${DIRECTORY} ${HOME_DIR}\"
bsub -P SJLIFE -J ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo ${DIRECTORY}/ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY}/ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -R "rusage[mem=$mem_req]" -q $queue -K "sh ${HOME_DIR}/ukbb_haplotype_extract3.sh $(echo ${POPULATION[*]/%/.indiv} | sed 's/ /,/g') \"${SNPS}\" \"${INFO}\" ${CHR} ${BP_RANGE[0]},${BP_RANGE[1]} \"${VCF_FILES}\" ${DIRECTORY} ${HOME_DIR}"
printf "\n"

echo Wait until haplotypes extracted:
echo bsub \-P SJLIFE \-J sleep_0.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo ${DIRECTORY}/sleep_0.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY}/sleep_0.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-w \"done\(ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}\)\" \-R \"rusage[mem=32]\" \-K \"sleep 10\"
bsub -P SJLIFE -J sleep_0.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo ${DIRECTORY}/sleep_0.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY}/sleep_0.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -w "done(ukbb_haplotype_extract3.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]})" -R "rusage[mem=32]" -K "sleep 10"
printf "\n"

# find haplotype_estimates files
test -f haplotype_estimates.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt && CASES_HAPLOTYPES=haplotype_estimates.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
test -f haplotype_estimates.${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt && CONTROLS_HAPLOTYPES=haplotype_estimates.${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
echo $CASES_HAPLOTYPES
test -f $CASES_HAPLOTYPES || (>&2 echo "haplotype_estimates.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt does not exist in ${DIRECTORY}"; exit 1)
echo $CONTROLS_HAPLOTYPES
test -f $CONTROLS_HAPLOTYPES || (>&2 echo "haplotype_estimates.${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt does not exist in ${DIRECTORY}"; exit 1)
printf "\n"

echo Keep only rs SNPs:
echo awk \'NR==1{printf \"%s\",\$1\; for \(i=2\;i\<NF\;i++\) {if \(\$i \~ /rs[0-9]*_[A-Z]\$/\) {included[i]=1\; printf \"\\t%s\",\$i} } if \(\$NF \~ /rs[0-9]*_[A-Z]\$/\) {included[NF]=1\; printf \"\\t%s\\n\",\$NF} else {printf \"\\n\"} } NR\>1{printf \"%s\",\$1\; for \(i=2\;i\<NF\;i++\) {if \(included[i]==1\) {printf \"\\t%s\",\$i} } if \(included[NF]==1\) {printf \"\\t%s\\n\",\$NF} else {printf \"\\n\"} }\' ${CASES_HAPLOTYPES} \> ${CASES_HAPLOTYPES%.*}.tmp
awk 'NR==1{printf "%s",$1; for (i=2;i<NF;i++) {if ($i ~ /rs[0-9]*_[A-Z]$/) {included[i]=1; printf "\t%s",$i} } if ($NF ~ /rs[0-9]*_[A-Z]$/) {included[NF]=1; printf "\t%s\n",$NF} else {printf "\n"} } NR>1{printf "%s",$1; for (i=2;i<NF;i++) {if (included[i]==1) {printf "\t%s",$i} } if (included[NF]==1) {printf "\t%s\n",$NF} else {printf "\n"} }' ${CASES_HAPLOTYPES} > ${CASES_HAPLOTYPES%.*}.tmp
printf "\n"

test -f ${CASES_HAPLOTYPES%.*}.tmp && echo mv ${CASES_HAPLOTYPES%.*}.tmp ${CASES_HAPLOTYPES} && printf "\n"
test -f ${CASES_HAPLOTYPES%.*}.tmp && mv ${CASES_HAPLOTYPES%.*}.tmp ${CASES_HAPLOTYPES}
printf "\n"

echo Included SNPs:
head -1 ${CASES_HAPLOTYPES}
awk 'NR==1{printf "%s SNP%s\n",NF-1,(NF-1==1?"":"s")}' ${CASES_HAPLOTYPES}
printf "\n"

echo awk \'NR==1{printf \"%s\",\$1\; for \(i=2\;i\<NF\;i++\) {if \(\$i \~ /rs[0-9]*_[A-Z]\$/\) {included[i]=1\; printf \"\\t%s\",\$i} } if \(\$NF \~ /rs[0-9]*_[A-Z]\$/\) {included[NF]=1\; printf \"\\t%s\\n\",\$NF} else {printf \"\\n\"} } NR\>1{printf \"%s\",\$1\; for \(i=2\;i\<NF\;i++\) {if \(included[i]==1\) {printf \"\\t%s\",\$i} } if \(included[NF]==1\) {printf \"\\t%s\\n\",\$NF} else {printf \"\\n\"} }\' ${CONTROLS_HAPLOTYPES} \> ${CONTROLS_HAPLOTYPES%.*}.tmp
awk 'NR==1{printf "%s",$1; for (i=2;i<NF;i++) {if ($i ~ /rs[0-9]*_[A-Z]$/) {included[i]=1; printf "\t%s",$i} } if ($NF ~ /rs[0-9]*_[A-Z]$/) {included[NF]=1; printf "\t%s\n",$NF} else {printf "\n"} } NR>1{printf "%s",$1; for (i=2;i<NF;i++) {if (included[i]==1) {printf "\t%s",$i} } if (included[NF]==1) {printf "\t%s\n",$NF} else {printf "\n"} }' ${CONTROLS_HAPLOTYPES} > ${CONTROLS_HAPLOTYPES%.*}.tmp
printf "\n"

test -f ${CONTROLS_HAPLOTYPES%.*}.tmp && echo mv ${CONTROLS_HAPLOTYPES%.*}.tmp ${CONTROLS_HAPLOTYPES} && printf "\n"
test -f ${CONTROLS_HAPLOTYPES%.*}.tmp && mv ${CONTROLS_HAPLOTYPES%.*}.tmp ${CONTROLS_HAPLOTYPES}

if [ ! -z $HAPLOTYPE ] && [ -f $(echo ${POPULATION[*]} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.haplotypes.hg38.vcf ]
then
  if [ ! -z $SNPS ]
  then
    SNPS=$(awk 'NR==FNR{snp[$1]; next} ($3 in snp && "'${HAPLOTYPE}'" !~ $3){printf "%s%s_%s",(found>0?",":""),$3,$5; found+=1} END{printf "\n"}' ${SNPS} $(echo ${POPULATION[*]} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.haplotypes.hg38.vcf) # create list rsid_allele from snps in SNPS file but not in HAPLOTYPE list
    echo Select SNPs $SNPS
    printf "\n"
  else
    echo No subset SNPs selected.
    printf "\n"
  fi

  echo Subset haplotype_estimates files based on haplotype $HAPLOTYPE
  echo sh ${HOME_DIR}/ukbb_conditional_haplotype_subset.sh ${CASES_HAPLOTYPES},${CONTROLS_HAPLOTYPES} ${HAPLOTYPE} ${SNPS}
  sh ${HOME_DIR}/ukbb_conditional_haplotype_subset.sh ${CASES_HAPLOTYPES},${CONTROLS_HAPLOTYPES} ${HAPLOTYPE} ${SNPS}
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
echo bsub \-P SJLIFE \-J haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo ${DIRECTORY}/haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY}/haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-R \"rusage[mem=$mem_req]\" \-q $queue \-K \"sh ${HOME_DIR}/ukbb_haplotype_transpose.sh ${CASES_HAPLOTYPES} 0 ${DIRECTORY}\"
bsub -P SJLIFE -J haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo ${DIRECTORY}/haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY}/haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -R "rusage[mem=$mem_req]" -q $queue -K "sh ${HOME_DIR}/ukbb_haplotype_transpose.sh ${CASES_HAPLOTYPES} 0 ${DIRECTORY}"
printf "\n"

echo Wait until cases haplotypes transposed:
echo bsub \-P SJLIFE \-J sleep_1.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-w \"done\(haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}\)\" \-R \"rusage[mem=32]\" \-oo ${DIRECTORY}/sleep_1.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY}/sleep_1.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-K \"sleep 10\"
bsub -P SJLIFE -J sleep_1.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -w "done(haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]})" -R "rusage[mem=32]" -oo ${DIRECTORY}/sleep_1.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY}/sleep_1.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -K "sleep 10"
printf "\n"

echo Transpose controls haplotypes file, no header:
echo bsub \-P SJLIFE \-J haplotypes_transpose.${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo ${DIRECTORY}/haplotypes_transpose.${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY}/haplotypes_transpose.${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-R \"rusage[mem=$mem_req]\" \-q $queue \-K \"sh ${HOME_DIR}/ukbb_haplotype_transpose.sh ${CONTROLS_HAPLOTYPES} 0 ${DIRECTORY}\"
bsub -P SJLIFE -J haplotypes_transpose.${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo ${DIRECTORY}/haplotypes_transpose.${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY}/haplotypes_transpose.${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -R "rusage[mem=$mem_req]" -q $queue -K "sh ${HOME_DIR}/ukbb_haplotype_transpose.sh ${CONTROLS_HAPLOTYPES} 0 ${DIRECTORY}"
printf "\n"

echo Wait until controls haplotypes transposed:
echo bsub \-P SJLIFE \-J sleep_2.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-w \"done\(haplotypes_transpose.${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}\)\" \-R \"rusage[mem=32]\" \-oo ${DIRECTORY}/sleep_2.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY}/sleep_2.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-K \"sleep 10\"
bsub -P SJLIFE -J sleep_2.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -w "done(haplotypes_transpose.${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]})" -R "rusage[mem=32]" -oo ${DIRECTORY}/sleep_2.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY}/sleep_2.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -K "sleep 10"
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
  echo python ${HOME_DIR}/ukbb_row_sample.py \-n $N_SAMPLES \-s $SEED \-f ${CASES_HAPLOTYPES}
  python ${HOME_DIR}/ukbb_row_sample.py -n $N_SAMPLES -s $SEED -f ${CASES_HAPLOTYPES} # sample rows using random seed
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
  echo bsub \-P SJLIFE \-J haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${ROWS} \-oo ${DIRECTORY}/haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${ROWS}.out \-eo ${DIRECTORY}/haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${ROWS}.err \-R \"rusage[mem=$mem_req]\" \-q $queue \-K \"sh ${HOME_DIR}/ukbb_haplotype_transpose.sh ${CASES_HAPLOTYPES%.*}${ROWS}.txt 0 ${DIRECTORY}\"
  bsub -P SJLIFE -J haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${ROWS} -oo ${DIRECTORY}/haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${ROWS}.out -eo ${DIRECTORY}/haplotypes_transpose.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}${ROWS}.err -R "rusage[mem=$mem_req]" -q $queue -K "sh ${HOME_DIR}/ukbb_haplotype_transpose.sh ${CASES_HAPLOTYPES%.*}${ROWS}.txt 0 ${DIRECTORY}"
  printf "\n"

  echo Examine sampled haplotypes file:
  echo head ${CASES_HAPLOTYPES%.*}${ROWS}.txt
  head ${CASES_HAPLOTYPES%.*}${ROWS}.txt
  echo $(cat ${CASES_HAPLOTYPES%.*}${ROWS}.txt | wc -l) line$( (( $(cat ${CASES_HAPLOTYPES%.*}${ROWS}.txt | wc -l) > 1 )) && echo "s" || echo "")
  printf "\n"
fi

echo Initial round of overlaps of all cases ${POPULATION[0]}.chr${CHR}:
echo bsub \-P SJLIFE \-J subjects_bsub2.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo ${DIRECTORY}/subjects_bsub2.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY}/subjects_bsub2.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-R \"rusage[mem=512]\" \-K \"sh ${HOME_DIR}/subjects_bsub2.sh ${CASES_HAPLOTYPES%.*}${ROWS}.txt $CHR ${BP_RANGE[0]},${BP_RANGE[1]} 1 100 50 ${POPULATION[0]} ${DIRECTORY} ${HOME_DIR}\"
bsub -P SJLIFE -J subjects_bsub2.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo ${DIRECTORY}/subjects_bsub2.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY}/subjects_bsub2.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -R "rusage[mem=512]" -K "sh ${HOME_DIR}/subjects_bsub2.sh ${CASES_HAPLOTYPES%.*}${ROWS}.txt $CHR ${BP_RANGE[0]},${BP_RANGE[1]} 1 100 50 ${POPULATION[0]} ${DIRECTORY} ${HOME_DIR}"
printf "\n"

echo Combine results from Iteration000:
echo bsub \-P SJLIFE \-J combine_Iteration000.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo combine_Iteration000.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo combine_Iteration000.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -R \"rusage[mem=10000]\" -q large_mem \-K \"sh ${HOME_DIR}/subjects_combine_tuples.sh ${POPULATION[0]} $CHR ${BP_RANGE[0]},${BP_RANGE[1]} 1 \\\"Iteration000\\\" $DIRECTORY $HOME_DIR\"
bsub -P SJLIFE -J combine_Iteration000.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo combine_Iteration000.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo combine_Iteration000.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -R "rusage[mem=10000]" -q large_mem -K "sh ${HOME_DIR}/subjects_combine_tuples.sh ${POPULATION[0]} $CHR ${BP_RANGE[0]},${BP_RANGE[1]} 1 \"Iteration000\" $DIRECTORY $HOME_DIR"
printf "\n"

echo Wait until combinations complete:
echo bsub \-P SJLIFE \-J sleep_3.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo ${DIRECTORY}/sleep_3.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY}/sleep_3.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-w \"done\(combine_Iteration000.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}\)\" \-R \"rusage[mem=32]\" \-K \"sleep 10\"
bsub -P SJLIFE -J sleep_3.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo ${DIRECTORY}/sleep_3.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY}/sleep_3.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -w "done(combine_Iteration000.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]})" -R "rusage[mem=32]" -K "sleep 10"
echo Done.
printf "\n"

combined_patterns_file=$(for file in Pattern_combined.Iteration000*
do
  echo $file | awk '!($0 ~ /[+]/){print $0}'
done | grep ${POPULATION[0]}) # get file with one bar group corresponding to cases only

test -z $ALPHA && ALPHA=0.00000001 # p value cutoff for carrying patterns forward

if [ -f $combined_patterns_file ];
then
  echo Evaluate pattern frequency differences:
  echo bsub \-P SJLIFE \-J pattern_difference.Iteration000.${POPULATION[0]}+${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo ${DIRECTORY}/pattern_difference.Iteration000.${POPULATION[0]}+${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY}/pattern_difference.Iteration000.${POPULATION[0]}+${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-R \"rusage[mem=10000]\" \-q large_mem \-K \"${HOME_DIR}/ukbb_pattern_difference_sub.sh ${CASES_HAPLOTYPES},${CONTROLS_HAPLOTYPES} $combined_patterns_file 50 \\\"Iteration000\\\" $ALPHA 0\"
  bsub -P SJLIFE -J pattern_difference.Iteration000.${POPULATION[0]}+${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo ${DIRECTORY}/pattern_difference.Iteration000.${POPULATION[0]}+${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY}/pattern_difference.Iteration000.${POPULATION[0]}+${POPULATION[1]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -R "rusage[mem=10000]" -q large_mem -K "${HOME_DIR}/ukbb_pattern_difference_sub.sh ${CASES_HAPLOTYPES},${CONTROLS_HAPLOTYPES} $combined_patterns_file 50 \"Iteration000\" $ALPHA 0" # filter patterns with p value < $ALPHA and remove ubiquitous alleles appearing in all patterns
  printf "\n"
fi

echo Begin iterations:
echo bsub \-P SJLIFE \-J pattern_overlap_iterate3.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo ${DIRECTORY}/pattern_overlap_iterate3.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY}/pattern_overlap_iterate3.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-R \"rusage[mem=512]\" \"sh ${HOME_DIR}/pattern_overlap_iterate3.sh ${CASES_HAPLOTYPES},${CONTROLS_HAPLOTYPES} $CHR ${BP_RANGE[0]},${BP_RANGE[1]} 100 50 2 2,j \\\"\\\" \\\"${POPULATION[0]}\\\" $DIRECTORY $HOME_DIR\"
# bsub -P SJLIFE -J pattern_overlap_iterate3.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo ${DIRECTORY}/pattern_overlap_iterate3.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY}/pattern_overlap_iterate3.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -R "rusage[mem=512]" -K "sh ${HOME_DIR}/pattern_overlap_iterate3.sh ${CASES_HAPLOTYPES},${CONTROLS_HAPLOTYPES} $CHR ${BP_RANGE[0]},${BP_RANGE[1]} 100 50 2 2,j \"\" \"${POPULATION[0]}\" $DIRECTORY $HOME_DIR"
# printf "\n"
#
# echo Wait until iterations have completed:
# echo bsub \-P SJLIFE \-J sleep_4.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.overlap \-oo ${DIRECTORY}/sleep_4.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.overlap.out \-eo ${DIRECTORY}/sleep_4.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.overlap.err \-w \"done\(pattern_overlap_iterate3.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}\)\" \-R \"rusage[mem=32]\" \-K \"sleep 10\"
# bsub -P SJLIFE -J sleep_4.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.overlap -oo ${DIRECTORY}/sleep_4.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.overlap.out -eo ${DIRECTORY}/sleep_4.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.overlap.err -w "done(pattern_overlap_iterate3.${POPULATION[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]})" -R "rusage[mem=32]" -K "sleep 10"
# printf "\n"
