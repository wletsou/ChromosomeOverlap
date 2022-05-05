#! /bin/bash
set -e

# bsub -P SJLIFE -J pattern_differences_sub.Iteration001.ukbb_bca_20200116_cases.chr11.98738371-98987273.REP_0 -oo /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.6/REP_0/pattern_differences_sub.Iteration001.ukbb_bca_20200116_cases.chr11.98738371-98987273.REP_0.out -eo /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.6/REP_0/pattern_differences_sub.Iteration001.ukbb_bca_20200116_cases.chr11.98738371-98987273.REP_0.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_pattern_difference_sub.sh haplotype_estimates.ukbb_bca_20200116_cases.chr11.98738371-98987273.txt,haplotype_estimates.ukbb_bca_20200116_controls.chr11.98738371-98987273.txt Pattern_combined.Iteration001.ukbb_bca_20200116_cases.chr11.98738371-98987273_2,j.txt 50 Iteration001"

# bsub -P SJLIFE -J pattern_differences_sub.Iteration000.ukbb_bca_cases+ukbb_bca_controls.chr11.68700000-69700000 -oo /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.10/pattern_differences_sub.Iteration000.ukbb_bca_cases+ukbb_bca_controls.chr11.68700000-69700000.out -eo /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.10/pattern_differences_sub.Iteration000.ukbb_bca_cases+ukbb_bca_controls.chr11.68700000-69700000.err -R "rusage[mem=128]" "sh /home/wletsou/scripts/ukbb_pattern_difference_sub.sh haplotype_estimates.ukbb_bca_cases.chr11.68700000-69700000.subset.txt,haplotype_estimates.ukbb_bca_controls.chr11.68700000-69700000.subset.txt Pattern_combined.Iteration000.population_0.chr11.68700000-69700000_2,j.txt 50 Iteration000 0.001"

# bsub -P SJLIFE -J pattern_differences_sub.Iteration003.ukbb_bca_cases+ukbb_bca_controls.chr11.1802846-1902305 -oo /research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/biobank/ChromosomeOverlap/UKBB_chr11.0/pattern_differences_sub.Iteration003.ukbb_bca_cases+ukbb_bca_controls.chr11.802846-1902305.out -eo /research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/biobank/ChromosomeOverlap/UKBB_chr11.0/pattern_differences_sub.Iteration003.ukbb_bca_cases+ukbb_bca_controls.chr11.802846-1902305.err -R "rusage[mem=128]" "/research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/biobank/ChromosomeOverlap/ukbb_pattern_difference_sub.sh haplotype_estimates.ukbb_bca_cases.chr11.1802846-1902305.txt,haplotype_estimates.ukbb_bca_controls.chr11.1802846-1902305.txt Pattern_combined.Iteration003.ukbb_bca_cases.chr11.1802846-1902305_2,j.txt 50 Iteration003 0.0001"

POPULATION=$1 # comma-separated list of two haplotypes files (subjectx x rsid with header and row names)
PATTERNS=$2 # two-column file of pattern count and patterns in the form column1_value1,column2_value2 in the first population
DELTA=$3 # optional number of patterns to evaluate at once
NAME=$4 # optional name to prepend to job subdirectories
ALPHA=$5 # optional p value cutoff for keeping patterns
REMOVE_UBIQUITOUS=$6 # whether to remove ubiquitous alleles from patterns after final list is formed
DIRECTORY=$7 #
HOME_DIR=$8

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

POPULATION=($(echo $POPULATION | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))

if [ -z $DELTA ] #number of haplotypes to evaluate per job
then
  DELTA=50
fi

MAX_JOBS=$(cat "/hpcf/lsf/lsf_prod/conf/lsbatch/hpcf_research_cluster/configdir/lsb.users" | grep "#wletsou" | awk '{print $2}')
if [ -z $MAX_JOBS ]
then
  MAX_JOBS=100 #revised limit
fi

n_patterns=$(cat $PATTERNS | wc -l)
n_jobs=$(awk 'BEGIN{printf "%0.25f\n",'$n_patterns'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
echo Number of jobs for $n_patterns pattern$( (($n_patterns>1)) && echo "s" || echo "" ) is $n_jobs.
echo Number of jobs is $n_jobs with DELTA = $DELTA pattern$( (($DELTA>1)) && echo "s" || echo "" ) per job.
while (( $n_jobs > $MAX_JOBS ))
do
  DELTA=$((2*DELTA))
  n_jobs=$(awk 'BEGIN{printf "%0.25f\n",'$n_patterns'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
  echo Revised number of jobs is $n_jobs with DELTA = $DELTA pattern$( (($DELTA>1)) && echo "s" || echo "" ) per job.
done
printf "\n"

if [ ! -z $NAME ]
then
  NAME=${NAME}"."
fi

for subdir in ${NAME}pattern_difference_* # remove previous directories for pattern segregation
do
  test -d $subdir && echo rm -r $subdir && printf "\n"
  test -d $subdir && rm -r $subdir
done

for ((i=1;i<=$n_jobs;i++)) # create subdirectories for evaluating haplotypes in groups of DELTA
do
  test -d ${NAME}pattern_difference_${i} && echo rm -r ${NAME}pattern_difference_${i}
  test -d ${NAME}pattern_difference_${i} && rm -r ${NAME}pattern_difference_${i}
  echo mkdir ${NAME}pattern_difference_${i}
  mkdir ${NAME}pattern_difference_${i}

  test -f ${POPULATION[0]} && echo cp ${POPULATION[0]} ${NAME}pattern_difference_${i}/${POPULATION[0]}
  test -f ${POPULATION[0]} && cp ${POPULATION[0]} ${NAME}pattern_difference_${i}/${POPULATION[0]}
  test -f ${POPULATION[1]} && echo cp ${POPULATION[1]} ${NAME}pattern_difference_${i}/${POPULATION[1]}
  test -f ${POPULATION[1]} && cp ${POPULATION[1]} ${NAME}pattern_difference_${i}/${POPULATION[1]}

  test -f $PATTERNS && echo cp $PATTERNS ${NAME}pattern_difference_${i}/${PATTERNS} && printf "\n"
  test -f $PATTERNS && cp $PATTERNS ${NAME}pattern_difference_${i}/${PATTERNS}

  printf "\n"
done

echo Submit haplotype jobs:
echo bsub \-P SJLIFE \-J \"myJob[1-$n_jobs]\" \-oo ${DIRECTORY}/${NAME}pattern_difference_%I/ukbb_pattern_difference.%I.out \-eo ${DIRECTORY}/${NAME}pattern_difference_%I/ukbb_pattern_difference.%I.err \-R \"rusage[mem=2048]\" \-R \"select[ut \< 0.8]\" \-R \"order[!ut]\" \"sh ${HOME_DIR}/ukbb_pattern_difference3.sh ${POPULATION[0]},${POPULATION[1]} $PATTERNS ${DELTA}.\\\$LSB_JOBINDEX $DIRECTORY/${NAME}pattern_difference_\\\$LSB_JOBINDEX $HOME_DIR\"
job_id=$(bsub -P SJLIFE -J "myJob[1-$n_jobs]" -oo ${DIRECTORY}/${NAME}pattern_difference_%I/ukbb_pattern_difference.%I.out -eo ${DIRECTORY}/${NAME}pattern_difference_%I/ukbb_pattern_difference.%I.err -R "rusage[mem=2048]" -R "select[ut < 0.8]" -R "order[!ut]" "sh ${HOME_DIR}/ukbb_pattern_difference3.sh ${POPULATION[0]},${POPULATION[1]} $PATTERNS ${DELTA}.\$LSB_JOBINDEX $DIRECTORY/${NAME}pattern_difference_\$LSB_JOBINDEX $HOME_DIR")
printf "\n"
echo $job_id
job_id=$(echo $job_id | awk 'b=gensub(/.*<([0-9]*)>.*/,"\\1","g",$0) {print b}') #extract job_id (number) from output

echo Wait until all jobs start:
echo bsub \-P SJLIFE \-J sleep.ukbb_pattern_difference \-w \"numrun\($job_id,\*\) \|\| numended\($job_id,\*\)\" \-R \"rusage[mem=32]\" \-oo ${DIRECTORY}/sleep.ukbb_pattern_difference.out \-eo ${DIRECTORY}/sleep.ukbb_pattern_difference.err \-K \"sleep 10\"
bsub -P SJLIFE -J sleep.ukbb_pattern_difference -w "numrun($job_id,*) || numended($job_id,*)" -R "rusage[mem=32]" -oo ${DIRECTORY}/sleep.ukbb_pattern_difference.out -eo ${DIRECTORY}/sleep.ukbb_pattern_difference.err -K "sleep 10"
printf "\n"

echo Wait until all jobs done:
job_array=($(bjobs -J myJob 2> /dev/null | awk 'NR>1{print $7}')) # get job name from 7th field in all non-header rows of bjobs
declare -p job_array
while (( ${#job_array[@]}>0 ))
do
  sleep 60
  job_array=($(bjobs -J myJob 2> /dev/null | awk 'NR>1{print $7}'))
done
printf "\n"

filename=${PATTERNS/combined/differences}
filename=${filename/.txt/.p_vals.txt} # p values from t-test or Fisher's exact test
test -f $filename && echo rm $filename && printf "\n"
test -f $filename && rm $filename
echo touch $filename && printf "\n"
touch $filename

for subdir in ${NAME}pattern_difference_*
do
  test -d $subdir && echo cd $subdir && printf "\n"
  test -d $subdir && cd $subdir
  if [ "${PWD##*/}" == "$subdir" ] # make sure PWD is subdir before deleting files
  then
    for file in fisher_exact*.txt
    do
      test -f $file && found+=1
      test -f $file && echo cat $file \>\> ${DIRECTORY}/$filename
      test -f $file && cat $file >> ${DIRECTORY}/$filename
    done
    test ! -z $found && (( $found>0 )) && found=0 && printf "\n"
    for file in haplotype_estimates*
    do
      test -f $file && found+=1
      test -f $file && echo rm $file
      test -f $file && rm $file
    done
    test ! -z $found && (( $found>0 )) && found=0 && printf "\n"
    test -f $PATTERNS && echo rm $PATTERNS && printf "\n"
    test -f $PATTERNS && rm $PATTERNS
  fi
  test -d $DIRECTORY && echo cd $DIRECTORY && printf "\n"
  test -d $DIRECTORY && cd $DIRECTORY
done

test -z $ALPHA && ALPHA=0.01 # p_value cutoff

# echo Get patterns with alpha \< $ALPHA:
# extract p-values from t-test results
# echo awk \'\(\$2\<\'$ALPHA\'\){print \$1}\' $filename \> ${PATTERNS%combined*}combined_alpha${ALPHA}${PATTERNS##*combined}
# awk '($2<'$ALPHA'){print $1}' $filename > ${PATTERNS%combined*}combined_alpha${ALPHA}${PATTERNS##*combined}
# printf "\n"

# extract p-values from Fisher's exact test results
# echo awk \'\(\(\$5+0\)\<\'$ALPHA\'\){print \$1}\' $filename \> ${PATTERNS%combined*}combined_alpha${ALPHA}${PATTERNS##*combined}
# awk '(($5+0)<'$ALPHA'){print $1}' $filename > ${PATTERNS%combined*}combined_alpha${ALPHA}${PATTERNS##*combined}
# printf "\n"
#
# echo Format filtered patterns in new Pattern_combined file as \"count,pattern\":
# echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{count[\$2]=\$1\; next} \(\$1 in count\){print count[\$1],\$1}\' $PATTERNS ${PATTERNS%combined*}combined_alpha${ALPHA}${PATTERNS##*combined} \> ${PATTERNS%*.txt}.tmp
# awk 'BEGIN{OFS="\t"} NR==FNR{count[$2]=$1; next} ($1 in count){print count[$1],$1}' $PATTERNS ${PATTERNS%combined*}combined_alpha${ALPHA}${PATTERNS##*combined} > ${PATTERNS%*.txt}.tmp
# printf "\n"

echo Filtered patterns by Fisher\'s exact test results and format in new Pattern_combined file as \"count,pattern\":
echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{count[\$2]=\$1\; next} \(\$1 in count\){print count[\$1],\$1}\' $PATTERNS \<\(awk \'\(\$5+0\<$ALPHA\){print \$1}\' $filename\) \> ${PATTERNS%*.txt}.tmp
awk 'BEGIN{OFS="\t"} NR==FNR{count[$2]=$1; next} ($1 in count){print count[$1],$1}' $PATTERNS <(awk '($5+0<'$ALPHA'){print $1}' $filename) > ${PATTERNS%*.txt}.tmp
printf "\n"

test -f ${PATTERNS} && echo mv ${PATTERNS} ${PATTERNS%combined*}combined_old${PATTERNS##*combined} && printf "\n"
test -f ${PATTERNS} && mv ${PATTERNS} ${PATTERNS%combined*}combined_old${PATTERNS##*combined}

test -f ${PATTERNS%*.txt}.tmp && echo mv ${PATTERNS%*.txt}.tmp $PATTERNS && printf "\n"
test -f ${PATTERNS%*.txt}.tmp && mv ${PATTERNS%*.txt}.tmp $PATTERNS

if [ -z $REMOVE_UBIQUITOUS ]
then
  REMOVE_UBIQUITOUS=0 # default is to not remove
fi

if [[ $REMOVE_UBIQUITOUS == "1" ]]
then
  echo Remove ubiquitous alleles:
  awk '{n=split($2,pattern,","); for (i=1;i<=n;i++) {array[pattern[i]]+=1} } END{for (i in array) {if (array[i]==NR) {print i} } }' $PATTERNS # find alleles appearing in all patterns
  printf "\n"
  awk 'BEGIN{OFS="\t"} NR==FNR{allele[$1]; next} NR!=FNR{b=$2; for (i in allele) {b=gensub("[,]*"i,"","g",b) }; $2=gensub("^,","","g",b); print $0 }' <(awk '{n=split($2,pattern,","); for (i=1;i<=n;i++) {array[pattern[i]]+=1} } END{for (i in array) {if (array[i]==NR) {print i} } }' $PATTERNS) $PATTERNS > ${PATTERNS%.txt}.tmp # remove ubiquitous alleles and squeeze out remaining commas
  test -f ${PATTERNS%.txt}.tmp && echo mv ${PATTERNS%.txt}.tmp $PATTERNS && printf "\n"
  test -f ${PATTERNS%.txt}.tmp && mv ${PATTERNS%.txt}.tmp $PATTERNS
fi

# echo Get pattern frequencies and counts:
# echo awk \'BEGIN{OFS=\"\\t\"} \(\$5+0\<$ALPHA\){seen[\$2]+=1} END{for \(i in seen\) {print i,seen[i]} }\' $filename \| sort \-gr \> ${PATTERNS/combined/combined_frequencies}
# awk 'BEGIN{OFS="\t"} ($5+0<'$ALPHA'){seen[$2]+=1} END{for (i in seen) {print i,seen[i]} }' $filename | sort -gr > ${PATTERNS/combined/combined_frequencies}
# printf "\n"
#
# signif_patterns=($(awk '{print $1}' ${PATTERNS/combined/combined_frequencies}))
# n_signif_patterns=${#signif_patterns[@]}
# char=$(( $(echo $n_signif_patterns | wc -c)-1 ))
#
# for ((i=0;i<${#signif_patterns[@]};i++))
# do
#   echo Find patterns with frequency ${signif_patterns[i]}:
#   echo awk \'\(\$2==${signif_patterns[i]} \&\& \$5+0\<$ALPHA\){print \$1}\' $filename \> ${PATTERNS%combined*}combined.temp${i}${PATTERNS##*combined}
#   awk '($2=='${signif_patterns[i]}' && $5+0<'$ALPHA'){print $1}' $filename > ${PATTERNS%combined*}combined.temp${i}${PATTERNS##*combined}
#   printf "\n"
#
#   echo Format patterns in new Pattern_combined file as \"count,pattern\":
#   echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{count[\$2]=\$1\; next} \(\$1 in count\){print count[\$1],\$1}\' $PATTERNS ${PATTERNS%combined*}combined.temp${i}${PATTERNS##*combined} \> ${PATTERNS%combined*}combined.${i}${PATTERNS##*combined}
#   awk 'BEGIN{OFS="\t"} NR==FNR{count[$2]=$1; next} ($1 in count){print count[$1],$1}' $PATTERNS ${PATTERNS%combined*}combined.temp${i}${PATTERNS##*combined} > ${PATTERNS%combined*}combined.$(printf "%0.${char}d" ${i})${PATTERNS##*combined}
#   printf "\n"
#
#   test -f ${PATTERNS%combined*}combined.temp${i}${PATTERNS##*combined} && echo rm ${PATTERNS%combined*}combined.temp${i}${PATTERNS##*combined} && printf "\n"
#   test -f ${PATTERNS%combined*}combined.temp${i}${PATTERNS##*combined} && rm ${PATTERNS%combined*}combined.temp${i}${PATTERNS##*combined}
# done
