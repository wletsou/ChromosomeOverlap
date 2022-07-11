#! /bin/bash
set -e

# Submits groups of $DELTA haplotypes in the $HAPLOTYPES file to ChromosomeOverlap_haplotype_count.sh and ChromosomeOverlap_fisher_exact.R to evaluate p-values for thesegregation between cases and controls

# sh HOME_DIR/ChromosomeOverlap_haplotype_count_sub.v2.sh cases_haplotypes,controls_haplotypes Pattern_combined.Iteration000.NAME.2,j.txt 50 "" "Iteration000.NAME" DIRECTORY HOME_DIR

POPULATION=$1 # comma-separated list of case,control haplotype_estimates files
HAPLOTYPES=$2 # haplotype frequency with name in last column and optional iteration in first
DELTA=$3 # number of haplotype to evaluate in a single job
ITERATION=$4 # optional value of overlap iteration (comma-separated list) for selecting patterns to partition; only use if HAPLOTYPES's first column is iteration number
NAME=$5 # optional name for output files
DIRECTORY=$6 # location to store output/where .indiv files are found
HOME_DIR=$7 # location of program files

module load R/3.6.1

if [ -z $HOME_DIR ];
then
  HOME_DIR="" # in case HOME_DIR is on your PATH
else
  HOME_DIR="${HOME_DIR}/"
fi
echo Home directory is $HOME_DIR

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY
echo Directory path is $DIRECTORY
printf "\n"

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

if [ ! -z $HAPLOTYPES ] && [ -f $HAPLOTYPES ]
then
  echo awk \'\(\$1 \~ /\\y"(${ITERATION//,/|})"\\y/\){print \$NF}\' $HAPLOTYPES \> ${HAPLOTYPES%.*}.tmp
  awk '($1 ~ /\y'"(${ITERATION//,/|})"'\y/){print $NF}' $HAPLOTYPES > ${HAPLOTYPES%.*}.tmp && printf "\n"
  n_haplotypes=$(cat ${HAPLOTYPES%.*}.tmp | wc -l)
  char=$(($( echo $n_haplotypes | wc -m)-1))

  n_jobs=$(awk 'BEGIN{printf "%0.25f\n",'$n_haplotypes'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
  echo Number of jobs for $n_haplotypes haplotype$( (($n_haplotypes!=1)) && echo "s" || echo "" ) is $n_jobs.
  echo Number of jobs is $n_jobs with DELTA = $DELTA haplotype$( (($DELTA>1)) && echo "s" || echo "" ) per job.
  while (( $n_jobs > $MAX_JOBS ))
  do
    DELTA=$((2*DELTA))
    n_jobs=$(awk 'BEGIN{printf "%0.25f\n",'$n_haplotypes'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
    echo Revised number of jobs is $n_jobs with DELTA = $DELTA haplotype$( (($DELTA!=1)) && echo "s" || echo "" ) per job.
  done
  printf "\n"

  # n_char=$(($(echo $n_jobs | wc -m)-1))
  # echo split \--numeric-suffixes=1 a ${n_char} \-l $DELTA ${HAPLOTYPES%.*}.tmp ${HAPLOTYPES%.*}. \--additional-suffix=\".txt\"
  # split --numeric-suffixes=1 -a ${n_char} -l $DELTA ${HAPLOTYPES%.*}.tmp ${HAPLOTYPES%.*}. --additional-suffix=".txt" && printf "\n" # https://stackoverflow.com/questions/3194349/how-do-i-split-a-file-into-n-no-of-parts

  # write to new file if remained when dividing line number by DELTA is 1; suffixes are not fixed width
  echo awk \'NR%${DELTA}==1{i++\; file=sprintf\(\"%s.%d.txt\",\"${HAPLOTYPES%.*}\",i\)} {print \> file}\' ${HAPLOTYPES%.*}.tmp
  awk 'NR%'${DELTA}'==1{i++; file=sprintf("%s.%d.txt","'${HAPLOTYPES%.*}'",i)} {print > file}' ${HAPLOTYPES%.*}.tmp && printf "\n" # https://stackoverflow.com/questions/14973815/how-to-split-a-file-using-a-numeric-suffix

  test -f ${HAPLOTYPES%.*}.tmp && echo rm ${HAPLOTYPES%.*}.tmp
  test -f ${HAPLOTYPES%.*}.tmp && rm ${HAPLOTYPES%.*}.tmp && printf "\n"
else
  test -z $HAPLOTYPES && (>&2 echo "HAPLOTYPES file must be supplied."; exit 1)
  test -f $HAPLOTYPES || (>&2 echo "HAPLOTYPES file not found."; exit 1)
fi

echo Submit haplotype jobs:
echo bsub \-P SJLIFE \-J \"myJob[1-$n_jobs]\" \-oo ${DIRECTORY}/ChromosomeOverlap_haplotype_count.v2.%I.out \-eo ${DIRECTORY}/ChromosomeOverlap_haplotype_count.v2.%I.err \-R \"rusage[mem=1024]\" \-R \"select[ut \< 0.9]\" \-R \"order[!ut]\" \"sh ${HOME_DIR}ChromosomeOverlap_haplotype_count.v2.sh ${POPULATION[0]},${POPULATION[1]} ${HAPLOTYPES%.*}.\\\$LSB_JOBINDEX.txt \\\"${DELTA}.\\\$LSB_JOBINDEX\\\" \\\"$NAME\\\" $DIRECTORY $HOME_DIR\"
job_id=$(bsub -P SJLIFE -J "myJob[1-$n_jobs]" -oo ${DIRECTORY}/ChromosomeOverlap_haplotype_count.v2.%I.out -eo ${DIRECTORY}/ChromosomeOverlap_haplotype_count.v2.%I.err -R "rusage[mem=1024]" -R "select[ut < 0.9]" -R "order[!ut]" "sh ${HOME_DIR}ChromosomeOverlap_haplotype_count.v2.sh ${POPULATION[0]},${POPULATION[1]} ${HAPLOTYPES%.*}.\$LSB_JOBINDEX.txt \"${DELTA}.\$LSB_JOBINDEX\" \"$NAME\" $DIRECTORY $HOME_DIR")
echo $job_id
job_id=$(echo $job_id | awk 'b=gensub(/.*<([0-9]*)>.*/,"\\1","g",$0) {print b}') #extract job_id (number) from output
printf "\n"

echo Wait until all jobs start:
echo bsub \-P SJLIFE \-J sleep.ChromosomeOverlap_haplotype_count_sub.v2 \-w \"numrun\($job_id,\*\) \|\| numended\($job_id,\*\)\" \-R \"rusage[mem=32]\" \-oo ${DIRECTORY}/sleep.ChromosomeOverlap_haplotype_count_sub.v2.out \-eo ${DIRECTORY}/sleep.ChromosomeOverlap_haplotype_count_sub.v2.err \-K \"sleep 10\"
bsub -P SJLIFE -J sleep.ChromosomeOverlap_haplotype_count_sub.v2 -w "numrun($job_id,*) || numended($job_id,*)" -R "rusage[mem=32]" -oo ${DIRECTORY}/sleep.ChromosomeOverlap_haplotype_count_sub.v2.out -eo ${DIRECTORY}/sleep.ChromosomeOverlap_haplotype_count_sub.v2.err -K "sleep 10"
printf "\n"

echo Wait until all jobs done:
job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}')) # get job name from 7th field (or 6th if no exectution host yet) in all non-header rows of bjobs
declare -p job_array
while (( ${#job_array[@]}>0 ))
do
  sleep 60
  job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}'))
  # declare -p job_array
done
printf "\n"

for file in ${HAPLOTYPES%.*}.*.txt
do
  test -f $file && echo rm $file || echo ""
  test -f $file && rm $file || echo ""
done
printf "\n"

if [ ! -z $NAME ]
then
  NAME=".$NAME"
fi
test -f fisher_exact${NAME}.tmp && echo rm fisher_exact${NAME}.tmp || printf ""
test -f fisher_exact${NAME}.tmp && rm fisher_exact${NAME}.tmp || printf ""
echo touch fisher_exact${NAME}.tmp
touch fisher_exact${NAME}.tmp && printf "\n"

echo Evaluate Fisher exact p values:
i=1 # job counter
for file in haplotype_segregation${NAME}.patterns*
do
  if [ -f $file ]
  then
    ll_temp=$(echo $file | awk '{b=gensub(".*patterns_([0-9]*)-([0-9]*).*","\\1","g",$0); print b}')
    ul_temp=$(echo $file | awk '{b=gensub(".*patterns_([0-9]*)-([0-9]*).*","\\2","g",$0); print b}')

    f_name=${file##*haplotype_segregation.}
    f_name=${f_name%*.txt}

    echo bsub \-P SJLIFE \-J myJob[$i] \-oo ${DIRECTORY}/ChromosomeOverlap_fisher_exact.%I.out \-eo ${DIRECTORY}/ChromosomeOverlap_fisher_exact.%I.err \-R \"rusage[mem=1024]\" \-R \"select[ut \< 0.9]\" \-R \"order[!ut]\" \"Rscript ${HOME_DIR}ChromosomeOverlap_fisher_exact.R \-f $file \-d ${DIRECTORY} \-o ${f_name}\"
    bsub -P SJLIFE -J myJob[$i] -oo ${DIRECTORY}/ChromosomeOverlap_fisher_exact.%I.out -eo ${DIRECTORY}/ChromosomeOverlap_fisher_exact.%I.err -R "rusage[mem=1024]" -R "select[ut < 0.9]" -R "order[!ut]" "Rscript ${HOME_DIR}ChromosomeOverlap_fisher_exact.R -f $file -d ${DIRECTORY} -o ${f_name}" && printf "\n"
    found+=1

    let i=$i+1
  fi
done
test ! -z $found && unset found && printf "\n" || printf ""

echo Wait until all jobs done:
job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}')) # get job name from 7th field (or 6th if no exectution host yet) in all non-header rows of bjobs
declare -p job_array
while (( ${#job_array[@]}>0 ))
do
  sleep 60
  job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}'))
  # declare -p job_array
done
printf "\n"

echo Remove counts files:
for file in haplotype_segregation${NAME}.patterns*
do
  if [ -f $file ]
  then
    found+=1
    echo rm $file
    rm $file
  fi
done
test ! -z $found && unset found && printf "\n" || printf ""

echo Collect Fisher exact p values:
ll=0 # initial lower limit of patterns
ul=0 # initial upper limit of patterns
for file in fisher_exact${NAME}.patterns*
do
  if [ -f $file ]
  then
    found+=1
    ll_temp=$(echo $file | awk '{b=gensub(".*patterns_([0-9]*)-([0-9]*).*","\\1","g",$0); print b}')
    ul_temp=$(echo $file | awk '{b=gensub(".*patterns_([0-9]*)-([0-9]*).*","\\2","g",$0); print b}')
    (($ll_temp<$ll)) && ll=$ll_temp
    (($ul_temp>$ul)) && ul=$ul_temp
    echo cat $file \>\> ${DIRECTORY}/fisher_exact${NAME}.tmp
    cat $file >> ${DIRECTORY}/fisher_exact${NAME}.tmp
    test -f $file && echo rm $file || printf ""
    test -f $file && rm $file || printf ""
  fi
done
test ! -z $found && unset found && printf "\n" || printf ""
ll=1
test -f ${DIRECTORY}/fisher_exact${NAME}.tmp && ul=$(cat ${DIRECTORY}/fisher_exact${NAME}.tmp | wc -l)
echo awk \'BEGIN{OFS=\"\\t\"\; print \"pattern\",\"cases_freq\",\"controls_freq\",\"OR\",\"p\"}\' \> fisher_exact${NAME}.patterns_$(printf "%0.${char}d-%0.${char}d" ${ll} ${ul}).txt
awk 'BEGIN{OFS="\t"; print "pattern","cases_freq","controls_freq","OR","p"}' > fisher_exact${NAME}.patterns_$(printf "%0.${char}d-%0.${char}d" ${ll_new} ${ul_new}).txt && printf "\n" || printf ""

test -f fisher_exact${NAME}.tmp && echo cat fisher_exact${NAME}.tmp \>\> fisher_exact${NAME}.patterns_$(printf "%0.${char}d-%0.${char}d" ${ll} ${ul}).txt || printf ""
test -f fisher_exact${NAME}.tmp && cat fisher_exact${NAME}.tmp >> fisher_exact${NAME}.patterns_$(printf "%0.${char}d-%0.${char}d" ${ll} ${ul}).txt && printf "\n" || printf ""

test -f fisher_exact${NAME}.tmp && echo rm fisher_exact${NAME}.tmp || printf ""
test -f fisher_exact${NAME}.tmp && rm fisher_exact${NAME}.tmp && printf "\n" || printf ""
