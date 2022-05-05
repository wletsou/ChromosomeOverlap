#! /bin/bash
set -e

# bsub -P SJLIFE -J ukbb_haplotype_partition_sub -oo ukbb_haplotype_partition_sub.out -eo ukbb_haplotype_partition_sub.err -R "rusage[mem=256]" sh /home/wletsou/scripts/ukbb_haplotype_partition_sub.sh haplotype_estimates.ukbb_bca_20200116_cases.chr17.5222961-5436196.txt,haplotype_estimates.ukbb_bca_20200116_controls.chr17.5222961-5436196.txt Core_patterns_names.ukbb_bca_20200116_cases.subset_9346.chr17.5222961-5436196_2,3,j.txt 17 5222961,5436196 50

# bsub -P SJLIFE -J ukbb_haplotype_partition_sub -oo ukbb_haplotype_partition_sub.out -eo ukbb_haplotype_partition_sub.err -R "rusage[mem=256]" sh /home/wletsou/scripts/ukbb_haplotype_partition_sub.sh haplotype_estimates.ukbb_bca_20200116_cases.chr11.98807358-98897224.txt,haplotype_estimates.ukbb_bca_20200116_controls.chr11.98807358-98897224.txt Core_patterns_names.ukbb_bca_20200116_cases.subset_9346.chr11.98807358-98897224_2,3,j.txt 11 98807358,98897224 50

POPULATION=$1 # comma-separated list of case,control haplotype_estimates files
HAPLOTYPES=$2 # haplotype frequency file with iteration in column 1, length or count in column 2, name in column 3
CHR=$3 # chromosome number
BP_RANGE=$4 # comma-separated list from-bp,to-bp (separate with semicolons or commas for multiple ranges; use quotes)
DELTA=$5 # number of haplotype to evaluate in a single job
ITERATION=$6 # optional value of overlap iteration for selecting patterns to partition
DIRECTORY=$7 # location to store output/where .indiv files are found
HOME_DIR=$8 # location of program files

module load R/3.6.1

printf "\n"
if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
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

if [ -f $HAPLOTYPES ]
then
  # if HAPLOTYPES supplied as the third column of a file
  if [ ! -z $ITERATION ]
  then
    HAPLOTYPES_LIST=($(awk '($1 == '$ITERATION' ){printf "%s ",$3}' $HAPLOTYPES)) # select only haplotypes found at iteration ITERATION
  else
    HAPLOTYPES_LIST=($(awk '{printf "%s ",$3}' $HAPLOTYPES)) # select haplotypes at all iterations
  fi
else
  # if HAPLOTYPES supplied as a colon- (:) separated list
  HAPLOTYPES_LIST=($(echo $HAPLOTYPES_LIST | perl -pne 's/([0-9A-Za-z_.,=-]+)[:;]*/$1 /g')) # separate snp sets by (:) or (;), but not (,)
fi

if [ ! -z $ITERATION ]
then
  Iteration=$(printf ".Iteration%03d" $ITERATION)
fi

n_haplotypes=${#HAPLOTYPES_LIST[@]}
n_jobs=$(awk 'BEGIN{printf "%0.25f\n",'$n_haplotypes'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
echo Number of jobs for $n_haplotypes haplotype$( (($n_haplotypes>1)) && echo "s" || echo "" ) is $n_jobs.
echo Number of jobs is $n_jobs with DELTA = $DELTA haplotype$( (($DELTA>1)) && echo "s" || echo "" ) per job.
while (( $n_jobs > $MAX_JOBS ))
do
  DELTA=$((2*DELTA))
  n_jobs=$(awk 'BEGIN{printf "%0.25f\n",'$n_haplotypes'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
  echo Revised number of jobs is $n_jobs with DELTA = $DELTA haplotype$( (($DELTA>1)) && echo "s" || echo "" ) per job.
done
printf "\n"

for subdir in fisher_exact${Iteration}_* # remove previous directories for pattern segregation
do
  test -d $subdir && echo rm -r $subdir && printf "\n"
  test -d $subdir && rm -r $subdir && printf "\n"
done

for ((i=1;i<=$n_jobs;i++)) # create subdirectories for evaluating haplotypes in groups of DELTA
do
  test -d fisher_exact${Iteration}_${i} && echo rm -r fisher_exact${Iteration}_${i}
  test -d fisher_exact${Iteration}_${i} && rm -r fisher_exact${Iteration}_${i}
  echo mkdir fisher_exact${Iteration}_${i}
  mkdir fisher_exact${Iteration}_${i}

  test -f ${POPULATION[0]} && cp ${POPULATION[0]} fisher_exact${Iteration}_${i}/${POPULATION[0]}
  test -f ${POPULATION[0]} && echo cp ${POPULATION[0]} fisher_exact${Iteration}_${i}/${POPULATION[0]}
  test -f ${POPULATION[1]} && cp ${POPULATION[1]} fisher_exact${Iteration}_${i}/${POPULATION[1]}
  test -f ${POPULATION[1]} && echo cp ${POPULATION[1]} fisher_exact${Iteration}_${i}/${POPULATION[1]}

  test -f $HAPLOTYPES && echo cp $HAPLOTYPES fisher_exact${Iteration}_${i}/${HAPLOTYPES}
  test -f $HAPLOTYPES && cp $HAPLOTYPES fisher_exact${Iteration}_${i}/${HAPLOTYPES}

  printf "\n"
done

echo Submit haplotype jobs:
echo bsub \-P SJLIFE \-J \"myJob[1-$n_jobs]\" \-oo ${DIRECTORY}/fisher_exact${Iteration}_%I/ukbb_haplotype_partition.%I.out \-eo ${DIRECTORY}/fisher_exact${Iteration}_%I/ukbb_haplotype_partition.%I.err \-R \"rusage[mem=256]\" \-R \"select[ut \< 0.9]\" \-R \"order[!ut]\" \"sh ${HOME_DIR}/ukbb_haplotype_partition.sh ${POPULATION[0]},${POPULATION[1]} $HAPLOTYPES $CHR ${BP_RANGE[0]},${BP_RANGE[1]} ${DELTA}.\\\$LSB_JOBINDEX "$ITERATION" $DIRECTORY/fisher_exact${Iteration}_\\\$LSB_JOBINDEX $HOME_DIR\"
job_id=$(bsub -P SJLIFE -J "myJob[1-$n_jobs]" -oo ${DIRECTORY}/fisher_exact${Iteration}_%I/ukbb_haplotype_partition.%I.out -eo ${DIRECTORY}/fisher_exact${Iteration}_%I/ukbb_haplotype_partition.%I.err -R "rusage[mem=256]" -R "select[ut < 0.9]" -R "order[!ut]" "sh ${HOME_DIR}/ukbb_haplotype_partition.sh ${POPULATION[0]},${POPULATION[1]} $HAPLOTYPES $CHR ${BP_RANGE[0]},${BP_RANGE[1]} ${DELTA}.\$LSB_JOBINDEX "$ITERATION" $DIRECTORY/fisher_exact${Iteration}_\$LSB_JOBINDEX $HOME_DIR")
printf "\n"

job_id=$(echo $job_id | awk 'b=gensub(/.*<([0-9]*)>.*/,"\\1","g",$0) {print b}') #extract job_id (number) from output

echo Wait until all haplotypes evaluated:
echo bsub \-P SJLIFE \-J sleep.ukbb_haplotype_partition \-w \"done\($job_id\)\" \-R \"rusage[mem=32]\" \-oo ${DIRECTORY}/sleep.ukbb_haplotype_partition.out \-eo ${DIRECTORY}/sleep.ukbb_haplotype_partition.err \-K \"sleep 10\"
bsub -P SJLIFE -J sleep.ukbb_haplotype_partition -w "done($job_id)" -R "rusage[mem=32]" -oo ${DIRECTORY}/sleep.ukbb_haplotype_partition.out -eo ${DIRECTORY}/sleep.ukbb_haplotype_partition.err -K "sleep 10"
printf "\n"

echo Collect Fisher exact p values:
ll=0 # initial lower limit of patterns
ul=0 # initial upper limit of patterns
touch fisher_exact.tmp
for subdir in fisher_exact${Iteration}_*
do
  cd $subdir
  if [ "${PWD##*/}" == "$subdir" ] # make sure PWD is subdir before deleting files
  then
    for file in fisher_exact.pattern*
    do
      if [ -f $file ]
      then
        found+=1
        ll_temp=$(echo $file | awk '{b=gensub(".*patterns_([0-9]*)-([0-9]*).*","\\1","g",$0); print b}')
        ul_temp=$(echo $file | awk '{b=gensub(".*patterns_([0-9]*)-([0-9]*).*","\\2","g",$0); print b}')
        (($ll_temp<$ll)) && ll=$ll_temp
        (($ul_temp>$ul)) && ul=$ul_temp
        echo cat $file \>\> ${DIRECTORY}/fisher_exact.tmp
        cat $file >> ${DIRECTORY}/fisher_exact.tmp
      fi
    done
    test ! -z $found && (( $found>0 )) && found=0 && printf "\n"
    for file in haplotype_estimates*
    do
      test -f $file && found+=1
      test -f $file && echo rm $file
      test -f $file && rm $file
    done
    test ! -z $found && (( $found>0 )) && found=0 && printf "\n"
    test -f $HAPLOTYPES && echo rm $HAPLOTYPES && printf "\n"
    test -f $HAPLOTYPES && rm $HAPLOTYPES
  fi
  cd $DIRECTORY
done
test -f fisher_exact.tmp && echo mv fisher_exact.tmp fisher_exact.patterns_${ll}-${ul}.txt
test -f fisher_exact.tmp && mv fisher_exact.tmp fisher_exact.patterns_${ll}-${ul}.txt
