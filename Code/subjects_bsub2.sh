#! /bin/bash
set -e

POPULATION=$1 # comma-separated list of haplotype_estimates files
CHR=$2 # chromosome where snps are located
BP_RANGE=$3 # comma-separated list of lower and upper bounds on chromosome CHR where snps come from
SIGMA=$4 # number of chromosomes to join next to each chromosome
STEP_SIZE=$5 #how many jobs to submit in one array
DELTA=$6 # number of combinations to do in one job
NAMES=$7 # optional comma-separated list of population names in output file
DIRECTORY=$8 # home directory for storing output
HOME_DIR=$9 # location of program files

POPULATION=($(echo $POPULATION | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))

MAX_JOBS=$(cat "/hpcf/lsf/lsf_prod/conf/lsbatch/hpcf_research_cluster/configdir/lsb.users" | grep "#wletsou" | awk '{print $2}')
if [ -z $MAX_JOBS ]
then
  MAX_JOBS=100 # revised limit
fi

BP_RANGE=($(echo $BP_RANGE | perl -pne 's/([0-9]+)[,]*/$1 /g'))

if [ ! -z $SIGMA ]; # if non-empty, use specified overlap at each step
then
  SIGMA=$SIGMA
else
  SIGMA=2 # default overlap
fi

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi

n_overlaps=()
for ((i=0;i<${#POPULATION[@]};i++))
do
  if [ -f ${POPULATION[i]} ]
  then
    #get number of overlaps from binomial coefficient (num0 choose SIMGA)
    num0=$(awk 'NR>1{ x+=1 } END{ print x }' ${POPULATION[i]}) # total number of haplotypes excluding header row
    if (($num0>=$SIGMA))
    then
      let num=1;
      let den=1
      for ((j=1; j<=$SIGMA; j++));
      do
        let num=$num*$((num0-j+1)) # error if negative
        let den=$den*$j
      done
      let n_tuples=($num/$den)
    else
      n_tuples=0
    fi
    n_overlaps+=($n_tuples) # add to array of counts https://stackoverflow.com/questions/1951506/add-a-new-element-to-an-array-without-specifying-the-index-in-bash
  fi
done

N_OVERLAPS=0
for n in ${n_overlaps[@]}
do
  (( $n > $N_OVERLAPS )) && N_OVERLAPS=$n
done
echo Maximum number of overlaps to be performed is $N_OVERLAPS.
# N_JOBS=$(echo $DELTA | awk -v x=$N_OVERLAPS '{print x/$1}') #http://evc-cit.info/cit052/pass_to_awk.html
# N_JOBS=$(echo $N_JOBS | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }') #https://unix.stackexchange.com/questions/131073/awk-printf-number-in-width-and-round-it-up
N_JOBS=$(awk 'BEGIN{printf "%0.25f\n",'$N_OVERLAPS'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
echo Number of jobs for $N_OVERLAPS overlap$( (($N_OVERLAPS>1)) && echo "s" || echo "" ) is $N_JOBS with DELTA = $DELTA overlap$( (($DELTA>1)) && echo "s" || echo "" ) per job.
while (( $N_JOBS > $MAX_JOBS ))
do
  DELTA=$((2*DELTA))
  N_JOBS=$(awk 'BEGIN{printf "%0.25f\n",'$N_OVERLAPS'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
  echo Revised number of jobs is $N_JOBS with DELTA = $DELTA overlap$( (($DELTA>1)) && echo "s" || echo "" ) per job.
done
printf "\n"
for ((i=1; i<=$((($((N_JOBS-1))/STEP_SIZE)+1)); i++))
do

  # Create subdirectory with necessary files
  char=$(echo $N_JOBS | wc -m) # one more than the number of character in the total number of jobs
  char=$(($char-1))
  SUBDIR=$(printf "Iteration000.Step%0.${char}d" $i)

  test -d $SUBDIR && echo rm -r $SUBDIR
  test -d $SUBDIR && rm -r $SUBDIR
  test ! -d $SUBDIR && echo mkdir $SUBDIR && printf "\n"
  test ! -d $SUBDIR && mkdir $SUBDIR

  test -d $SUBDIR && echo cd $SUBDIR && printf "\n"
  test -d $SUBDIR && cd $SUBDIR

  ll=$(($((i-1))*STEP_SIZE+1)); # jobs lower limit
  ul=$((i*STEP_SIZE)); # jobs upper limit
  if (($N_JOBS < $ul));
  then
    ul=$N_JOBS
  fi

  # copy files to SUBDIR, one haplotypes_estimates_transpose file for each job.
  echo Copy haplotype files to ${SUBDIR}:
  for ((j=$ll; j<=$ul; j++));
  do
    files=""
    for ((i=0;i<${#POPULATION[@]};i++))
    do
      filename=${POPULATION[i]/estimates/estimates_transpose}
      test -f ${DIRECTORY}/${filename} || (>&2 echo "transposed haplotype estimates file ${filename} does not exist in ${DIRECTORY}"; exit 1)
      test -f ${DIRECTORY}/${filename} && echo cp ${DIRECTORY}/${filename} ${filename%%.*}.job${j}.${filename#*.} && printf "\n"
      test -f ${DIRECTORY}/${filename} && cp ${DIRECTORY}/${filename} ${filename%%.*}.job${j}.${filename#*.} # rename as haplotype_estimates_transpose.job_${j}....

      files=${files}$(test ! -z $files && echo "," || echo "")${filename%%.*}".job\\\$LSB_JOBINDEX."${filename#*.} # string of haplotype_estimates_transpose files with placeholder for job index; escape slashes so $LSB_JOBINDEX is not expanded by shell in eval command
    done

  done

  echo Perform overlap jobs ${ll}-${ul}:
  str=$(echo "bsub -P SJLIFE -J \"myJob[$ll-$ul]\" -eo $DIRECTORY/${SUBDIR}/overlaps.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.${DELTA}.%I.err -oo $DIRECTORY/${SUBDIR}/overlaps.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.${DELTA}.%I.out -R \"rusage[mem=256]\" -n 1 -R \"span[hosts=1]\" -R \"select[ut < 0.8]\" -R \"order[!ut]\" \"sh ${HOME_DIR}/subjects_overlap_sub2.sh $files $CHR ${BP_RANGE[0]},${BP_RANGE[1]} $SIGMA ${DELTA}.\\\$LSB_JOBINDEX \\\"$NAMES\\\" $DIRECTORY $HOME_DIR\"") # haplotype_estimates_transpose files in DIRECTORY/SUBDIR and job submitted DIRECTORY/SUBDIR in but output Pattern files go to DIRECTORY
  printf "\n"
  echo $str
  eval $str
  printf "\n"

  test -d $DIRECTORY && echo cd $DIRECTORY && printf "\n"
  test -d $DIRECTORY && cd $DIRECTORY
done

# while jobs are running, requeue those that are suspended
echo Checking for suspended jobs...
unset job_array
job_array=($(bjobs -J myJob 2> /dev/null | awk '{print $7}')) # array of "myJob" if myJob is being executed
n_running_jobs=${#job_array[@]}
while (( $n_running_jobs>0 ))
do
  requeue_array=($(bjobs -J myJob 2> /dev/null | awk '$3 ~ /SSUSP/{print $7}')) # array of "myJob" that are in the SSUSP state
  if (( ${#requeue_array[@]} > 0 ))
  then
    for ((j=0;j<${#requeue_array[@]};j++))
    do
      echo Requeuing ${requeue_array[j]}:
      echo brequeue -J ${requeue_array[j]}
      brequeue -J ${requeue_array[j]}
      printf "\n"
      while [[ "$(bjobs -J ${requeue_array[j]} | awk 'NR>1{if ($3 ~ /SSUSP/) {print 1} }')" == "1" ]] # wait while job is requeued
      do
        if [[ ! "$(bjobs -J ${requeue_array[j]} | awk 'NR>1{if ($3 ~ /SSUSP/) {print 1} }')" == "1" ]] # continue when successful
        then
          continue
        fi
      done
    done
  fi
  unset job_array
  unset requeue_array
  job_array=($(bjobs -J myJob 2> /dev/null | awk '{print $7}'))
  n_running_jobs=${#job_array[@]}
done

# for ((i=1; i<=$(($((N_JOBS/STEP_SIZE))+1)); i++))
# do
#   ll=$(($((i-1))*STEP_SIZE+1));
#   char=$(echo $N_JOBS | wc -m) #one more than the number of character in the total number of jobs
#   char=$(($char-1))
#   SUBDIR=$(printf "%s/Iteration000.Step%0.${char}d" ${DIRECTORY} $i)
#   #get JOBID of array starting with index $ll
#
#   array1=($(bjobs -J myJob[$ll] | grep ^JOBID)) #first line of headers
#   array2=($(bjobs -J myJob[$ll] | grep ^[0-9])) #second line of values
#   if ((${#array1} > 0)) #test if JOBID is the first header
#   then
#     JOBID=${array2[0]}
#     str=$(echo "bsub -P SJLIFE -R \"rusage[mem=32]\" -eo ${SUBDIR}/sleep.err -oo ${SUBDIR}/sleep.out -g /done_overlaps${DIRECTORY////__} -w \"done(${JOBID})\" -K \"sleep 10\"")
#     printf "\n"
#     echo $str
#     eval $str
#     printf "\n"
#   fi
#
#   ll=$(($((i-1))*STEP_SIZE+1)); #jobs lower limit
#   ul=$((i*STEP_SIZE)); #jobs upper limit
#   if (($N_JOBS < $ul));
#   then
#     ul=$N_JOBS
#   fi
# done
