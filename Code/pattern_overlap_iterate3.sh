#! /bin/bash
set -e

POPULATION=$1 # optional list of haplotype_estimates files for comparing pattern frequencies
CHR=$2 # chromosome where snps are located
BP_RANGE=$3 # comma-separated list of lower and upper bounds on chromosome CHR where snps come from
STEP_SIZE=$4 # how many jobs to submit in one array
DELTA=$5 # number of combinations to do in one job
SIGMA=$6 # number of patterns to be overlapped in one comparison job
PATTERN=$7 # specify a suffix to test a specfic combination, e.g., 2,3,j or 2+3+j; else uses all
ALPHA=$8 # optional p value cutoff at each step
NAMES=$9 # optional comma-separated list of population names in output file
DIRECTORY=${10} # location of haplotype matrices (snps x subjects); a full file path or relative ./
HOME_DIR=${11} # location of program files, subject_generate_haplotypes.sh; a full file path

module load R/3.6.1

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
echo Working directory is $DIRECTORY
printf "\n"
cd $DIRECTORY

# determine cases and controls
POPULATION=($(echo $POPULATION | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))

if [ ! -z $NAMES ]
then
  NAMES=($(echo $NAMES | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))
else
  for ((i=0;i<${#POPULATION[@]};i++))
  do
    NAMES+=("population_${i}")
  done
fi
test -z $NAMES && NAMES="population_0" # in case no NAMES or POPULATION variables provided

patterns=()
for file in Pattern_combined.Iteration0* # Pattern files from 0th iteration, each has a pattern suffix e.g., 2,3,j
do
  idx=1 # initially pattern suffix is new
  pattern_temp=${file##*_} # trim pattern suffix after last _ in filename
  for pattern_str in ${patterns[@]}
  do
    if [ ${pattern_temp%.txt} == $pattern_str ] # check if pattern already in array
    then
      idx=0; # pattern suffix is old
      break
    fi
  done
  if (($idx==1))
  then
    patterns+=(${pattern_temp%.txt}) # array of unique pattern types
  fi
done

if [ -z $PATTERN ]; # input chromosome pattern
then
  PATTERN=""
else
  found=0 # haven't seen a pattern file with the input pattern suffix yet
  for pattern_str in ${patterns[*]}
  do
    if [ "$PATTERN" == "$pattern_str" ];
    then
      found=$((found+=1))
      break
    fi
  done
  if (($found==0));
  then
    echo Pattern $PATTERN not found.  Using all.
    PATTERN=""
    printf "\n"
  else
    echo Found pattern $PATTERN.
    printf "\n"
  fi
fi

BP_RANGE=($(echo $BP_RANGE | perl -pne 's/([0-9]+)[,]*/$1 /g'))

if [ -z $SIGMA ];
then
  # SIGMA=1 # number of patterns in an overlap
  SIGMA=2 # number of patterns in an overlap
fi

MAX_JOBS=$(cat "/hpcf/lsf/lsf_prod/conf/lsbatch/hpcf_research_cluster/configdir/lsb.users" | grep "#wletsou" | awk '{print $2}')
if [ -z $MAX_JOBS ]
then
  MAX_JOBS=100 #revised limit
fi

for ((j=0;j<${#NAMES[@]};j++))
do
  ls Pattern_combined.Iteration000.${NAMES[j]}*${PATTERN}* 1> /dev/null 2> /dev/null || (>&2 echo "Pattern_combined file for ${NAMES[j]} is missing from $DIRECTORY"; exit 1)
done

iteration=0; # only 0th round of overlaps has been completed
echo Pattern_combined files found for ${PATTERN}:
for file in Pattern_combined.*${PATTERN}.txt
do
  test -f $file && echo $file && found+=1
done
test ! -z $found && (( $found>0 )) && found=0 && printf "\n"

echo Initial iteration string:
file_str=$(printf "Iteration%0.3d" $iteration)
echo $file_str
printf "\n"

n_overlaps=() #initialize array of overlaps to be performed, based on number of lines in Pattern_combined files for each population
for file in Pattern_combined.${file_str}.*${PATTERN}.txt
do
  echo Determine number of $SIGMA\-overlaps for $file
  num0=$(awk 'END{print NR}' $file)
  if (($num0>=$SIGMA))
  then
    let num=1;
    let den=1
    for ((i=1; i<=$SIGMA; i++));
    do
      let num=$num*$((num0-i+1)) # error if negative
      let den=$den*$i
    done
    let n_tuples=($num/$den)
  else
    n_tuples=0
  fi
  let n_tuples=($num/$den)
  n_overlaps+=($n_tuples) #add to array of counts https://stackoverflow.com/questions/1951506/add-a-new-element-to-an-array-without-specifying-the-index-in-bash
  n_input_files=${#n_overlaps[@]} #total number of SNP files in. Should be the same from iteration to iteration
done
printf "\n"
N_OVERLAPS=0
echo Array of tuples to be performed in each file:
declare -p n_overlaps
printf "\n"
#get maximum in array n_overlaps
for n in ${n_overlaps[@]}
do
  if (( $n > $N_OVERLAPS ));
  then
    N_OVERLAPS=$n #replace old maximum
  fi
done
echo Maximum number of overlaps to be performed is $N_OVERLAPS.
printf "\n"

N_JOBS=$(awk 'BEGIN{printf "%0.25f\n",'$N_OVERLAPS'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
echo Requested number of jobs for $N_OVERLAPS pattern$( (($N_OVERLAPS>1)) && echo "s" || echo "" ) is $N_JOBS with DELTA = $DELTA pattern$( (($DELTA>1)) && echo "s" || echo "" ) per job.
while (( $N_JOBS > $MAX_JOBS ))
do
  DELTA=$((2*DELTA))
  N_JOBS=$(awk 'BEGIN{printf "%0.25f\n",'$N_OVERLAPS'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
  echo Revised number of jobs is $N_JOBS with DELTA = $DELTA haplotype$( (($DELTA>1)) && echo "s" || echo "" ) per job.
done
printf "\n"

# Get memory rusage estimate from max of previous iteration
mem_array=($(for folder in ${file_str}*; # Loop through directories starting with Iteration000
do
  test -d $folder && cd $folder
  for file in overlaps*.out
  do
    if [ -f $file ]
    then
      mem=$(head -100 $file | grep "Max Memory " | tr -dc '0-9');
      echo $mem
    fi
  done
  cd ..
done))

max_mem=256; # used at most this much in last round
for mem in ${mem_array[@]}
do
  if (( $mem > $max_mem ))
  then
    max_mem=$mem
  fi
done
mem_req=$((2*max_mem)) # increase to new limit
# echo Memory request per job is $mem_req MB.
# printf "\n"

while (($N_OVERLAPS>0)) # continue as long as there is more than one pattern to be overlapped with another
do
  file_str=$(printf "Iteration%0.3d" $iteration) #from the last round
  iteration=$(($iteration+1)) #starts at 0, so first round creates Iteration001.Step001
  file_str_next=$(printf "Iteration%0.3d" $iteration) #in this round
  for ((i=1; i<=$((($((N_JOBS-1))/STEP_SIZE)+1)); i++)); # looping over number of submission steps; because of integer floor, only get another step if NJOBS > STEP_SIZE
  do
    # Create subdirectory with necessary files
    char=$(($(echo $N_JOBS | wc -m)-1)) # the number of characters in the total number of jobs
    SUBDIR=$(printf "%s/Iteration%0.3d.Step%0.${char}d" ${DIRECTORY} $iteration $i) # full path to SUBDIR
    ITERATION=$(printf "Iteration%0.3d" $iteration) #incremented by 1 from variable file_str
    STEP=$(printf "Step%0.${char}d" $i)

    test -d $SUBDIR && echo rm -r $SUBDIR
    test -d $SUBDIR && rm -r $SUBDIR
    test ! -d $SUBDIR && echo mkdir $SUBDIR && printf "\n"
    test ! -d $SUBDIR && mkdir $SUBDIR

    test -d $SUBDIR && echo cd $SUBDIR && printf "\n"
    test -d $SUBDIR && cd $SUBDIR

    ll=$(((i-1)*STEP_SIZE+1)); #jobs lower limit
    ul=$((i*STEP_SIZE)); #jobs upper limit
    if (($N_JOBS < $ul));
    then
      ul=$N_JOBS
    fi
    if (($ll>$ul))
    then
      continue
    fi
    echo Jobs to be performed: $ll to $ul.
    printf "\n"

    echo Moving copies of files Pattern_combined.${file_str}* to $SUBDIR
    for ((j=$ll; j<=$ul; j++));
    do
      # Move SNP_FILES for cases, controls, and combined to SUBDIR, collecting all patterns
      for file in ${DIRECTORY}/Pattern_combined.${file_str}*chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}*${PATTERN}.txt
      do
        echo cp $file ${SUBDIR}/job${j}.${file##*/}
        cp $file ${SUBDIR}/job${j}.${file##*/}
      done
      printf "\n"
    done

    mem_req=512
    if (( $mem_req>10000 )) # if requested memory is greater than 10 GB, use large_mem queue
    then
      queue="large_mem"
    else
      queue="standard"
    fi

    # cd $DIRECTORY
    # cd $SUBDIR

    # str=$(echo "bsub -P SJLIFE -J \"myJob[$ll-$ul]\" -eo ${SUBDIR}/overlaps.${file_str}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.${DELTA}.%I.err -oo ${SUBDIR}/overlaps.${file_str}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.${DELTA}.%I.out -R \"rusage[mem=$mem_req]\" -q $queue -n 1 -R \"span[hosts=1]\" -R \"select[ut < 0.8]\" -R \"order[!ut]\" \"sh ${HOME_DIR}/pattern_overlap_loop2.sh ${DELTA}.\\\$LSB_JOBINDEX $ITERATION job\\\$LSB_JOBINDEX $DIRECTORY $HOME_DIR\"")
    # echo $str
    # eval $str
    for file in ${DIRECTORY}/Pattern_combined.${file_str}*chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}*${PATTERN}.txt
    do
      str=$(echo "bsub -P SJLIFE -J \"myJob[$ll-$ul]\" -eo ${SUBDIR}/overlaps.${file##*/}.%I.err -oo ${SUBDIR}/overlaps.${file##*/}.%I.out -R \"rusage[mem=5000]\" -q large_mem -n 1 -R \"status=='ok'\" -R \"span[hosts=1]\" -R \"select[ut < 0.8]\" -R \"order[!ut]\" \"sh ${HOME_DIR}/pattern_overlap_loop4.sh ${SUBDIR}/job\\\$LSB_JOBINDEX.${file##*/} ${DELTA}.\\\$LSB_JOBINDEX 2500000 $SIGMA $ITERATION $DIRECTORY $HOME_DIR\"")
      echo $str
      eval $str
      printf "\n"
    done

    cd $DIRECTORY
  done

  cd $DIRECTORY

  echo Wait until all jobs done:
  job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}')) # get job name from 7th field (or 6th if no exectution host yet) in all non-header rows of bjobs
  declare -p job_array
  while (( ${#job_array[@]}>0 ))
  do
    sleep 20
    job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}'))
  done
  printf "\n"

  str=$(echo "bsub -P SJLIFE -J sleep.${SUBDIR##*/} -R \"rusage[mem=32]\" -R \"status=='ok'\" -eo ${SUBDIR}/sleep.err -oo ${SUBDIR}/sleep.out -K \"sleep 10\"")
  printf "\n"
  echo $str
  eval $str
  printf "\n"

  #delete extra files once all jobs complete
  # for ((i=1;i<=$((($((N_JOBS-1))/STEP_SIZE)+1));i++));
  # do
  #   ll=$(($((i-1))*STEP_SIZE+1));
  #   ul=$((i*STEP_SIZE)); #jobs upper limit
  #   if (($N_JOBS < $ul));
  #   then
  #     ul=$N_JOBS
  #   fi
  #   if (($ll>$ul))
  #   then
  #     continue
  #   fi
  #   char=$(echo $N_JOBS | wc -m) #one more than the number of character in the total number of jobs
  #   char=$(($char-1))
  #   SUBDIR=$(printf "%s/Iteration%0.3d.Step%0.${char}d" ${DIRECTORY} $iteration $i)
  #   #get JOBID of array starting with index $ll
  #
  #   array1=($(bjobs -J myJob[$ll] | grep ^JOBID)) #first line of headers
  #   array2=($(bjobs -J myJob[$ll] | grep ^[0-9])) #second line of values
  #   if ((${#array1} > 0)) #test if JOBID is the first header
  #   then
  #     JOBID=${array2[0]}
  #     # str=$(echo "bsub -P SJLIFE -J remove.${SUBDIR##*/} -g /done_overlaps${DIRECTORY////__} -R \"rusage[mem=32]\" -eo ${SUBDIR}/remove.err -oo ${SUBDIR}/remove.out -w \"done(${JOBID})\" \"${HOME_DIR}/remove.sh $SUBDIR\"")
  #     str=$(echo "bsub -P SJLIFE -J sleep.${SUBDIR##*/} -g /done_overlaps${DIRECTORY////__} -R \"rusage[mem=32]\" -eo ${SUBDIR}/sleep.err -oo ${SUBDIR}/sleep.out -w \"done(${JOBID})\" -K \"sleep 10\"")
  #     echo $str
  #     eval $str
  #     printf "\n"
  #   fi
  #
  #   ll=$(((i-1)*STEP_SIZE+1)); #jobs lower limit
  #   ul=$((i*STEP_SIZE)); #jobs upper limit
  #   if (($N_JOBS < $ul));
  #   then
  #     ul=$N_JOBS
  #   fi
  # done

  # while pattern_overlap_loop jobs are running, requeue those that are suspended
  # echo Checking for suspended jobs...
  # unset job_array
  # job_array=($(bjobs -J myJob 2> /dev/null | awk '{print $7}')) # array of "myJob" if myJob is being executed
  # n_running_jobs=${#job_array[@]}
  # while (( $n_running_jobs>0 ))
  # do
  #   requeue_array=($(bjobs -J myJob 2> /dev/null | awk '$3 ~ /SSUSP/{print $7}')) # array of "myJob" that are in the SSUSP state
  #   if (( ${#requeue_array[@]} > 0 ))
  #   then
  #     for ((j=0;j<${#requeue_array[@]};j++))
  #     do
  #       echo Requeuing ${requeue_array[j]}:
  #       echo brequeue -J ${requeue_array[j]}
  #       brequeue -J ${requeue_array[j]}
  #       printf "\n"
  #       while [[ "$(bjobs -J ${requeue_array[j]} | awk 'NR>1{if ($3 ~ /SSUSP/) {print 1} }')" == "1" ]] # wait while job is requeued
  #       do
  #         if [[ ! "$(bjobs -J ${requeue_array[j]} | awk 'NR>1{if ($3 ~ /SSUSP/) {print 1} }')" == "1" ]] # continue when successful
  #         then
  #           continue
  #         fi
  #       done
  #     done
  #   fi
  #   unset job_array
  #   unset requeue_array
  #   job_array=($(bjobs -J myJob 2> /dev/null | awk '{print $7}'))
  #   n_running_jobs=${#job_array[@]}
  # done

  #Suspend for 1 min after all jobs have completed before combining
  # str=$(echo "bsub -P SJLIFE -J sleep.${file_str_next} -w \"numpend(/done_overlaps${DIRECTORY////__},0)\" -R \"rusage[mem=32]\" -oo sleep.${file_str_next}.out -eo sleep.${file_str_next}.err -K \"sleep 10\"" )
  # echo $str
  # eval $str
  # printf "\n"
  #
  # str=$(echo "bsub -P SJLIFE -J sleep.${file_str_next} -w \"numrun(/done_overlaps${DIRECTORY////__},0)\" -R \"rusage[mem=32]\" -oo sleep.${file_str_next}.out -eo sleep.${file_str_next}.err -K \"sleep 10\"" )
  # echo $str
  # eval $str
  # printf "\n"

  str=$(echo "bsub -P SJLIFE -J sleep.${file_str_next} -R \"rusage[mem=32]\" -R \"status=='ok'\" -oo sleep.${file_str_next}.out -eo sleep.${file_str_next}.err -K \"sleep 10\"" )
  echo $str
  eval $str
  printf "\n"

  # combine output files in DIRECTORY after all of the "myJob" array have finished.  Combine by population.  Using new file_str for next iteration
  suffix_array=() # initiate array of pattern types
  for file in Overlap_tuples.${file_str_next}*.txt
  do
    echo $file
    if [ -f $file ];
    then
      suffix=${file##*_}
      pattern_str=${suffix%.*} # get pattern type after last _ and before .txt of file name
      y=1 # whether to include pattern_str in suffix array
      for str in ${suffix_array[@]}
      do
        if [ $str == $pattern_str ] # only get new pattern types
        then
          y=0
          break
        fi
      done
      if (( $y==1 ))
      then
        suffix_array+=($pattern_str)
      fi
    fi
  done
  printf "\n"

  echo Patterns found among overlaps at iteration ${iteration}:
  declare -p suffix_array
  printf "\n"

  mem_array=()
  for pattern_str in ${suffix_array[@]}
  do
    x=0
    for file in Overlap_tuples.${file_str_next}*${pattern_str}.txt
    do
      if [ -f $file ];
      then
        x=$((x+$(du -b $file | cut -f1))) # size of each Overlap_tuples file in bytes
      fi
    done
    mem_array+=($x) # store total size of all Overlap_tuples files for this pattern type
  done
  f_limit=128000000 # find the maximum total size (in b)
  for x in ${mem_array[@]}
  do
    (( $x > $f_limit )) && f_limit=$x
  done
  f_limit=$((2*f_limit/1000000)) # allocate 2 times the maximum total size (in mb) for sorting the combined files
  echo Reserving $f_limit mb for combining and sorting overlaps.

  if (( $f_limit>10000 )) # if requested memory is greater than 10 GB, use large_mem queue
  then
    queue="large_mem"
  else
    queue="standard"
  fi

  for ((j=0;j<${#NAMES[@]};j++))
  do
    # str=$(echo "bsub -P SJLIFE -J combine.${NAMES[j]}.${file_str_next} -w \"numpend(/done_overlaps${DIRECTORY////__},0)\" -q $queue -R \"rusage[mem=$f_limit]\" -oo combine.${NAMES[j]}.${file_str_next}.out -eo combine.${NAMES[j]}.${file_str_next}.err -K \"${HOME_DIR}/pattern_combine2.sh ${NAMES[j]} $CHR ${BP_RANGE[0]},${BP_RANGE[1]} $file_str_next \\\"$PATTERN\\\" $DIRECTORY $HOME_DIR\"" )
    str=$(echo "bsub -P SJLIFE -J combine.${NAMES[j]}.${file_str_next} -q $queue -R \"rusage[mem=$f_limit]\" -R \"status=='ok'\" -oo combine.${NAMES[j]}.${file_str_next}.out -eo combine.${NAMES[j]}.${file_str_next}.err -K \"${HOME_DIR}/pattern_combine2.sh ${NAMES[j]} $CHR ${BP_RANGE[0]},${BP_RANGE[1]} $file_str_next \\\"$PATTERN\\\" $DIRECTORY $HOME_DIR\"" )
    echo $str
    eval $str
    printf "\n"
  done

  if (( ${#POPULATION[@]}>1 )) # whether there is a second haplotype_estimates file with which to compare pattern frequencies
  then
    # test -z $ALPHA && ALPHA=0.000000001
    # test -z $ALPHA && ALPHA=0.0001
    if [ ! -z $ALPHA ]
    then
      echo Keep only patterns with p-value less than $ALPHA:
      echo bsub \-P SJLIFE \-J pattern_differences_sub.${file_str_next}.${NAMES[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} \-oo ${DIRECTORY}/pattern_differences_sub.${file_str_next}.${NAMES[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out \-eo ${DIRECTORY}/pattern_differences_sub.${file_str_next}.${NAMES[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err \-R \"rusage[mem=256]\" \-R \"status==\'ok\'\" \-K \"${HOME_DIR}/ukbb_pattern_difference_sub.sh ${POPULATION[0]},${POPULATION[1]} Pattern_combined.${file_str_next}.${NAMES[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt 50 ${file_str_next} $ALPHA 0 $DIRECTORY $HOME_DIR\"
      bsub -P SJLIFE -J pattern_differences_sub.${file_str_next}.${NAMES[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]} -oo ${DIRECTORY}/pattern_differences_sub.${file_str_next}.${NAMES[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.out -eo ${DIRECTORY}/pattern_differences_sub.${file_str_next}.${NAMES[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.err -R "rusage[mem=256]" -R "status=='ok'" -K "${HOME_DIR}/ukbb_pattern_difference_sub.sh ${POPULATION[0]},${POPULATION[1]} Pattern_combined.${file_str_next}.${NAMES[0]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt 50 ${file_str_next} $ALPHA 0 $DIRECTORY $HOME_DIR" # filter patterns with p value < $ALPHA and DO NOT (=0) remove ubiquitous alleles appearing in all patterns
      printf "\n"
    fi
  fi

  # get updated number of jobs
  # get maximum number of patterns in one of the Pattern_combined files
  n_overlaps=()
  echo Looking for $n_input_files file$((($n_input_files!=1)) && echo "s" || echo ""). # number corresponds to number of populations
  printf "\n"
  while ((${#n_overlaps[@]} < $n_input_files)) # wait until all input files have been combined, or count the remaining files again
  do
    for file in Pattern_combined.${file_str_next}.*${PATTERN}.txt
    do
      echo $file
      # if file exists, get number of overlaps from binomial coefficient (num0 choose SIMGA)
      if [ -f $file ]
      then
        echo Evaluating number of overlaps:
        echo awk \'END{print NR}\' $file
        echo $(awk 'END{print NR}' $file) line$( (($(awk 'END{print NR}' $file)!=1)) && echo "s" || echo "").
        num0=$(awk 'END{print NR}' $file)
        if (($num0>=$SIGMA))
        then
          let num=1;
          let den=1
          for ((i=1; i<=$SIGMA; i++));
          do
            let num=$num*$((num0-i+1)) # error if negative
            let den=$den*$i
          done
          let n_tuples=($num/$den)
        else
          n_tuples=0
        fi
        echo $n_tuples tuple$((($n_tuples!=1)) && echo "s" || echo "").
        printf "\n"
        n_overlaps+=($n_tuples) # add to array of counts https://stackoverflow.com/questions/1951506/add-a-new-element-to-an-array-without-specifying-the-index-in-bash
      fi
    done
    echo Array of tuples to be performed in each file:
    declare -p n_overlaps
    printf "\n"
    N_OVERLAPS=0
    for n in ${n_overlaps[@]}
    do
      (( $n > $N_OVERLAPS )) && N_OVERLAPS=$n
    done
  done

  # Find closed patterns not appearing in the next iteration
  echo Array of patterns:
  declare -p patterns
  printf "\n"
  for pattern_str in ${patterns[*]}
  do
    if ls Pattern_combined.${file_str_next}*chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt 1> /dev/null 2> /dev/null # only evaluate pattern_str (e.g., 2,j) if it has associated Pattern_combined files
    then
      n_bar_groups=$(($(echo ${pattern_str} | tr -c -d "+" | wc -c)+1))
      echo Number of groups is $n_bar_groups

      # Get the second column (snp pattern) of the new file (file_str_next) and print lines of the first file that do not have a match in the second column https://stackoverflow.com/questions/15251188/find-the-difference-between-two-files
      for ((j=0;j<${#NAMES[@]};j++))
      do
        echo Pattern_combined.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
        printf "\n"
        if [ -f Pattern_combined.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt ];
        then
          f_limit=$((2*$(du -b Pattern_combined.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt | cut -f1)/1000000)) #twice file size (in mb) of file of new patterns found at iteration i

          if (( $f_limit>10000 )) # if requested memory is greater than 10 GB, use large_mem queue
          then
            queue="large_mem"
          else
            queue="standard"
          fi

          N_LINES=$(awk 'END{print NR}' Pattern_combined.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt)
          echo Total number of lines in Pattern_combined.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt is $N_LINES.
          printf "\n"
          DELTA_LINES=20000

          COMPARISON_JOBS=$(echo $DELTA_LINES | awk -v x=$N_LINES '{printf "%0.25f\n", x/$1}') #exact number of comparison jobs
          COMPARISON_JOBS=$(echo $COMPARISON_JOBS | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
          if (($COMPARISON_JOBS < 1));
          then
            COMPARISON_JOBS=1
          fi
          echo Number of jobs, based on $DELTA_LINES line$( (($DELTA_LINES > 1)) && echo "s" || echo "" ) per job, is $COMPARISON_JOBS.
          printf "\n"

          # patterns in the old iteration (TEST_FILE, input $1) not in the new iteration (REF_FILE, input $2), CORE PATTERNS
          str=$(echo "bsub -P SJLIFE -J \"myJob[1-$COMPARISON_JOBS]\" -oo ${file_str}/core_pattern_comparisons.%I.${NAMES[j]}.${file_str}.out -eo ${file_str}/core_pattern_comparisons.%I.${NAMES[j]}.${file_str}.err -R \"rusage[mem=$f_limit]\" -R \"status=='ok'\" -q $queue \"sh ${HOME_DIR}/file_compare.sh Pattern_combined.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt Pattern_combined.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt \\\$LSB_JOBINDEX $DELTA_LINES Core_patterns.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str} $((${#n_bar_groups}+1))\"")
          echo $str
          printf "\n"
          out=$(eval "$str")
          job_id=$(echo $out | awk '{print gensub(/(.+<)([0-9]*)(>.+)/,"\\2","g",$0)}')

          str=$(echo "bsub -P SJLIFE -J core_pattern_comparisons.${NAMES[j]}.sleep -w \"done(${job_id})\" -R \"rusage[mem=128]\" -R \"status=='ok'\" -oo ${file_str}/core_pattern_comparisons.${NAMES[j]}.sleep.out -eo ${file_str}/core_pattern_comparisons.${NAMES[j]}.sleep.err -K \"sleep 10\"")
          echo $str
          eval "$str"
          printf "\n"

          # conatenate comparison jobs into a single Core_patterns file
          echo Write ${NAMES[j]} core patterns at level $((iteration-1)) in Core_patterns.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          printf "\n"
          f_ind=1
          for file in Core_patterns.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.*.txt
          do
            echo Found file $file
            if (($f_ind==1))
            then
              if [ -f $file ];
              then
                echo "cat $file > Core_patterns.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt"
                cat $file > Core_patterns.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
                f_ind=$((f_ind+1))
                echo rm $file
                rm $file
                echo New file index is $f_ind
                printf "\n"
              else
                echo Create empty file Core_patterns.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
                touch Core_patterns.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
              fi
            else
              if [ -f $file ];
              then
                echo "cat $file >> Core_patterns.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt"
                cat $file >> Core_patterns.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
                f_ind=$((f_ind+1))
                echo rm $file
                rm $file
                echo New file index is $f_ind
                printf "\n"
              fi
            fi
          done

          f_limit=$((20*$(du -b Pattern_combined.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt | cut -f1)/1000000))
          if (( $f_limit>10000 )) # if requested memory is greater than 10 GB, use large_mem queue
          then
            queue="large_mem"
          else
            queue="standard"
          fi

          N_LINES=$(awk 'END{print NR}' Pattern_combined.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt)
          echo Total number of lines in Pattern_combined.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt is $N_LINES.
          DELTA_LINES=200000
          COMPARISON_JOBS=$(echo $DELTA_LINES | awk -v x=$N_LINES '{printf "%0.25f\n", x/$1}')
          COMPARISON_JOBS=$(echo $COMPARISON_JOBS | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
          if (($COMPARISON_JOBS < 1));
          then
            COMPARISON_JOBS=1
          fi
          echo Number of jobs, based on $DELTA_LINES line$( (($DELTA_LINES > 1)) && echo "s" || echo "" ) per job, is $COMPARISON_JOBS.
          printf "\n"

          # Patterns appearing in the new iteration (TEST_FILE, input $1) not included in the previous iteration (REF_FILE, input $2), CLOSED PATTERNS
          str=$(echo "bsub -P SJLIFE -J \"myJob[1-$COMPARISON_JOBS]\" -oo ${file_str_next}/closed_pattern_comparisons.%I.${NAMES[j]}.${file_str_next}.out -eo ${file_str_next}/closed_pattern_comparisons.%I.${NAMES[j]}.${file_str_next}.err -R \"rusage[mem=$f_limit]\" -R \"status=='ok'\" -q $queue \"sh ${HOME_DIR}/file_compare.sh Pattern_combined.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt Pattern_combined.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt \\\$LSB_JOBINDEX $DELTA_LINES Closed_patterns.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str} $((${#n_bar_groups}+1))\"")
          echo $str
          out=$(eval "$str")
          echo $out
          job_id=$(echo $out | awk '{print gensub(/(.+<)([0-9]*)(>.+)/,"\\2","g",$0)}')
          printf "\n"

          str=$(echo "bsub -P SJLIFE -J closed_pattern_comparisons.${NAMES[j]}.sleep -w \"done(${job_id})\" -R \"rusage[mem=128]\" -R \"status=='ok'\" -oo ${file_str_next}/closed_pattern_comparisons.${NAMES[j]}.sleep.out -eo ${file_str_next}/closed_pattern_comparisons.${NAMES[j]}.sleep.err -K \"sleep 10\"")
          echo $str
          eval "$str"
          printf "\n"

          # conatenate comparison jobs into a single Closed_patterns file
          echo Write ${NAMES[j]} closed patterns at level $((iteration)) in Closed_patterns.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          printf "\n"
          f_ind=1
          for file in Closed_patterns.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.*.txt
          do
            echo Found file $file
            if (($f_ind==1))
            then
              if [ -f $file ];
              then
                echo "cat $file > Closed_patterns.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt"
                cat $file > Closed_patterns.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
                f_ind=$((f_ind+1))
                echo rm $file
                rm $file
                echo New file index is $f_ind
                printf "\n"
              else
                echo Create empty file Closed_patterns.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
                touch Closed_patterns.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
                printf "\n"
              fi
            else
              if [ -f $file ];
              then
                echo "cat $file >> Closed_patterns.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt"
                cat $file >> Closed_patterns.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
                f_ind=$((f_ind+1))
                echo rm $file
                rm $file
                echo New file index is $f_ind
                printf "\n"
              fi
            fi
          done
        fi
      done

      if (($iteration==1)) # initialization of tables of closed and core patterns across all iterations in the first iteration
      then
        # Order in which populations are sampled
        idx_closed=1; # whether "Closed_patterns" files exist for each population
        idx_core=1; # whether "Core_patterns" files exist for each population
        if ((${#NAMES[@]}>1)) # only create stats file for combined population if more than one population
        then
          for ((i=0;i<${#NAMES[@]};i++))
          do
            if [ ! -f Closed_patterns.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt ]
            then
              echo Closed_patterns.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt missing for population ${NAMES[i]}
              idx_closed=0 # at least one population missing "Closed_patterns" for this pattern
              printf "\n"
            fi
            if [ ! -f Core_patterns.${file_str}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt ]
            then
              echo Core_patterns.${file_str}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt missing for population ${NAMES[i]}
              idx_core=0 # at least one population missing "Closed_patterns" for this pattern
              printf "\n"
            fi

            echo idx_closed = $idx_closed.  File Closed_patterns.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt does$( (($idx_closed == 0)) && echo " not" || echo "") exist.
            echo idx_core = $idx_core.  File Core_patterns.${file_str}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt does$( (($idx_core == 0)) && echo " not" || echo "") exist.
            printf "\n"
          done
        fi

        for ((j=0;j<${#NAMES[@]};j++)) # create _pattern_stats files for each population
        do
          (($idx_closed==1)) && echo Initiate file Closed_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt:
          echo test -f Closed_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt \&\& rm Closed_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          test -f Closed_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt && rm Closed_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          (($idx_closed==1)) && echo touch Closed_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          (($idx_closed==1)) && touch Closed_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          printf "\n"

          (($idx_core==1)) && echo Initiate file Core_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          echo test -f Core_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt \&\& rm Core_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          test -f Core_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt && rm Core_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          (($idx_core==1)) && echo touch Core_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          (($idx_core==1)) && touch Core_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          printf "\n"
        done
      fi

      # Check counts of controls closed patterns among cases and controls, stratified by iteration
      if ls Closed_patterns_stats.*.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt > /dev/null 2>&1 # https://stackoverflow.com/questions/6363441/check-if-a-file-exists-with-wildcard-in-shell-script
      then
        # instances of closed patterns
        for ((i=0;i<${#NAMES[@]};i++)) # loop over all possibilities for pop1 and append to list of counts of patterns in pop1 in all populations
        do
          if [ -s Closed_patterns.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt ] #if file exists and is non-empty https://stackoverflow.com/questions/9964823/how-to-check-if-a-file-is-empty-in-bash
          then
            idx_closed=1 # check whether Closed_patterns file exists for this population and pattern
          else
            idx_closed=0 # or not
          fi
          echo idx_closed=$idx_closed. File Closed_patterns.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt $( (($idx_closed == 1)) && echo " exists and is non-empty." || echo "does not exist.")
          printf "\n"
        done

        lines_array=()
        mem_array=()
        for pop in ${NAMES[@]}
        do
          if [ -f Closed_patterns.${file_str_next}.${pop}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt ];
          then
            lines_array+=($(awk 'END{print NR}' Closed_patterns.${file_str_next}.${pop}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt))
            mem_array+=($((2*$(du -b Closed_patterns.${file_str_next}.${pop}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt | cut -f1)/1000000)))
          fi
        done
        echo Lines in Closed_patterns files:
        declare -p lines_array
        printf "\n"

        N_LINES=0
        f_limit=10
        for n in ${lines_array[@]} #looping through line counts in Closed_patterns files for different populations
        do
          (( $n > $N_LINES )) && N_LINES=$n
        done
        for f in ${mem_array[@]}
        do
          (( $f > $f_limit )) && f_limit=$f
        done

        echo Maximum number of lines in the Closed_patterns files is $N_LINES.
        printf "\n"
        DELTA_LINES=200000
        COMPARISON_JOBS=$(echo $DELTA_LINES | awk -v x=$N_LINES '{printf "%0.25f\n", x/$1}') # exact number of comparison jobs
        COMPARISON_JOBS=$(echo $COMPARISON_JOBS | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
        if (($COMPARISON_JOBS < 1));
        then
          COMPARISON_JOBS=1
        fi
        echo Number of jobs, based on $DELTA_LINES line$( (($DELTA_LINES > 1)) && echo "s" || echo "" ) per job, is $COMPARISON_JOBS.
        printf "\n"

        # uncomment to evaluate haplotype counts of each closed pattern
        # if [ -f ${POPULATION[j]} ]
        # then
        #   str=$(echo "bsub -P SJLIFE -J \"myJob[1-$COMPARISON_JOBS]\" -oo ${file_str_next}/closed_pattern_counts.%I.out -eo ${file_str_next}/closed_pattern_counts.%I.err -R \"rusage[mem=$f_limit]\" \"sh ${HOME_DIR}/population_haplotype_counts.sh ${NAMES[j]} ${POPULATION[j]} Closed_patterns.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt \\\$LSB_JOBINDEX $DELTA_LINES $iteration\"")
        #   echo $str
        #   printf "\n"
        #   out=$(eval "$str")
        #   echo $out
        #   job_id=$(echo $out | awk '{print gensub(/(.+<)([0-9]*)(>.+)/,"\\2","g",$0)}')
        #   printf "\n"
        #   str=$(echo "bsub -P SJLIFE -w \"done(${job_id})\" -R \"rusage[mem=128]\" -oo ${file_str_next}/closed_pattern_counts.sleep.out -eo ${file_str_next}/closed_pattern_counts.sleep.err -K \"sleep 10\"")
        #   echo $str
        #   eval "$str"
        #   printf "\n"
        # fi

        # Comment out if evaluating Closed pattern counts in previous block
        for ((j=0;j<${#NAMES[@]};j++))
        do
          # instead of counting occurences of pattern, print 0 in the second field
          echo Write iteration $iteration patterns to Closed_patterns_stats for ${NAMES[j]}:
          echo awk \'BEGIN{OFS=\"\\t\"} {print ${iteration},0,\$0}\' Closed_patterns.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt \>\>  Closed_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          awk 'BEGIN{OFS="\t"} {print '${iteration}',0,$0}' Closed_patterns.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt >>  Closed_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          printf "\n"
        done

        # Get closed patterns for each population in a separate file
        for ((j=0;j<${#NAMES[@]};j++))
        do
          for file in Closed_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.*.txt
          do
            echo $file
            if [ -f $file ];
            then
              echo "cat $file >> Closed_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt"
              cat $file >> Closed_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
              echo rm $file
              rm $file
              printf "\n"
            fi
          done
        done
      fi

      if ls Core_patterns_stats.*.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt > /dev/null 2>&1
      then
        #Instances of core patterns
        for ((i=0;i<${#NAMES[@]};i++)) # loop over all possibilities for pop1 and append to list of counts of patterns in pop1 in all populations
        do
          if [ -s Core_patterns.${file_str}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt ] # if file exists and is non-empty https://stackoverflow.com/questions/9964823/how-to-check-if-a-file-is-empty-in-bash
          then
            idx_core=1 # check whether Core_patterns file exists for this population and pattern
          else
            idx_core=0 # or not
          fi
          echo idx_core=$idx_core. File Core_patterns.${file_str}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt $( (($idx_core == 1)) && echo " exists and is non-empty." || echo "does not exist.")
          printf "\n"
        done

        #get maximum number of lines in the Core_patterns files for this iteration across all populations
        lines_array=()
        mem_array=()
        for pop in ${NAMES[@]}
        do
          if [ -f Core_patterns.${file_str}.${pop}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt ];
          then
            lines_array+=($(awk 'END{print NR}' Core_patterns.${file_str}.${pop}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt))
            mem_array+=($((20*$(du -b Core_patterns.${file_str}.${pop}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt | cut -f1)/1000000)))
          fi
        done
        echo Lines in Core_patterns files:
        declare -p lines_array
        printf "\n"

        N_LINES=0
        f_limit=10
        for n in ${lines_array[@]} #looping through line counts in Closed_patterns files for different populations
        do
          (( $n > $N_LINES )) && N_LINES=$n
        done
        for f in ${mem_array[@]}
        do
          (( $f > $f_limit )) && f_limit=$f
        done

        echo Maximum number of lines in the Core_patterns files is $N_LINES.
        printf "\n"
        DELTA_LINES=20000
        COMPARISON_JOBS=$(echo $DELTA_LINES | awk -v x=$N_LINES '{printf "%0.25f\n", x/$1}') #exact number of comparison jobs
        COMPARISON_JOBS=$(echo $COMPARISON_JOBS | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
        if (($COMPARISON_JOBS < 1));
        then
          COMPARISON_JOBS=1
        fi
        echo Number of jobs, based on $DELTA_LINES line$( (($DELTA_LINES > 1)) && echo "s" || echo "" ) per job, is $COMPARISON_JOBS.
        printf "\n"

        # uncomment to evaluate haplotype counts of each core pattern
        # if [ -f ${POPULATION[j]} ]
        # then
        #   str=$(echo "bsub -P SJLIFE -J \"myJob[1-$COMPARISON_JOBS]\" -oo ${file_str}/core_pattern_counts.%I.out -eo ${file_str}/core_pattern_counts.%I.err -R \"rusage[mem=$f_limit]\" \"sh ${HOME_DIR}/population_haplotype_counts.sh ${NAMES[j]} ${POPULATION[j]} Core_patterns.${file_str_next}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt \\\$LSB_JOBINDEX $DELTA_LINES $iteration\"")
        #   echo $str
        #   printf "\n"
        #   out=$(eval "$str")
        #   echo $out
        #   job_id=$(echo $out | awk '{print gensub(/(.+<)([0-9]*)(>.+)/,"\\2","g",$0)}')
        #   str=$(echo "bsub -P SJLIFE -w \"done(${job_id})\" -R \"rusage[mem=128]\" -oo ${file_str}/core_pattern_counts.sleep.out -eo ${file_str}/core_pattern_counts.sleep.err -K \"sleep 10\"")
        #   echo $str
        #   eval "$str"
        #   printf "\n"
        # fi

        # Comment out if evaluating Core pattern counts in previous block
        for ((j=0;j<${#NAMES[@]};j++))
        do
          # instead of counting occurences of pattern, print 0 in the second field
          echo Write iteration $iteration patterns to Core_patterns_stats for ${NAMES[j]}:
          echo awk \'BEGIN{OFS=\"\\t\"} {print $((iteration-1)),0,\$0}\' Core_patterns.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt \>\> Core_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          awk 'BEGIN{OFS="\t"} {print '$((iteration-1))',0,$0}' Core_patterns.${file_str}.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt >> Core_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          printf "\n"
        done

        # Get core patterns for cases only in a single file
        for ((j=0;j<${#NAMES[@]};j++))
        do
          for file in Core_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
          do
            echo Found file $file
            if [ -f $file ];
            then
              if ((${#POPULATION[@]}>1))
              then
                echo Compute Fisher\'s exact test for patterns:
                echo bsub \-P SJLIFE \-J ${NAMES[j]}.ukbb_haplotype_partition_sub.${file_str} \-oo ${NAMES[j]}.ukbb_haplotype_partition_sub.${file_str}.out \-eo ${NAMES[j]}.ukbb_haplotype_partition_sub.${file_str}.err \-R \"rusage[mem=256]\" \-R \"status==\'ok\'\" \-K \"sh ${HOME_DIR}/ukbb_haplotype_partition_sub.sh ${POPULATION[0]},${POPULATION[1]} Core_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt ${CHR} ${BP_RANGE[0]},${BP_RANGE[1]} 50 $((iteration-1))\"
                job_id=$(bsub -P SJLIFE -J ${NAMES[j]}.ukbb_haplotype_partition_sub.${file_str} -oo ${NAMES[j]}.ukbb_haplotype_partition_sub.${file_str}.out -eo ${NAMES[j]}.ukbb_haplotype_partition_sub.${file_str}.err -R "rusage[mem=256]" -R "status=='ok'" -K "sh ${HOME_DIR}/ukbb_haplotype_partition_sub.sh ${POPULATION[0]},${POPULATION[1]} Core_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt ${CHR} ${BP_RANGE[0]},${BP_RANGE[1]} 50 $((iteration-1))")
                echo $job_id
                job_id=$(echo $job_id | awk 'b=gensub(/.*<([0-9]*)>.*/,"\\1","g",$0) {print b}') #extract job_id (number) from output
                printf "\n"

                echo Wait until all p values calculated:
                echo bsub \-P SJLIFE \-J sleep.${NAMES[j]}.ukbb_haplotype_partition_sub.${file_str} \-w \"done\($job_id\)\" \-R \"rusage[mem=32]\" \-R \"status==\'ok\'\" \-oo ${DIRECTORY}/sleep.${NAMES[j]}.ukbb_haplotype_partition_sub.${file_str}.out \-eo ${DIRECTORY}/sleep.${NAMES[j]}.ukbb_haplotype_partition_sub.${file_str}.err \-K \"sleep 10\"
                bsub -P SJLIFE -J sleep.${NAMES[j]}.ukbb_haplotype_partition_sub.${file_str} -w "done($job_id)" -R "rusage[mem=32]" -R "status=='ok'" -oo ${DIRECTORY}/sleep.${NAMES[j]}.ukbb_haplotype_partition_sub.${file_str}.out -eo ${DIRECTORY}/sleep.${NAMES[j]}.ukbb_haplotype_partition_sub.${file_str}.err -K "sleep 10"
                printf "\n"

                echo awk \'\$1==$((iteration-1)){count+=1} END{print count-1}\' Core_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
                total_patterns=$(awk '$1=='$((iteration-1))'{count+=1} END{print count-1}' Core_patterns_stats.${NAMES[j]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt)
                test -f fisher_exact.patterns_0-${total_patterns}.txt && echo mv fisher_exact.patterns_0-${total_patterns}.txt fisher_exact.${NAMES[j]}.${file_str}.txt && printf "\n"
                test -f fisher_exact.patterns_0-${total_patterns}.txt && mv fisher_exact.patterns_0-${total_patterns}.txt fisher_exact.${NAMES[j]}.${file_str}.txt

                echo Plot results for ${file_str}:
                echo Rscript ${HOME_DIR}/ukbb_fisher_exact_plot.R file=fisher_exact.${NAMES[j]}.${file_str}.txt
                Rscript ${HOME_DIR}/ukbb_fisher_exact_plot.R file=fisher_exact.${NAMES[j]}.${file_str}.txt
              fi
            fi
          done
        done
      fi
    fi
  done

  # update counts for next iteration
  DELTA=$5 # restart DELTA at its original value
  N_JOBS=$(awk 'BEGIN{printf "%0.25f\n",'$N_OVERLAPS'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
  echo Requested number of jobs for $N_OVERLAPS pattern$( (($N_OVERLAPS>1)) && echo "s" || echo "" ) is $N_JOBS with DELTA = $DELTA pattern$( (($DELTA>1)) && echo "s" || echo "" ) per job.
  while (( $N_JOBS > $MAX_JOBS ))
  do
    DELTA=$((2*DELTA))
    N_JOBS=$(awk 'BEGIN{printf "%0.25f\n",'$N_OVERLAPS'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
    echo Revised number of jobs is $N_JOBS with DELTA = $DELTA haplotype$( (($DELTA>1)) && echo "s" || echo "" ) per job.
  done
  printf "\n"

  mem_array=($(for folder in ${file_str}*; #Loop through directories starting with Iteration000
  do
    cd $folder
    for file in overlaps*.out
    do
      if [ -f $file ]
      then
        mem=$(head -100 $file | grep "Max Memory " | tr -dc '0-9');
        echo $mem
      fi
    done
    cd ..
  done))

  max_mem=256; # used at most this much in last round
  for mem in ${mem_array[@]}
  do
    if (( $mem > $max_mem ))
    then
      max_mem=$mem
    fi
  done
  mem_req=$((2*max_mem)) #increase to new limit

  # assumes each step of each job uses 10% of size (in mb) of input Pattern file just written
  mem_array=($(for file in Pattern_combined.${file_str_next}*
  do
    if [ -f $file ];
    then
      echo $((($DELTA/10)*$(du -b $file | cut -f1)/1000000));
    fi
  done))
  declare -p mem_array
  for mem in ${mem_array[@]}
  do
    if (( $mem > $max_mem ))
    then
      max_mem=$mem
    fi
  done
  mem_req=$((max_mem))

  echo Memory request per job is $mem_req MB
  printf "\n"
done

echo $file_str_next
printf "\n"
for ((i=0;i<${#NAMES[@]};i++))
do
  echo ${NAMES[i]}
  printf "\n"
  for pattern_str in ${suffix_array[@]}
  do
    echo $pattern_str
    printf "\n"
    if [ -f Pattern_combined.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt ]
    then
      echo cat Pattern_combined.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt \| wc -l
      cat Pattern_combined.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt | wc -l
      printf "\n"
      if (( $(cat Pattern_combined.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt | wc -l)==1 ))
      then
        echo Get results for last pattern:
        echo Last pattern is a core pattern:
        echo awk \'BEGIN{OFS=\"\\t\"} {for \(j=2\;j\<=NF\;j++\) {printf \"%s%s\",\$j,\(j==NF?\"\\n\":\"\\t\"\)} }\' Pattern_combined.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt \> Core_patterns.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
        awk 'BEGIN{OFS="\t"} {for (j=2;j<=NF;j++) {printf "%s%s",$j,(j==NF?"\n":"\t")} }' Pattern_combined.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt > Core_patterns.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
        printf "\n"

        echo Write iteration $iteration patterns to Core_patterns_stats for ${NAMES[i]}:
        echo awk \'BEGIN{OFS=\"\\t\"} {print $iteration,0,\$0}\' Core_patterns.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt \>\> Core_patterns_stats.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
        awk 'BEGIN{OFS="\t"} {print '$iteration',0,$0}' Core_patterns.${file_str_next}.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt >> Core_patterns_stats.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
        printf "\n"

        if [ -f Core_patterns_stats.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt ];
        then
          if ((${#POPULATION[@]}>1))
          then
            echo Compute Fisher\'s exact test for patterns:
            echo bsub \-P SJLIFE \-J ${NAMES[i]}.ukbb_haplotype_partition_sub.${file_str_next} \-oo ${NAMES[i]}.ukbb_haplotype_partition_sub.${file_str_next}.out \-eo ${NAMES[i]}.ukbb_haplotype_partition_sub.${file_str_next}.err \-R \"rusage[mem=256]\" \-R \"status==\'ok\'\" \-K \"sh ${HOME_DIR}/ukbb_haplotype_partition_sub.sh ${POPULATION[0]},${POPULATION[1]} Core_patterns_stats.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt ${CHR} ${BP_RANGE[0]},${BP_RANGE[1]} 50 $iteration\"
            job_id=$(bsub -P SJLIFE -J ${NAMES[i]}.ukbb_haplotype_partition_sub.${file_str_next} -oo ${NAMES[j]}.ukbb_haplotype_partition_sub.${file_str_next}.out -eo ${NAMES[i]}.ukbb_haplotype_partition_sub.${file_str_next}.err -R "rusage[mem=256]" -R "status=='ok'" -K "sh ${HOME_DIR}/ukbb_haplotype_partition_sub.sh ${POPULATION[0]},${POPULATION[1]} Core_patterns_stats.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt ${CHR} ${BP_RANGE[0]},${BP_RANGE[1]} 50 $iteration")
            echo $job_id
            job_id=$(echo $job_id | awk 'b=gensub(/.*<([0-9]*)>.*/,"\\1","g",$0) {print b}') #extract job_id (number) from output
            printf "\n"

            echo Wait until all p values calculated:
            echo bsub \-P SJLIFE \-J sleep.${NAMES[i]}.ukbb_haplotype_partition_sub.${file_str_next} \-w \"done\($job_id\)\" \-R \"rusage[mem=32]\" \-R \"status==\'ok\'\" \-oo ${DIRECTORY}/sleep.${NAMES[i]}.ukbb_haplotype_partition_sub.${file_str_next}.out \-eo ${DIRECTORY}/sleep.${NAMES[i]}.ukbb_haplotype_partition_sub.${file_str_next}.err \-K \"sleep 10\"
            bsub -P SJLIFE -J sleep.${NAMES[i]}.ukbb_haplotype_partition_sub.${file_str_next} -w "done($job_id)" -R "rusage[mem=32]" -R "status=='ok'" -oo ${DIRECTORY}/sleep.${NAMES[i]}.ukbb_haplotype_partition_sub.${file_str_next}.out -eo ${DIRECTORY}/sleep.${NAMES[i]}.ukbb_haplotype_partition_sub.${file_str_next}.err -K "sleep 10"
            printf "\n"

            echo awk \'$1==$iteration{count+=1} END{print count-1}\' Core_patterns_stats.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt
            total_patterns=$(awk '$1=='$iteration'{count+=1} END{print count-1}' Core_patterns_stats.${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${pattern_str}.txt)
            test -f fisher_exact.patterns_0-${total_patterns}.txt && echo mv fisher_exact.patterns_0-${total_patterns}.txt fisher_exact.${NAMES[i]}.${file_str_next}.txt && printf "\n"
            test -f fisher_exact.patterns_0-${total_patterns}.txt && mv fisher_exact.patterns_0-${total_patterns}.txt fisher_exact.${NAMES[i]}.${file_str_next}.txt

            echo Plot results for ${file_str_next}:
            echo Rscript ${HOME_DIR}/ukbb_fisher_exact_plot.R file=fisher_exact.${NAMES[i]}.${file_str_next}.txt
            Rscript ${HOME_DIR}/ukbb_fisher_exact_plot.R file=fisher_exact.${NAMES[i]}.${file_str_next}.txt
          fi
        fi
      fi
    fi
  done
done
