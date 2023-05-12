#! /bin/bash
set -e

# Iteratively overlaps the new closed patterns in Pattern_combined.Iteration000.NAME_PATTERN.txt until no new patterns have appeared

# sh HOME_DIR/ChromosomeOverlap_iteration_sub_parallel.v2.sh chr11.69231642-69431642 2 2,j "50" "DIRECTORY" "HOME_DIR"

NAME=$1 # optional identifier in Pattern_combined file
SIGMA=$2 # number of patterns to be overlapped in one comparison job
PATTERN=$3 # specify a suffix to test a specfic combination, e.g., 2,3,j or 2+3+j; else uses all
DELTA=$4 # number of combinations to do in one job
ITERATION=$5 # integer 0,1,...: iteration to start at
DIRECTORY=$6 # location of haplotype matrices (snps x subjects); a full file path or relative ./
HOME_DIR=$7 # location of program files, subject_generate_haplotypes.sh; a full file path

if [ -z $HOME_DIR ]
then
  unset HOME_DIR
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
echo Working directory is $DIRECTORY
printf "\n"
cd $DIRECTORY

if [ ! -z $NAME ]
then
  NAME=($(echo $NAME | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))
else
  (>&2 echo "No population name supplied."; exit 1)
fi

patterns=()
for file in Pattern_combined.Iteration$(printf "%0.3d" $ITERATION).* # Pattern files from 0th ITERATION, each has a pattern suffix e.g., 2,3,j
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

if [ -z $SIGMA ];
then
  SIGMA=2 # number of patterns in an overlap
fi

if [ -z $MAX_JOBS ]
then
  MAX_JOBS=2000 #revised limit
fi

if [ -z $DELTA ]
then
  DELTA=50
fi

if [ -z $ITERATION ]
then
  ITERATION=0; # only 0th round of overlaps has been completed
fi
echo Initial Iteration string:
file_str=$(printf "Iteration%0.3d" $ITERATION)
echo $file_str
printf "\n"

echo Pattern_combined files found for ${PATTERN}:
for file in Pattern_combined.${file_str}.*${NAME}*${PATTERN}.txt
do
  test -f $file && echo $file && found+=1
done
test ! -z $found && (( $found>0 )) && found=0 && printf "\n" || (>&2 echo "No Iteration$(printf "%0.3d" $ITERATION) file found."; exit 1)

n_overlaps=() # initialize array of overlaps to be performed, based on number of lines in Pattern_combined files for each population
for file in Pattern_combined.${file_str}.*${NAME}*${PATTERN}.txt
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
  n_overlaps+=($n_tuples) # add to array of counts
  n_input_files=${#n_overlaps[@]} # total number of SNP files. Should be the same from iteration to iteration
done
printf "\n"
N_OVERLAPS=0
echo Array of tuples to be performed in each file:
declare -p n_overlaps
printf "\n"
# get maximum in array n_overlaps
for n in ${n_overlaps[@]}
do
  if (( $n > $N_OVERLAPS ));
  then
    N_OVERLAPS=$n # replace old maximum
  fi
done
echo Maximum number of overlaps to be performed is $N_OVERLAPS.
printf "\n"

N_JOBS=$(awk 'BEGIN{printf "%0.25f\n",'$N_OVERLAPS'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
echo Requested number of jobs for $N_OVERLAPS combination$( (($N_OVERLAPS>1)) && echo "s" || echo "" ) is $N_JOBS with DELTA = $DELTA pattern$( (($DELTA>1)) && echo "s" || echo "" ) per job.
(( $N_JOBS > 1000000000 )) && MAX_JOBS=$((MAX_JOBS*2)) || MAX_JOBS=2000
while (( $N_JOBS > $MAX_JOBS ))
do
  DELTA=$((2*DELTA))
  N_JOBS=$(awk 'BEGIN{printf "%0.25f\n",'$N_OVERLAPS'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
  echo Revised number of jobs is $N_JOBS with DELTA = $DELTA pattern$( (($DELTA>1)) && echo "s" || echo "" ) per job.
done
MAX_JOBS=2000
printf "\n"

while (($N_OVERLAPS>0)) # continue as long as there is more than one pattern to be overlapped with another
do
  file_str=$(printf "Iteration%0.3d" $ITERATION) # from the last round
  ITERATION=$(($ITERATION+1)) # starts at 0, so first round creates iteration.Step001
  file_str_next=$(printf "Iteration%0.3d" $ITERATION) # in this round

  # Create subdirectory with necessary files
  char=$(($(echo $N_JOBS | wc -m)-1)) # the number of characters in the total number of jobs
  SUBDIR=$(printf "%sIteration%0.3d.Step%0.${char}d" ${DIRECTORY/%//} $ITERATION 1) # full path to SUBDIR
  # ITERATION=$(printf "Iteration%0.3d" $ITERATION) #incremented by 1 from variable file_str
  STEP=$(printf "Step%0.${char}d" 1)

  test -d $SUBDIR && echo rm -r $SUBDIR
  test -d $SUBDIR && rm -r $SUBDIR
  test ! -d $SUBDIR && echo mkdir $SUBDIR && printf "\n"
  test ! -d $SUBDIR && mkdir $SUBDIR

  echo Moving copy of Pattern_combined.${file_str}* to $SUBDIR
  for file in ${DIRECTORY/%//}Pattern_combined.${file_str}*${NAME}*${PATTERN}.txt
  do
    echo cp $file ${SUBDIR/%//}${file##*/}
    cp $file ${SUBDIR/%//}${file##*/} && printf "\n"

    echo Transpose Patterns_combined file and move to $SUBDIR:
    TRANSPOSE_FILE=${file%.*}.transpose.${file##*.}

    mem_req=$(du $file | cut -f1 | awk '{sum+=$1} END{print int(3*sum/1024)}') # 1.5 times total size, in MB, of Overlap_tuples files
    if (( $((mem_req+0)) < 1024 ))
    then
      mem_req=5000
    fi
    if (( $((mem_req+0)) > 1000000 ))
    then
      mem_req=1000000
    fi
    echo bsub \-P SJLIFE \-J transpose.${file_str} \-oo transpose.${file_str}.out \-eo transpose.${file_str}.err \-R \"rusage[mem=$mem_req]\" \-K \"awk \'BEGIN{OFS=\\\"\\\\t\\\"}\; {for\(j=1\;j\<=NF\;j++\) {a[NR,j]=\\\$j\; n_rows=NR\; n_cols=\(n_cols\<NF?NF:n_cols\)} } END{for \(j=1\;j\<=n_cols\;j++\) {for \(i=1\;i\<=n_rows\;i++\) {printf \\\"%s%s\\\",a[i,j],\(i==n_rows?\\\"\\\\n\\\":\\\"\\\\t\\\"\)} } }\' $file \> ${SUBDIR/%//}${TRANSPOSE_FILE##*/}\"
    bsub -P SJLIFE -J transpose.${file_str} -oo transpose.${file_str}.out -eo transpose.${file_str}.err -R "rusage[mem=$mem_req]" -K "awk 'BEGIN{OFS=\"\\t\"}; {for(j=1;j<=NF;j++) {a[NR,j]=\$j; n_rows=NR; n_cols=(n_cols<NF?NF:n_cols)} } END{for (j=1;j<=n_cols;j++) {for (i=1;i<=n_rows;i++) {printf \"%s%s\",a[i,j],(i==n_rows?\"\\n\":\"\\t\")} } }' $file > ${SUBDIR/%//}${TRANSPOSE_FILE##*/}" && printf "\n"

    mem_req=5000
  done

  for file in ${SUBDIR/%//}Pattern_combined.${file_str}*${NAME}*${PATTERN}.txt
  do

    echo bsub \-P SJLIFE \-J \"myJob[1-$N_JOBS]\" \-eo ${SUBDIR/%//}ChromosomeOverlap_iteration.v2.${file##*/}.%I.err \-oo ${SUBDIR/%//}ChromosomeOverlap_iteration.v2.${file##*/}.%I.out \-R \"rusage[mem=5000]\" \-R \"order[!ut]\" \-R \"select[ut \< 0.9 \&\& status==\'ok\']\" \"sh ${HOME_DIR/%//}ChromosomeOverlap_iteration.v2.sh \\\"\\\" ${SUBDIR/%//}${file##*/} ${SUBDIR/%//}${TRANSPOSE_FILE##*/} ${DELTA}.\\\$LSB_JOBINDEX 2500000 $SIGMA ${file_str_next} $DIRECTORY $HOME_DIR\"
    job_id=$(bsub -P SJLIFE -J "myJob[1-$N_JOBS]" -eo ${SUBDIR/%//}ChromosomeOverlap_Iteration.v2.${file##*/}.%I.err -oo ${SUBDIR/%//}ChromosomeOverlap_iteration.v2.${file##*/}.%I.out -R "rusage[mem=5000]" -R "order[!ut]" -R "select[ut < 0.9 && status=='ok']" "sh ${HOME_DIR/%//}ChromosomeOverlap_iteration.v2.sh \"\" ${SUBDIR/%//}${file##*/} ${SUBDIR/%//}${TRANSPOSE_FILE##*/} ${DELTA}.\$LSB_JOBINDEX 2500000 $SIGMA ${file_str_next} $DIRECTORY $HOME_DIR")
    echo $job_id && printf "\n"
    job_id=$(echo $job_id | awk 'b=gensub(/.*<([0-9]*)>.*/,"\\1","g",$0) {print b}') #extract job_id (number) from output

    mem_req=5000

  done

  echo Wait until all jobs start:
  echo bsub \-P SJLIFE \-J sleep.ChromosomeOverlap_iteration.v3.${file_str_next} \-w \"numrun\($job_id,*\) \|\| numended\($job_id,*\)\" \-R \"rusage[mem=32]\" \-R \"status==\'ok\'\" \-oo ${DIRECTORY/%//}sleep.ChromosomeOverlap_iteration.v3.${file_str_next}.out \-eo ${DIRECTORY/%//}sleep.ChromosomeOverlap_iteration.v3.${file_str_next}.err \-K \"sleep 10\"
  bsub -P SJLIFE -J sleep.ChromosomeOverlap_iteration.v3.${file_str_next} -w "numrun($job_id,*) || numended($job_id,*)" -R "rusage[mem=32]" -R "status=='ok'" -oo ${DIRECTORY/%//}sleep.ChromosomeOverlap_iteration.v3.${file_str_next}.out -eo ${DIRECTORY/%//}sleep.ChromosomeOverlap_iteration.v3.${file_str_next}.err -K "sleep 10"
  printf "\n"

  echo Wait until all jobs done:
  job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}')) # get job name from 7th field (or 6th if no exectution host yet) in all non-header rows of bjobs
  declare -p job_array
  while (( ${#job_array[@]}>0 ))
  do
    sleep 20
    job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}'))
  done
  printf "\n"

  str=$(echo "bsub -P SJLIFE -J sleep.${file_str_next} -R \"rusage[mem=32]\" -R \"status=='ok'\" -eo ${SUBDIR/%//}sleep.err -oo ${SUBDIR/%//}sleep.out -K \"sleep 10\"")
  printf "\n"
  echo $str
  eval $str
  printf "\n"

  echo cd $DIRECTORY
  cd $DIRECTORY && printf "\n"

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

  echo Patterns found among overlaps at iteration ${ITERATION}:
  declare -p suffix_array
  printf "\n"

  mem_req=$(du Overlap_tuples* | cut -f1 | awk '{sum+=$1} END{print int(1.5*sum/1024)}') # 1.5 times total size, in MB, of Overlap_tuples files
  echo $mem_req Mb && printf "\n"
  if (( $((mem_req+0)) < 5000 ))
  then
    mem_req=5000
  fi
  if (( $((mem_req+0)) > 1000000 ))
  then
    mem_req=1000000
  fi
  # Combine patterns from differnt jobs and move to DIRECTORY as new Pattern_combined file
  str=$(echo "bsub -P SJLIFE -J ChromosomeOverlap_iteration_combine.${NAME}.${file_str_next} -R \"rusage[mem=$mem_req]\" -R \"status=='ok'\" -oo combine.${NAME}.${file_str_next}.out -eo combine.${NAME}.${file_str_next}.err -K \"${HOME_DIR/%//}ChromosomeOverlap_iteration_combine.sh ${NAME} $file_str_next "$PATTERN" $DIRECTORY\"" )
  echo $str
  eval $str
  printf "\n"

  mem_req=5000


  test -f $TRANSPOSE_FILE && echo rm $TRANSPOSE_FILE || printf ""
  test -f $TRANSPOSE_FILE && rm $TRANSPOSE_FILE && printf "\n" || printf ""

  # Find closed patterns not appearing in the next iteration
  echo Array of patterns:
  declare -p patterns
  printf "\n"
  for pattern_str in ${patterns[*]}
  do
    if ls Pattern_combined.${file_str_next}*${NAME}*_${pattern_str}.txt 1> /dev/null 2> /dev/null # only evaluate pattern_str (e.g., 2,j) if it has associated Pattern_combined files
    then
      n_bar_groups=$(($(echo ${pattern_str} | tr -c -d "+" | wc -c)+1))
      echo Number of groups is $n_bar_groups

      # Get the second column (snp pattern) of the new file (file_str_next) and print lines of the first file that do not have a match in the second column

      echo Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt
      printf "\n"
      if [ -f Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt ];
      then

        N_LINES=$(awk 'END{print NR}' Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt)
        echo Total number of lines in Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt is $N_LINES.
        DELTA_LINES=2000000
        COMPARISON_JOBS=$(echo $DELTA_LINES | awk -v x=$N_LINES '{printf "%0.25f\n", x/$1}')
        COMPARISON_JOBS=$(echo $COMPARISON_JOBS | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
        if (($COMPARISON_JOBS < 1));
        then
          COMPARISON_JOBS=1
        fi
        echo Number of jobs, based on $DELTA_LINES line$( (($DELTA_LINES > 1)) && echo "s" || echo "" ) per job, is $COMPARISON_JOBS.
        printf "\n"



        for ((i=$5;i<${ITERATION};i++)) # from the starting iteration up to the last iteration
        do
          file_str_temp=$(printf "Iteration%0.3d" $i)

          # Patterns appearing in the new iteration (TEST_FILE, input $1) not included in the previous iteration (REF_FILE, input $2)
          mem_req=$(awk 'BEGIN{print int(1.5*'$(du Pattern_combined.${file_str_temp}.${NAME/%/_}${pattern_str}.txt | cut -f1)'/1024)}') # 1.5 times the size of REF_FILE, in MB
          if (( $((mem_req+0)) < 1024 ))
          then
            mem_req=5000
          fi
          # Create Pattern_combined file with only newly appearing patterns
          str=$(echo "bsub -P SJLIFE -J myJob[1-$COMPARISON_JOBS] -oo compare.${file_str_next}-${file_str_temp}.${NAME/%/_}${pattern_str}.%I.out -eo compare.${file_str_next}-${file_str_temp}.${NAME/%/_}${pattern_str}.%I.err -R \"rusage[mem=$mem_req]\" \"sh ${HOME_DIR/%//}file_compare.v2.sh Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt Pattern_combined.${file_str_temp}.${NAME/%/_}${pattern_str}.txt \\\$LSB_JOBINDEX $DELTA_LINES Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str} $((${#n_bar_groups}+1))\"") # delete patterns in TEST_FILE that have appeared in each earlier REF_FILE
          echo $str
          eval $str && printf "\n"
          mem_req=5000

          echo Wait until all jobs done:
          job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}')) # get job name from 7th field (or 6th if no exectution host yet) in all non-header rows of bjobs
          declare -p job_array
          while (( ${#job_array[@]}>0 ))
          do
            sleep 20
            job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}'))
          done
          printf "\n"
          sleep 10

          # conatenate comparison jobs into a single Closed_patterns file
          echo Write ${NAME} closed patterns appearing at iteration $((ITERATION)) in Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt
          printf "\n"
          f_ind=1
          for file in Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.*.txt
          do
            echo Found file $file
            if (($f_ind==1))
            then
              if [ -f $file ];
              then
                echo wc -l $file
                wc -l $file && printf "\n"
                echo "cat $file > Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt"
                cat $file > Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt && printf "\n"
                f_ind=$((f_ind+1))
                echo rm $file
                rm $file && printf "\n"
                echo New file index is $f_ind
                printf "\n"
              else
                echo Create empty file Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt
                touch Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt
                printf "\n"
              fi
            else
              if [ -f $file ];
              then
                echo wc -l $file
                wc -l $file && printf "\n"
                echo "cat $file >> Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt"
                cat $file >> Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt && printf "\n"
                f_ind=$((f_ind+1))
                echo rm $file
                rm $file && printf "\n"
                echo New file index is $f_ind
                printf "\n"
              fi
            fi
          done
        done
      fi

      # get updated number of jobs
      # get maximum number of patterns in one of the Pattern_combined files
      n_overlaps=()
      for file in Pattern_combined.${file_str_next}.*${NAME}*${PATTERN}.txt
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
          n_overlaps+=($n_tuples) # add to array of counts
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

      if (($ITERATION==1)) # initialization of tables of closed patterns across all itrerations in the first iteration
      then

        echo Initiate file Closed_patterns_stats.${NAME/%/_}${pattern_str}.txt:
        echo test -f Closed_patterns_stats.${NAME/%/_}${pattern_str}.txt \&\& rm Closed_patterns_stats.${NAME/%/_}${pattern_str}.txt
        test -f Closed_patterns_stats.${NAME/%/_}${pattern_str}.txt && rm Closed_patterns_stats.${NAME/%/_}${pattern_str}.txt
        echo touch Closed_patternse_stats.${NAME/%/_}${pattern_str}.txt
        touch Closed_patterns_stats.${NAME/%/_}${pattern_str}.txt
        printf "\n"

        if [ -f Pattern_combined.${file_str}.${NAME/%/_}${pattern_str}.txt ]
        then
          echo Write iteration $((ITERATION-1)) patterns to Closed_patterns_stats for ${NAME}:
          echo awk \'BEGIN{OFS=\"\\t\"} {print $((ITERATION-1)),\$0}\' Pattern_combined.${file_str}.${NAME/%/_}${pattern_str}.txt \>\> Closed_patterns_stats.${NAME/%/_}${pattern_str}.txt
          awk 'BEGIN{OFS="\t"} {print '$((ITERATION-1))',$0}' Pattern_combined.${file_str}.${NAME/%/_}${pattern_str}.txt >> Closed_patterns_stats.${NAME/%/_}${pattern_str}.txt
          printf "\n"
        fi

      fi

      # Check counts of controls closed patterns among cases and controls, stratified by iteration
      if ls Closed_patterns_stats.*_${pattern_str}.txt > /dev/null 2>&1
      then
        if [ -f Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt ]
        then
          echo Write iteration $ITERATION patterns to Closed_patterns_stats for ${NAME}:
          echo awk \'BEGIN{OFS=\"\\t\"} {print ${ITERATION},\$0}\' Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt \>\> Closed_patterns_stats.${NAME/%/_}${pattern_str}.txt
          awk 'BEGIN{OFS="\t"} {print '${ITERATION}',$0}' Pattern_combined.${file_str_next}.${NAME/%/_}${pattern_str}.txt >> Closed_patterns_stats.${NAME/%/_}${pattern_str}.txt
          printf "\n"
        fi
      fi
    fi
  done

  # update counts for next iteration
  if [ ! -z $4 ]
  then
    DELTA=$4 # restart DELTA at its original value
  else
    DELTA=50
  fi
  N_JOBS=$(awk 'BEGIN{printf "%0.25f\n",'$N_OVERLAPS'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
  echo Requested number of jobs for $N_OVERLAPS pattern$( (($N_OVERLAPS>1)) && echo "s" || echo "" ) is $N_JOBS with DELTA = $DELTA pattern$( (($DELTA>1)) && echo "s" || echo "" ) per job.
  (( $N_JOBS > 1000000000 )) && MAX_JOBS=$((MAX_JOBS*2)) || MAX_JOBS=2000
  while (( $N_JOBS > $MAX_JOBS ))
  do
    DELTA=$((2*DELTA))
    N_JOBS=$(awk 'BEGIN{printf "%0.25f\n",'$N_OVERLAPS'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
    echo Revised number of jobs is $N_JOBS with DELTA = $DELTA haplotype$( (($DELTA>1)) && echo "s" || echo "" ) per job.
  done
  MAX_JOBS=2000
  printf "\n"
  if (($ITERATION==1)) # stop (and regroup) if last iteration was the initiation step 000
  then
    exit
  fi
done
