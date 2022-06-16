#! /bin/bash
set -e

# Iteratively overlaps the patterns in Pattern_combined.Iteration000.NAME_PATTERN.txt until all patterns have disappeared

# sh HOME_DIR/ChromosomeOverlap_iteration_sub.sh chr11.69231642-69431642 2 2,j "50" "DIRECTORY" "HOME_DIR"

NAME=$1 # optional identifier in Pattern_combined file
SIGMA=$2 # number of patterns to be overlapped in one comparison job
PATTERN=$3 # specify a suffix to test a specfic combination, e.g., 2,3,j or 2+3+j; else uses all
DELTA=$4 # number of combinations to do in one job
DIRECTORY=$5 # location of haplotype matrices (snps x subjects); a full file path or relative ./
HOME_DIR=$6 # location of program files, subject_generate_haplotypes.sh; a full file path

if [ -z $HOME_DIR ];
then
  HOME_DIR="" # in case HOME_DIR is on your PATH
else
  HOME_DIR="${HOME_DIR}/"
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

if [ -z $SIGMA ];
then
  # SIGMA=1 # number of patterns in an overlap
  SIGMA=2 # number of patterns in an overlap
fi

if [ -z $MAX_JOBS ]
then
  MAX_JOBS=100 #revised limit
fi

if [ -z $DELTA ]
then
  DELTA=50
fi

iteration=0; # only 0th round of overlaps has been completed
echo Pattern_combined files found for ${PATTERN}:
for file in Pattern_combined.*${NAME}*${PATTERN}.txt
do
  test -f $file && echo $file && found+=1
done
test ! -z $found && (( $found>0 )) && found=0 && printf "\n"

echo Initial iteration string:
file_str=$(printf "Iteration%0.3d" $iteration)
echo $file_str
printf "\n"

n_overlaps=() #initialize array of overlaps to be performed, based on number of lines in Pattern_combined files for each population
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
echo Requested number of jobs for $N_OVERLAPS combination$( (($N_OVERLAPS>1)) && echo "s" || echo "" ) is $N_JOBS with DELTA = $DELTA pattern$( (($DELTA>1)) && echo "s" || echo "" ) per job.
while (( $N_JOBS > $MAX_JOBS ))
do
  DELTA=$((2*DELTA))
  N_JOBS=$(awk 'BEGIN{printf "%0.25f\n",'$N_OVERLAPS'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
  echo Revised number of jobs is $N_JOBS with DELTA = $DELTA haplotype$( (($DELTA>1)) && echo "s" || echo "" ) per job.
done
printf "\n"

while (($N_OVERLAPS>0)) # continue as long as there is more than one pattern to be overlapped with another
do
  file_str=$(printf "Iteration%0.3d" $iteration) # from the last round
  iteration=$(($iteration+1)) # starts at 0, so first round creates Iteration001.Step001
  file_str_next=$(printf "Iteration%0.3d" $iteration) #in this round

  # Create subdirectory with necessary files
  char=$(($(echo $N_JOBS | wc -m)-1)) # the number of characters in the total number of jobs
  SUBDIR=$(printf "%s/Iteration%0.3d.Step%0.${char}d" ${DIRECTORY} $iteration 1) # full path to SUBDIR
  ITERATION=$(printf "Iteration%0.3d" $iteration) #incremented by 1 from variable file_str
  STEP=$(printf "Step%0.${char}d" 1)

  test -d $SUBDIR && echo rm -r $SUBDIR
  test -d $SUBDIR && rm -r $SUBDIR
  test ! -d $SUBDIR && echo mkdir $SUBDIR && printf "\n"
  test ! -d $SUBDIR && mkdir $SUBDIR

  echo Moving copy of Pattern_combined.${file_str}* to $SUBDIR
  for file in ${DIRECTORY}/Pattern_combined.${file_str}*${NAME}*${PATTERN}.txt
  do
    echo cp $file ${SUBDIR}/${file##*/}
    cp $file ${SUBDIR}/${file##*/} && printf "\n"

    echo Transpose Patterns_combined file: # used by pattern_overlap_loop5.sh
    TRANSPOSE_FILE=${file%.*}.transpose.${file##*.}
    echo awk \'BEGIN{OFS=\"\\t\"}\; {for\(j=1\;j\<=NF\;j++\) {a[NR,j]=\$j\; n_rows=NR\; n_cols=\(n_cols\<NF?NF:n_cols\)} } END{for \(j=1\;j\<=n_cols\;j++\) {for \(i=1\;i\<=n_rows\;i++\) {printf \"%s%s\",a[i,j],\(i==n_rows?\"\\n\":\"\\t\"\)} } }\' $file \> $TRANSPOSE_FILE
    awk 'BEGIN{OFS="\t"}; {for(j=1;j<=NF;j++) {a[NR,j]=$j; n_rows=NR; n_cols=(n_cols<NF?NF:n_cols)} } END{for (j=1;j<=n_cols;j++) {for (i=1;i<=n_rows;i++) {printf "%s%s",a[i,j],(i==n_rows?"\n":"\t")} } }' $file > $TRANSPOSE_FILE # transpose to counts in first row, snp patterns by bar groups in subsequent rows
  done

  for file in ${SUBDIR}/Pattern_combined.${file_str}*${NAME}*${PATTERN}.txt
  do
    # or submit jobs sequentially
    for ((i=1;i<=$N_JOBS;i++))
    do
      echo sh ${HOME_DIR}ChromosomeOverlap_iteration.sh ${SUBDIR}/${file##*/} ${DELTA}.$i 1000000 $SIGMA $ITERATION $DIRECTORY $HOME_DIR
      sh ${HOME_DIR}ChromosomeOverlap_iteration.sh ${SUBDIR}/${file##*/} ${DELTA}.$i 1000000 $SIGMA $ITERATION $DIRECTORY $HOME_DIR && printf "\n"
    done

  done
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

  echo Patterns found among overlaps at iteration ${iteration}:
  declare -p suffix_array
  printf "\n"

  # Combine patterns from differnt jobs and move to DIRECTORY as new Pattern_combined file
  echo ${HOME_DIR}ChromosomeOverlap_iteration_combine.sh ${NAME} $file_str_next \"$PATTERN\" $DIRECTORY
  sh ${HOME_DIR}ChromosomeOverlap_iteration_combine.sh ${NAME} $file_str_next "$PATTERN" $DIRECTORY
  printf "\n"

  # get updated number of jobs
  # get maximum number of patterns in one of the Pattern_combined files
  n_overlaps=()
  echo Looking for $n_input_files file$((($n_input_files!=1)) && echo "s" || echo ""). # number corresponds to number of populations
  printf "\n"
  while ((${#n_overlaps[@]} < $n_input_files)) # wait until all input files have been combined, or count the remaining files again
  do
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

  echo Remove subdirectories
  test -d $SUBDIR && echo rm -r $SUBDIR
  test -d $SUBDIR && rm -r $SUBDIR || echo No directories removed.
  printf "\n"

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

      # Get the second column (snp pattern) of the new file (file_str_next) and print lines of the first file that do not have a match in the second column https://stackoverflow.com/questions/15251188/find-the-difference-between-two-files

      echo Pattern_combined.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
      printf "\n"
      if [ -f Pattern_combined.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt ];
      then

        N_LINES=$(awk 'END{print NR}' Pattern_combined.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt)
        echo Total number of lines in Pattern_combined.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt is $N_LINES.
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

        # patterns in the old iteration (TEST_FILE, input $1) not in the new iteration (REF_FILE, input $2)
        for ((i=1;i<=$COMPARISON_JOBS;i++))
        do
          echo sh ${HOME_DIR}file_compare.sh Pattern_combined.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt Pattern_combined.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt $i $DELTA_LINES Closed_patterns_post.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str} $((${#n_bar_groups}+1))
          sh ${HOME_DIR}file_compare.sh Pattern_combined.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt Pattern_combined.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt $i $DELTA_LINES Closed_patterns_post.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str} $((${#n_bar_groups}+1)) && printf "\n"
        done

        # conatenate comparison jobs into a single Closed_patterns_post file
        echo Write ${NAME} closed patterns disappearing at iteration $((iteration-1)) in Closed_patterns_post.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        printf "\n"
        f_ind=1
        for file in Closed_patterns_post.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.*.txt
        do
          echo Found file $file
          if (($f_ind==1))
          then
            if [ -f $file ];
            then
              echo "cat $file > Closed_patterns_post.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt"
              cat $file > Closed_patterns_post.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
              f_ind=$((f_ind+1))
              echo rm $file
              rm $file
              echo New file index is $f_ind
              printf "\n"
            else
              echo Create empty file Closed_patterns_post.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
              touch Closed_patterns_post.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
            fi
          else
            if [ -f $file ];
            then
              echo "cat $file >> Closed_patterns_post.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt"
              cat $file >> Closed_patterns_post.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
              f_ind=$((f_ind+1))
              echo rm $file
              rm $file
              echo New file index is $f_ind
              printf "\n"
            fi
          fi
        done

        N_LINES=$(awk 'END{print NR}' Pattern_combined.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt)
        echo Total number of lines in Pattern_combined.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt is $N_LINES.
        DELTA_LINES=200000
        COMPARISON_JOBS=$(echo $DELTA_LINES | awk -v x=$N_LINES '{printf "%0.25f\n", x/$1}')
        COMPARISON_JOBS=$(echo $COMPARISON_JOBS | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
        if (($COMPARISON_JOBS < 1));
        then
          COMPARISON_JOBS=1
        fi
        echo Number of jobs, based on $DELTA_LINES line$( (($DELTA_LINES > 1)) && echo "s" || echo "" ) per job, is $COMPARISON_JOBS.
        printf "\n"

        # Patterns appearing in the new iteration (TEST_FILE, input $1) not included in the previous iteration (REF_FILE, input $2)

        for ((i=1;i<=$COMPARISON_JOBS;i++))
        do
          echo sh ${HOME_DIR}file_compare.sh Pattern_combined.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt Pattern_combined.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt $i $DELTA_LINES Closed_patterns_pre.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str} $((${#n_bar_groups}+1))
          sh ${HOME_DIR}file_compare.sh Pattern_combined.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt Pattern_combined.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt $i $DELTA_LINES Closed_patterns_pre.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str} $((${#n_bar_groups}+1)) && printf "\n"
        done

        # conatenate comparison jobs into a single Closed_patterns_pre file
        echo Write ${NAME} closed patterns appearing at iteration $((iteration)) in Closed_patterns_pre.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        printf "\n"
        f_ind=1
        for file in Closed_patterns_pre.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.*.txt
        do
          echo Found file $file
          if (($f_ind==1))
          then
            if [ -f $file ];
            then
              echo "cat $file > Closed_patterns_pre.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt"
              cat $file > Closed_patterns_pre.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
              f_ind=$((f_ind+1))
              echo rm $file
              rm $file
              echo New file index is $f_ind
              printf "\n"
            else
              echo Create empty file Closed_patterns_pre.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
              touch Closed_patterns_pre.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
              printf "\n"
            fi
          else
            if [ -f $file ];
            then
              echo "cat $file >> Closed_patterns_pre.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt"
              cat $file >> Closed_patterns_pre.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
              f_ind=$((f_ind+1))
              echo rm $file
              rm $file
              echo New file index is $f_ind
              printf "\n"
            fi
          fi
        done
      fi

      if (($iteration==1)) # initialization of tables of closed patterns across all iterations in the first iteration
      then
        # Order in which populations are sampled
        idx_pre=1; # whether "Closed_patterns_pre" files exist for each population
        idx_post=1; # whether "Closed_patterns_post" files exist for each population

        echo Closed patterns from initial step, iteration $((iteration-1)):
        for ((i=1;i<=$COMPARISON_JOBS;i++))
        do
          echo sh ${HOME_DIR}file_compare.sh Pattern_combined.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt /dev/null $i $DELTA_LINES Closed_patterns_pre.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str} $((${#n_bar_groups}+1))
          sh ${HOME_DIR}file_compare.sh Pattern_combined.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt /dev/null $i $DELTA_LINES Closed_patterns_pre.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str} $((${#n_bar_groups}+1)) && printf "\n"  # lines in Pattern_combined.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt that are not in the empty file /dev/null
        done

        echo Write ${NAME} closed patterns appearing at iteration $((iteration-1)) in Closed_patterns_pre.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        printf "\n"
        f_ind=1
        for file in Closed_patterns_pre.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.*.txt
        do
          echo Found file $file
          if (($f_ind==1))
          then
            if [ -f $file ];
            then
              echo "cat $file > Closed_patterns_pre.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt"
              cat $file > Closed_patterns_pre.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
              f_ind=$((f_ind+1))
              echo rm $file
              rm $file
              echo New file index is $f_ind
              printf "\n"
            else
              echo Create empty file Closed_patterns_pre.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
              touch Closed_patterns_pre.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
              printf "\n"
            fi
          else
            if [ -f $file ];
            then
              echo "cat $file >> Closed_patterns_pre.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt"
              cat $file >> Closed_patterns_pre.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
              f_ind=$((f_ind+1))
              echo rm $file
              rm $file
              echo New file index is $f_ind
              printf "\n"
            fi
          fi
        done

        (($idx_pre==1)) && echo Initiate file Closed_patterns_pre_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt:
        echo test -f Closed_patterns_pre_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt \&\& rm Closed_patterns_pre_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        test -f Closed_patterns_pre_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt && rm Closed_patterns_pre_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        (($idx_pre==1)) && echo touch Closed_patterns_pre_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        (($idx_pre==1)) && touch Closed_patterns_pre_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        printf "\n"

        if [ -f Closed_patterns_pre.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt ]
        then
          echo Write iteration $((iteration-1)) patterns to Closed_patterns_pre_stats for ${NAME}:
          echo awk \'BEGIN{OFS=\"\\t\"} {print $((iteration-1)),0,\$0}\' Closed_patterns_pre.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt \>\>  Closed_patterns_pre_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
          awk 'BEGIN{OFS="\t"} {print '$((iteration-1))',0,$0}' Closed_patterns_pre.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt >>  Closed_patterns_pre_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
          printf "\n"
        fi

        (($idx_post==1)) && echo Initiate file Closed_patterns_post_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        echo test -f Closed_patterns_post_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt \&\& rm Closed_patterns_post_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        test -f Closed_patterns_post_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt && rm Closed_patterns_post_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        (($idx_post==1)) && echo touch Closed_patterns_post_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        (($idx_post==1)) && touch Closed_patterns_post_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        printf "\n"
      fi

      # Check counts of controls closed patterns among cases and controls, stratified by iteration
      if ls Closed_patterns_pre_stats.*_${pattern_str}.txt > /dev/null 2>&1 # https://stackoverflow.com/questions/6363441/check-if-a-file-exists-with-wildcard-in-shell-script
      then
        # instances of closed patterns
        if [ -s Closed_patterns_pre.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt ] # if file exists and is non-empty https://stackoverflow.com/questions/9964823/how-to-check-if-a-file-is-empty-in-bash
        then
          idx_pre=1 # check whether Closed_patterns_pre file exists for this population and pattern
        else
          idx_pre=0 # or not
        fi
        echo idx_pre=$idx_pre. File Closed_patterns_pre.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt $( (($idx_pre == 1)) && echo " exists and is non-empty." || echo "does not exist.")
        printf "\n"

        echo Write iteration $iteration patterns to Closed_patterns_pre_stats for ${NAME}:
        echo awk \'BEGIN{OFS=\"\\t\"} {print ${iteration},0,\$0}\' Closed_patterns_pre.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt \>\>  Closed_patterns_pre_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        awk 'BEGIN{OFS="\t"} {print '${iteration}',0,$0}' Closed_patterns_pre.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt >>  Closed_patterns_pre_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        printf "\n"
      fi

      if ls Closed_patterns_post_stats.*_${pattern_str}.txt > /dev/null 2>&1
      then
        # Instances of closed patterns
        if [ -s Closed_patterns_post.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt ] # if file exists and is non-empty https://stackoverflow.com/questions/9964823/how-to-check-if-a-file-is-empty-in-bash
        then
          idx_post=1 # check whether Closed_patterns_post file exists for this population and pattern
        else
          idx_post=0 # or not
        fi
        echo idx_post=$idx_post. File Closed_patterns_post.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt $( (($idx_post == 1)) && echo " exists and is non-empty." || echo "does not exist.")
        printf "\n"

        # instead of counting occurences of pattern, print 0 in the second field
        echo Write iteration $((iteration-1)) patterns to Closed_patterns_post_stats for ${NAME}:
        echo awk \'BEGIN{OFS=\"\\t\"} {print $((iteration-1)),0,\$0}\' Closed_patterns_post.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt \>\> Closed_patterns_post_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        awk 'BEGIN{OFS="\t"} {print '$((iteration-1))',0,$0}' Closed_patterns_post.${file_str}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt >> Closed_patterns_post_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
        printf "\n"
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
  while (( $N_JOBS > $MAX_JOBS ))
  do
    DELTA=$((2*DELTA))
    N_JOBS=$(awk 'BEGIN{printf "%0.25f\n",'$N_OVERLAPS'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
    echo Revised number of jobs is $N_JOBS with DELTA = $DELTA haplotype$( (($DELTA>1)) && echo "s" || echo "" ) per job.
  done
  printf "\n"
done

for pattern_str in ${suffix_array[@]}
do
  if [ -f Pattern_combined.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt ]
  then
    echo cat Pattern_combined.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt \| wc -l
    cat Pattern_combined.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt | wc -l
    printf "\n"
    if (( $(cat Pattern_combined.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt | wc -l)==1 ))
    then
      echo Get results for last pattern:
      echo Last pattern is a closed pattern:
      echo awk \'BEGIN{OFS=\"\\t\"} {for \(j=2\;j\<=NF\;j++\) {printf \"%s%s\",\$j,\(j==NF?\"\\n\":\"\\t\"\)} }\' Pattern_combined.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt \> Closed_patterns_post.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
      awk 'BEGIN{OFS="\t"} {for (j=2;j<=NF;j++) {printf "%s%s",$j,(j==NF?"\n":"\t")} }' Pattern_combined.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt > Closed_patterns_post.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
      printf "\n"

      echo Write iteration $iteration patterns to Closed_patterns_post_stats for ${NAME}:
      echo awk \'BEGIN{OFS=\"\\t\"} {print $iteration,0,\$0}\' Closed_patterns_post.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt \>\> Closed_patterns_post_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
      awk 'BEGIN{OFS="\t"} {print '$iteration',0,$0}' Closed_patterns_post.${file_str_next}.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt >> Closed_patterns_post_stats.$( [ -z $NAME ] && echo "" || echo "${NAME}_")${pattern_str}.txt
      printf "\n"
    fi
  fi
done
