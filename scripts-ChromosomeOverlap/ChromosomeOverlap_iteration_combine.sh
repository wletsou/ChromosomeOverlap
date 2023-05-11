#! /bin/bash
set -e

# Combines Overlap_tuples files from ITERATION

# sh HOME_DIR/ChromosomeOverlap_iteration_combine.sh chr11.69231642-69431642 "Iteration001" "2,j" DIRECTORY

NAME=$1 # optional identifier for combined file
ITERATION=$2 # string "Iteration00x" used to form SUBDIR
PATTERNS=$3 # optional colon-separated list of pattern types to overlap.  If empty, use all patterns found
DIRECTORY=$4 # where output combined files will be stored

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
echo Directory is $DIRECTORY
printf "\n"

if [ -z $ITERATION ];
then
  (>&2 echo "Iteration string must be supplied"; exit)
fi

if [ -z $NAME ]
then
  unset NAME
fi

unset pattern_array
if [ -z $PATTERNS ]; # get all patterns if no pattern supplied
then
  echo Currently in $PWD
  echo Combining files:
  for file in Overlap_tuples.${ITERATION}.*${NAME}*.txt
  do
    if [ -f $file ]
    then
      echo $file
    fi
  done
  # trim pattern from end file and store in array
  pattern_array=($(for file in Overlap_tuples.${ITERATION}.*${NAME}*.txt
  do
    if [ -f $file ];
    then
      pattern=${file##*_}
      echo ${pattern%.txt}
    fi
  done | uniq))
else
  pattern_array=($(echo $PATTERNS | perl -pne 's/([^:]+)[:]*/$1 /g'))
fi
declare -p pattern_array && printf "\n"

# combine files of the same pattern type
for ((i=0; i<${#pattern_array[@]}; i++))
do
  n_bars=$(echo ${pattern_array[i]} | tr -cd "+" | wc -c)
  let n_groups=$n_bars+2 # first group is in second field, and last is $n_bar fields after
  # initiate combined file, and change name from "Overlap" to "Pattern"
  test -f Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt && echo rm Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt
  test -f Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt && rm Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt
  echo touch Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt
  touch Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt && printf "\n"
  for file in Overlap_tuples.${ITERATION}.*.${NAME/%/_}${pattern_array[i]}.txt
  do
    if [ -f $file ];
    then
      echo cat $file \>\> Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt
      cat $file >> Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt
      echo rm $file
      rm $file
    fi
  done
  printf "\n"

  echo Combined file Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt has $(awk '{total+=$1} END{print total}' Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt) entries in $(awk 'END{print NR}' Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt) lines.

  str=$(echo "awk 'BEGIN{OFS=\"\t\"} { array[sprintf(\"%s")
  for ((k=3; k<=$n_groups; k++));
  do
    str=${str}$(echo "\t%s")
  done
  str=${str}$(echo "\"")
  for ((k=2; k<=$n_groups; k++));
  do
    str=${str}$(echo ",\$$k")
  done
  SORT="TRUE" # whether to sort output by pattern
  SORT="FALSE" # whether to sort output by pattern
  if [ $SORT == "TRUE" ];
  then
    sort_str="unsorted."
  else
    sort_str=""
  fi
  str=${str}$(echo ")]+=\$1 } END{ for (row in array) { print array[row],row } }' ${DIRECTORY}/Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt > ${DIRECTORY}/Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.${sort_str}tmp ")
  echo $str
  eval "$str"
  printf "\n"

  test -f ${DIRECTORY}/Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.${sort_str}tmp && echo mv ${DIRECTORY}/Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.${sort_str}tmp ${DIRECTORY}/Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.${sort_str}txt
  test -f ${DIRECTORY}/Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.${sort_str}tmp && mv ${DIRECTORY}/Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.${sort_str}tmp ${DIRECTORY}/Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.${sort_str}txt && printf "\n"


  if [ $SORT == "TRUE" ];
  then
    if [ -f Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.${sort_str}txt ]
    then
      sort -k2 Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.${sort_str}txt > Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt
    fi
    echo Output file Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt has $(awk '{total+=$1} END{print total}' Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt) entries in $(awk 'END{print NR}' Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.txt) lines after sorting.

    for file in Pattern_combined.${ITERATION}.${NAME}*${sort_str}txt
    do
      if [ -f $file ];
      then
        echo rm $file
        rm $file
      fi
    done
  else
    test -f ${DIRECTORY}/Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.${sort_str}txt && echo Output file Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.${sort_str}txt has $(awk '{total+=$1} END{print total}' Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.${sort_str}txt) entries in $(awk 'END{print NR}' Pattern_combined.${ITERATION}.${NAME/%/_}${pattern_array[i]}.${sort_str}txt) lines with no sorting.
  fi
done

# remove unsorted combined file
