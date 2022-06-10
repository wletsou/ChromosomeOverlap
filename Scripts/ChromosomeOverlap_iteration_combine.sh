#! /bin/bash

POPULATION=$1 # single population name or combined, not a list
#SUBDIR=$4 #Index the iteration round by "iteration00x" to keep track of which Pattern_combined.iteration00x files are being processed
ITERATION=$2 # string "Iteration00x" used to form SUBDIR
PATTERNS=$3 # optional colon-separated list of pattern types to overlap.  If empty, use all patterns found
DIRECTORY=$4 # where output combined files will be stored

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
echo Directory is $DIRECTORY
printf "\n"

# create a subdir if it doesn't already exist for storing Pattern files after they have been combined
if [ -z $ITERATION ];
then
  SUBDIR="Iteration001"
else
  SUBDIR=${ITERATION}
fi
# remove previous instances of the Pattern files associated with this population
if [ ! -d $SUBDIR ];
then
  mkdir_ind=1
  echo mkdir ${DIRECTORY}/$SUBDIR
  mkdir ${DIRECTORY}/$SUBDIR
fi
test ! -z $mkdir_ind && printf "\n"

if [ "${PWD##*/}" == "$SUBDIR" ];
then
  # remove combined files from other calls to script, but not combined files from last iteration
  for file in Pattern_combined.${ITERATION}.*.txt
  do
  file_ind=1
    if [ -f $file ];
    then
      echo rm $file
      rm $file
    fi
  done
  test ! -z $file_ind && printf "\n"
fi
echo Subdirectory is $SUBDIR

unset pattern_array
if [ -z $PATTERNS ]; # get all patterns if no pattern supplied
then
  echo Currently in $PWD
  echo Combining files:
  for file in Overlap_tuples.${ITERATION}.*.${POPULATION}*.txt
  do
    if [ -f $file ]
    then
      echo $file
    fi
  done
  # trim pattern from end file and store in array
  pattern_array=($(for file in Overlap_tuples.${ITERATION}.*.${POPULATION}*.txt
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
declare -p pattern_array

# combine files of the same pattern type
for ((i=0; i<${#pattern_array[@]}; i++))
do
  n_bars=$(echo ${pattern_array[i]} | tr -cd "+" | wc -c)
  let n_groups=$n_bars+2 # first group is in second field, and last is $n_bar fields after
  # initiate combined file, and change name from "Overlap" to "Pattern"
  touch Pattern_combined.${ITERATION}.${POPULATION}_${pattern_array[i]}.txt
  for file in Overlap_tuples.${ITERATION}.*.${POPULATION}_${pattern_array[i]}.txt
  do
    if [ -f $file ];
    then
      echo cat $file \>\> Pattern_combined.${ITERATION}.${POPULATION}_${pattern_array[i]}.txt
      cat $file >> Pattern_combined.${ITERATION}.${POPULATION}_${pattern_array[i]}.txt
      mv $file ${DIRECTORY}/${SUBDIR}/${file} # store Overlap_tuples files in SUBDIR after combining in Pattern_combined
    fi
  done
  printf "\n"

  echo Combined file Pattern_combined.${ITERATION}.${POPULATION}_${pattern_array[i]}.txt has $(awk '{total+=$1} END{print total}' Pattern_combined.${ITERATION}.${POPULATION}_${pattern_array[i]}.txt) entries in $(awk 'END{print NR}' Pattern_combined.${ITERATION}.${POPULATION}_${pattern_array[i]}.txt) lines.

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
  if [ $SORT == "TRUE" ];
  then
    sort_str="unsorted."
  else
    sort_str=""
  fi
  str=${str}$(echo ")]+=\$1 } END{ for (row in array) { print array[row],row } }' ${DIRECTORY}/Pattern_combined.${ITERATION}.${POPULATION}_${pattern_array[i]}.txt > ${DIRECTORY}/Pattern_combined.${ITERATION}.${POPULATION}_${pattern_array[i]}.${sort_str}txt ")
  echo $str
  eval "$str"
  printf "\n"


  if [ $SORT == "TRUE" ];
  then
    if [ -f Pattern_combined.${ITERATION}.${POPULATION}_${pattern_array[i]}.${sort_str}txt ]
    then
      sort -k2 Pattern_combined.${ITERATION}.${POPULATION}_${pattern_array[i]}.${sort_str}txt > Pattern_combined.${ITERATION}.${POPULATION}_${pattern_array[i]}.txt
    fi
  fi

  echo Output file Pattern_combined.${ITERATION}.${POPULATION}_${pattern_array[i]}.txt has $(awk '{total+=$1} END{print total}' Pattern_combined.${ITERATION}.${POPULATION}_${pattern_array[i]}.txt) entries in $(awk 'END{print NR}' Pattern_combined.${ITERATION}.${POPULATION}_${pattern_array[i]}.txt) lines after sorting.
done

# remove unsorted combined file
for file in Pattern_combined.${ITERATION}.${POPULATION}*${sort_str}txt
do
  if [ -f $file ];
  then
    echo rm $file
    rm $file
  fi
done
