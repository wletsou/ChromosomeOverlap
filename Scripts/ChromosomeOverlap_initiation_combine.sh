#! /bin/bash

# Combine files of the same overlap type and total the number of unique patterns

# /home/wletsou/scripts/ChromosomeOverlap_initiation_combine.sh chr11.69231642-69431642 1 "Iteration000" "" /scratch_space/wletsou/sjlife/GWAS/ChromosomeOverlapTest

POPULATION=$1 # single population name whose tuples are to be combined
COMBINE_GROUPS=$2 # whether to merge files with the same number of bar groups, e.g., 2,3+j 2,j+3 and 2+3,j form a two-group pattern
NAME=$3 # Default iteration name, Iteration000
DIRECTORY=$4

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
echo Output directory is $DIRECTORY
# Do NOT cd $DIRECTORY, because POPULATION files are in the submission directory, not the output directory

if [ -z $NAME ];
then
  NAME="Iteration000" # label for Pattern_combined files
fi
# create folder to store all files with the same pattern
if [ ! -z $NAME ];
then
  SUBDIR=${NAME}
  NAME=".${NAME}"
else
  SUBDIR="Iteration000"
  NAME=".${SUBDIR}"
fi

# create subdirectory for original files to be moved after combinations complete
test ! -d $SUBDIR && echo mkdir ${DIRECTORY}/${SUBDIR} && printf "\n"
test ! -d $SUBDIR && mkdir ${DIRECTORY}/${SUBDIR}

test -d $SUBDIR && echo cd $SUBDIR && printf "\n"
test -d $SUBDIR && cd $SUBDIR
if [ "${PWD##*/}" == "${SUBDIR}" ]; # if current directory (after removing up to rightmost slash) is the same as the subdirectory, i.e., if successfully in the SUBDIR, then remove other Pattern files
then
  for file in Pattern*.${POPULATION}.*.txt
  do
    if [ -f $file ]
    then
      test -f $file && found+=1 # counter whether any files have been removed
      test -f $file && echo rm $file
      test -f $file && rm $file
    fi
  done
  test ! -z $found && (( $found>0 )) && found=0 && printf "\n"
fi
test -d $DIRECTORY && echo cd $DIRECTORY && printf "\n"
test -d $DIRECTORY && cd $DIRECTORY

for file in Pattern_combined${NAME}.${POPULATION}* # careful only to delete Pattern_combined files on Iteration000
do
  if [ -f $file ]
  then
    test -f $file && found+=1 # counter whether any files have been removed
    test -f $file && echo rm $file
    test -f $file && rm $file
  fi
done
test ! -z $found && (( $found>0 )) && found=0 && printf "\n"

if [ -z $COMBINE_GROUPS ]
then
  COMBINE_GROUPS=1 # default is to merge files with the same number of bar groups
fi

# find reference pattern, i.e., comma-separated sorted list of pattern with smallest index at each position
unset ref_pattern
unset
echo Files:
ls Pattern.${POPULATION}*.txt
printf "\n"
for file in ${DIRECTORY}/Pattern.${POPULATION}*.txt
do
  pattern=${file##*_} # greedy matching, delete to the left of the last "_" to extract pattern
  echo $pattern
  printf "\n"
  echo $(echo ${pattern%.txt} | perl -pne 's/[+]/,/g')
  pattern_array=($(echo ${pattern%.txt} | perl -pne 's/[+]/,/g' | awk 'END{ split($0,array,","); m=asort(array); for (i=1;i<=m;i++) { printf("%s ",array[i]) } }')) # delete .txt from "pattern"; change all +'s to ,'s; split into a sorted array
  declare -p pattern_array
  printf "\n"

  if [ -z $ref_pattern ];
  then
    ref_pattern=("${pattern_array[@]}") # https://stackoverflow.com/questions/19417015/how-to-copy-an-array-in-bash
  fi
  let n=${#pattern_array[@]} # number of elements in pattern, e.g., 2,3,j has three
  for ((i=0;i<$n;i++))
  do
    if ((${pattern_array[$i]}>${ref_pattern[$i]}));
    then
      break
    fi
    if (($((i+1))==$n)); # smallest indexes possible
    then
      echo New reference pattern ${pattern_array[@]}
      ref_pattern=("${pattern_array[@]}")
      touch Pattern_combined${NAME}.${POPULATION}_${pattern%.txt}.tmp
    fi
  done
done
echo Reference pattern is ${ref_pattern[@]}
printf "\n"

# rename files with the reference pattern
# find reference pattern for files with the same number of bar groups
if (($COMBINE_GROUPS==1))
then
  n_char=${#ref_pattern[@]}
  unset ref_group
  ref_groups=()
  for ((n_pluses=0;n_pluses<$n_char;n_pluses++)) # indicator of number of bars in pattern (one fewer than number of groups)
  do
    for file in ${DIRECTORY}/Pattern.${POPULATION}*.txt
    do
      pattern=${file##*_} # greedy matching, delete to the left of the last "_" to extract pattern
      pattern=${pattern%.txt} # lazy matching, delete first .txt starting from the right
      if (($(echo ${pattern%.txt} | tr -cd "+" | wc -c)==$n_pluses))
      then
        pattern=${pattern%.txt}
        pattern_array=($(echo ${pattern%.txt} | perl -pne 's/[+]/,/g' | awk 'END{ split($0,array,","); m=asort(array); for (i=1;i<=m;i++) { printf("%s ",array[i]) } }'))
        if [ -z $ref_group ] # first pass, no reference group identified
        then
          ref_group=$pattern
          for ((j=0;j<$((${#pattern_array[@]}-1));j++)) # loop through all numeric characters in ref_pattern array, excluding the last which is the letter "j"
          do
            # square of difference between value and position of character ${ref_pattern[j]}
            # https://unix.stackexchange.com/questions/153339/how-to-find-a-position-of-a-character-using-grep/153342
            ref_sum=$(($ref_sum+((${pattern_array[j]}-$(echo $ref_group | grep -aob "\b${pattern_array[$j]}\b" | grep -oP '^[0-9]+')-1))**2))
          done
          # square of difference between length of pattern and position of last "+"
          # https://serverfault.com/questions/197123/getting-the-last-match-in-a-file-using-grep
          ref_sum=$(($ref_sum+((${#pattern}-$(echo "+${ref_group}" | grep -aob "+" | grep -oP '^[0-9]+' | tail -1)))**2))
          echo Reference $((${n_pluses}+1))\-group pattern is $ref_group with total deviation $ref_sum.
          printf "\n"
        fi
        position_sum=0 # after ref_group has been suggested and ref_sum defined
        for ((j=0;j<$((${#pattern_array[@]}-1));j++)) # loop through all numeric characters in pattern_array, excluding the last which is the letter "j"
        do
          # square of difference between value and position of character ${ref_pattern[j]}
          # https://unix.stackexchange.com/questions/153339/how-to-find-a-position-of-a-character-using-grep/153342
          position_sum=$(($position_sum+((${pattern_array[j]}-$(echo $pattern | grep -aob "\b${pattern_array[$j]}\b" | grep -oP '^[0-9]+')-1))**2))
        done
        # square of difference between length of pattern and position of last "+"
        # https://serverfault.com/questions/197123/getting-the-last-match-in-a-file-using-grep
        position_sum=$(($position_sum+((${#pattern}-$(echo "+${pattern}" | grep -aob "+" | grep -oP '^[0-9]+' | tail -1)))**2))
        if (($position_sum<$ref_sum)) # replace ref_(n pluses) if the new sum of square differences of position deviations is smaller than the old
        then
          ref_group=$pattern
          ref_sum=$position_sum
          echo New reference $((${n_pluses}+1))\-group pattern is $ref_group with total deviation $ref_sum.
          printf "\n"
        fi
        unset pattern
        unset pattern_array
      fi
    done
    ref_groups+=($ref_group)
    unset ref_group
  done
fi

for file in ${DIRECTORY}/Pattern.${POPULATION}*.txt
do
  # get indices of pattern and sort; now directly comparable to ref_pattern
  pattern=${file##*_} # remove longest matching string from left (i.e. t left of *)
  pattern_array=($(echo ${pattern%.txt} | perl -pne 's/[+]/,/g' | awk 'END{ split($0,array,","); m=asort(array); for (i=1;i<=m;i++) { printf("%s ",array[i]) } }'))

  # replace in pattern each index based on the ordered comparison between pattern_array and ref_pattern
  pattern=${pattern%.txt}
  echo ${pattern_array[@]}
  printf "\n"
  for ((i=0;i<$n;i++))
  do
    old=$(echo "(\b${pattern_array[i]}\b)(?![A-Za-z0-9])([,+]*)")
    new=${ref_pattern[i]}$(echo "x\$2")
    str=$(echo "echo $pattern | perl -pne 's/$old/$new/g'") # replace the ith element of the pattern the ith element of the reference pattern (with an "x" appended to indicate replacement)
    pattern=$(eval "$str")
  done
  pattern=$(echo $pattern | perl -pne 's/x//g') # delete all x's

  if (($COMBINE_GROUPS==1))
  then
    n_pluses=$(echo ${pattern%.txt} | tr -cd "+" | wc -c)
    echo The ${n_pluses}\-group pattern $pattern belongs to reference group ${ref_groups[$n_pluses]}.
    pattern=${ref_groups[$n_pluses]}
    printf "\n"
  fi

  echo File ${file} becomes ${file%_*}_${pattern}.tmp. Add to Pattern_combined${NAME}.${POPULATION}_${pattern}.tmp  # remove shortest matching string from right (i.e. to right of *)
  test -f ${file} && echo cp ${file} ${file%_*}_${pattern}.tmp && printf "\n"
  test -f ${file} && cp ${file} ${file%_*}_${pattern}.tmp

  test -f ${file%_*}_${pattern}.tmp && echo cat ${file%_*}_${pattern}.tmp \>\> Pattern_combined${NAME}.${POPULATION}_${pattern}.tmp && printf "\n"
  test -f ${file%_*}_${pattern}.tmp && cat ${file%_*}_${pattern}.tmp >> Pattern_combined${NAME}.${POPULATION}_${pattern}.tmp

  test -f ${file%_*}_${pattern}.tmp && echo rm ${file%_*}_${pattern}.tmp && printf "\n"
  test -f ${file%_*}_${pattern}.tmp && rm ${file%_*}_${pattern}.tmp

  echo mv $file ${DIRECTORY}/${SUBDIR}/${file##*/} #extract file name and move to subdirectory
  mv $file ${DIRECTORY}/${SUBDIR}/${file##*/} #extract file name and move to subdirectory
  printf "\n"
done

# sum the first field for all unique occurences of fields 2 to 2 + $n_groups - 1
for file in Pattern_combined${NAME}.${POPULATION}_*.tmp
do
  PATTERN=${file##*_}
  PATTERN=${PATTERN%.tmp*}
  echo Pattern for file $file is $PATTERN.
  printf "\n"

  n_groups=$(($(echo $PATTERN | tr -cd "+" | wc -c)+1)) # total number of groups in pattern is one more than the number of "+'s"

  echo Combined file $file has $(awk '{total+=$1} END{print total}' $file) entries in $(awk 'END{print NR}' $file) lines.
  printf "\n"

  str=$(echo "awk 'BEGIN{OFS=\"\t\"} { array[sprintf(\"%s")
  for ((k=3; k<=$((n_groups+1)); k++));
  do
    str=${str}$(echo "\t%s")
  done
  str=${str}$(echo "\"")
  for ((k=2; k<=$((n_groups+1)); k++));
  do
    str=${str}$(echo ",\$$k")
  done
  str=${str}$(echo ")]+=\$1 } END{ for (row in array) { print array[row],row } }' ${DIRECTORY}/Pattern_combined${NAME}.${POPULATION}_${PATTERN}.tmp > ${DIRECTORY}/Pattern_combined${NAME}.${POPULATION}_${PATTERN}.unsorted.txt ")
  echo $str
  eval "$str"
  printf "\n"

  echo sort -k2 Pattern_combined${NAME}.${POPULATION}_${PATTERN}.unsorted.txt \> Pattern_combined${NAME}.${POPULATION}_${PATTERN}.txt
  sort -k2 Pattern_combined${NAME}.${POPULATION}_${PATTERN}.unsorted.txt > Pattern_combined${NAME}.${POPULATION}_${PATTERN}.txt
  printf "\n"

  echo Output file Pattern_combined${NAME}.${POPULATION}_${PATTERN}.txt has $(awk '{total+=$1} END{print total}' Pattern_combined${NAME}.${POPULATION}_${PATTERN}.txt) entries in $(awk 'END{print NR}' Pattern_combined${NAME}.${POPULATION}_${PATTERN}.txt) lines after sorting.
  printf "\n"
done

# remove combined file
for file in Pattern_combined${NAME}.${POPULATION}*unsorted.txt ];
do
  if [ -f $file ]
  then
    test -f $file && found+=1
    test -f $file && echo rm $file
    test -f $file && rm $file
  fi
done
test ! -z $found && (( $found>0 )) && found=0 && printf "\n"

#remove temporary unique file
for file in Pattern_combined${NAME}.${POPULATION}*tmp ];
do
  if [ -f $file ]
  then
    test -f $file && found+=1
    test -f $file && echo rm $file
    test -f $file && rm $file
  fi
done
test ! -z $found && (( $found>0 )) && found=0 && printf "\n"
