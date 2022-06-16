#! /bin/bash
set -e

# loops through all SIGMA-tuples of chromosomes in the columns of HAPLOTYPES_FILE and look for SNPs that are the same across the group.  Output the unique patterns and their counts for each group type.  Relies on subject_overlap.sh to do overlaps and comparisons

# HOME_DIR/ChromosomeOverlap_initiation_sub.sh haplotype_estimates_transpose.ukbb_bca_cases.chr11.69231642-69431642.1-100.txt 1 chr11.69231642-69431642 "" DIRECTORY HOME_DIR

HAPLOTYPES_FILE=$1 # transposed haplotypes file of snps x chromosomes
SIGMA=$2 # (One fewer than) the number of chromosomes to overlap
OUTPUT=$3 # Base file name for output Pattern.${OUTPUT}_*.txt
RANGE=$4 # (optional) range of the tuples to use, indexed as "from,to" or "step_size.step_number" starting from 0 and ending at (2n choose sigma + 1) - 1; if empty, computes all combinations
DIRECTORY=$5 # Folder to store output Pattern file
HOME_DIR=$6 # location of program files

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
echo Pattern files stored in $DIRECTORY
# Do NOT cd $DIRECTORY, because POPULATION files are in the submission directory, not the output directory

file_ID=${HAPLOTYPES_FILE%.*} # name of haplotype_estimates transpose file minus file extension

n_samples=$(awk 'BEGIN{nf=0} {if (NF>nf) {nf=NF} } END{print nf-1}' $HAPLOTYPES_FILE) # columns (does not include rsid column)
n_snps=$(awk 'END{print NR}' $HAPLOTYPES_FILE)

if [ -z $OUTPUT ]; # Optional file extension to name population, chromosome, and region
then
  OUTPUT=""
else
  OUTPUT=.${OUTPUT}
fi

if [ -z $SIGMA ]
then
  SIGMA=1 #default is to join 1 chromosomes to every chromosome in HAPLOTYPES_FILE, creating pairs
fi

num0=$(($n_samples-1)) # total number of haplotypes excluding rsid column in $HAPLOTYPES_FILE
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

if [ ! -z $RANGE ]; # user-provided comma-separated list lower,upper limits on combinations to perform in this round
then
  RANGE=($(echo $RANGE | perl -pne 's/([0-9]+)[.]+/$1 /g'))
  if ((${#RANGE[@]}>1)) # length=1 if TUPLE_RANGE is a comma-separated list
  then
    # For the case that RANGE is of the form STEP_SIZE.STEP
    ll=$(($((RANGE[1]-1))*${RANGE[0]})) # (i - 1) * step_size
    ul=$((RANGE[1]*RANGE[0]-1)) # i * step_size
  elif ((${#RANGE[@]}==1))
  then
    # For the case that RANGE is of the form LOWER,UPPER
    RANGE=($(echo $RANGE | perl -pne 's/([0-9]+)[,]+/$1 /g'))
    ll=${RANGE[0]} # starting tuple index
    ul=${RANGE[1]} # unding tuple index
  else
    (>&2 echo "Invalid range."; exit 1)
  fi
  (( $ll > $((n_tuples-1)) )) && ll=$((n_tuples-1))
  (( $ul > $((n_tuples-1)) )) && ul=$((n_tuples-1))
  (( $ll > $ul )) && (>&2 echo "Invalid range."; exit 1)
  RANGE=(${ll} ${ul})
else
  RANGE=(0 $((n_tuples-1))) # compute all combinations if none supplied
fi
echo Supplied range is \(${RANGE[0]},${RANGE[1]}\)
printf "\n"

echo sh ${HOME_DIR}index2combo2.sh ${RANGE[0]} $n_samples $SIGMA
out=($(sh ${HOME_DIR}index2combo2.sh ${RANGE[0]} $n_samples $SIGMA)) #get the tuple corresponding to combination ${RANGE[0]} in the list of combinations from 0 to (n_rows choose SIGMA) - 1
echo First joining tuple \(of rows\) is:
declare -p out
str="let k=${RANGE[0]}; for i1 in \`seq $((${out[0]}+1)) $n_samples\`;" # add 1 since 1st column is subject IDs, i.e., every column is shifted over by 1
for i in `seq 2 ${#out[@]}` # so i-1 is the current entry in "out"
do
  let j=$i-1 # index number of previous i
  str=${str}$(echo " do (( i$j > $((${out[$((j-1))]}+1)) )) && ll=\$((i$j+1)) || ll=$((${out[$((i-1))]}+1)); for i$i in \`seq \$ll $n_samples\`;") # if previous index i$j is not at its starting value, begin counting the current index i$((j+1)) from i$j+1; else use the initial value given in the out array
done
str=${str}$(echo " do INDEX[\$k]=\$i1") # increment INDEX key k and assign the value i1,i2,...,iSIGMA
for i in `seq 2 $SIGMA`
do
  str=${str}$(echo ",\$i$i")
done
let j=$n_samples+1
str=${str}$(echo ",$j; k=\$((k+1));")
str=${str}$(echo " if (( \$k=="$((${RANGE[1]}+1))" )); then break $SIGMA; fi;") # stop counting when reaching upper level of RANGE
for i in `seq 1 $SIGMA`
do
  str=${str}$(echo " done;")
done
echo $str
eval "$str" || echo error
printf "\n"

echo Array of joining tuples \(of rows\) is:
declare -p INDEX
n_combos=$((${#INDEX[@]}-1)) # number of combinations of subject of length SIGMA, corrected for 0-based indexing
printf "\n"

pattern_array=()
# loop through index combinations
for ((i=${RANGE[0]}; i<=${RANGE[1]}; i++))
do
  unset files_array
  error=1 #check that the only eror message is `$' treated as plain `$', else repeat call to subject_overlap.sh
  while (($error!=0))
  do
    echo Now running "sh ${HOME_DIR}ChromosomeOverlap_initiation.sh $HAPLOTYPES_FILE ${INDEX[$i]} 2> subject_overlap_error_${file_ID}_${INDEX[$i]}.err"
    str=$(echo "sh ${HOME_DIR}ChromosomeOverlap_initiation.sh $HAPLOTYPES_FILE ${INDEX[$i]} 2> subject_overlap_error_${file_ID}_${INDEX[$i]}.err")
    printf "\n"
    code=$(echo $str)
    start=$(date +%s.%N)
    eval "$code" || echo error
    finish=$(date +%s.%N)
    awk 'BEGIN{duration='$finish'-'$start';printf "Overlap completed in %0.2f s.\n",duration}'
    printf "\n"

    echo Check for errors:
    cat subject_overlap_error_${file_ID}_${INDEX[$i]}.err
    error=$(cat subject_overlap_error_${file_ID}_${INDEX[$i]}.err | grep "[^\`\$\' treated as plain \`\$']$" | wc -l)

    test -f subject_overlap_error_${file_ID}_${INDEX[$i]}.err && rm subject_overlap_error_${file_ID}_${INDEX[$i]}.err

    INDEX_SEP=($(echo ${INDEX[$i]} | perl -pne 's/,/ /g'))
    declare -p INDEX_SEP
    printf "\n"
    # Find the maximum index in the input; it will be replaced with "j" in combinations
    max=0;
    for j in ${INDEX_SEP[*]}
    do
      (( $j > $max )) && max=$j
    done
    INDEX_NEW=$(echo ${INDEX[$i]} | perl -pne "s/$max/j/g") # index with maximum number replaced by j
    echo $INDEX_NEW
    INDEX_ARRAY=($(echo ${INDEX[$i]} | perl -pne "s/,/ /g")) # remove commas in INDEX to create array
    printf "\n"
    # Get array of file names with indicated SIGMA-tuple

    # See if ref_pattern is successfully printed to array on the first iteration
    if (( i==${RANGE[0]} ));
    then
      echo Testing for reference pattern:
      files_array=($(for file in snps_unique_${file_ID}_${INDEX_NEW}*.txt
      do
        echo ${file%_*}
      done | sort -u))
      declare -p files_array # list of snps_unique files for each pattern type
      printf "\n"

      n_groups=$((${#files_array[@]}-1)) # number of types of bar groups, using 0-based indexing

      # get stamdard ordering in reference case from the 0-bar group
      for j in `seq 0 $n_groups`
      do
        n_commas=$(echo ${files_array[$j]##*_} | tr -cd "," | wc -c) # trim to left of last remaining underscore _
        echo Found $n_commas comma$( (( $n_commas!=1 )) && echo "s" || echo "") in file name ${files_array[$j]##*_}
        printf "\n"

        # find 0-bar group by maximum number of commas, equalling SIGMA
        if (( $n_commas==$SIGMA ));
        then
          ref_pattern=${files_array[$j]##*_}
          ref_pattern=($(echo $ref_pattern | perl -pne 's/([A-Za-z0-9]+)[,]*/$1 /g'))
          echo Pattern with no bars is the reference pattern:
          declare -p ref_pattern
          printf "\n"
        fi
      done

      echo Determine expected number of groups:
      echo sh ${HOME_DIR}StirlingSum.sh ${#INDEX_ARRAY[@]}
      S=$(sh ${HOME_DIR}StirlingSum.sh ${#INDEX_ARRAY[@]}) # total number of partitions of INDEX[$i]
      echo Expected $S file$((($S==1)) && echo "" || echo "s"), found ${#files_array[@]}
      printf "\n"
      if ((${#files_array[@]}!=$S))
      then
        echo error=\$\(\($error+1\)\)
        error=$(($error+1))
      fi
      echo Errors: $error

      if [ -z $ref_pattern ]
      then
        echo Reference pattern does not exist.
        error=$(($error+1))
      fi
    fi
    echo There $((($error==1)) && echo "was" || echo "were") $error error$((($error==1)) && echo "" || echo "s")
    # remove output files from last run of subject_overlap.sh if there was an error
    if (($error!=0))
    then
      for file in snps_unique_${file_ID}_${INDEX_NEW}_*.txt
      do
        test -f $file && found+=1
        test -f $file && echo rm $file
        test -f $file && rm $file
      done
      test ! -z $found && (( $found>0 )) && found=0 && printf "\n"
    else
      echo Proceeding to append files to ${ref_pattern[*]}
      printf "\n"
    fi
  done

  # Get files again after testing
  files_array=($(for file in snps_unique_${file_ID}_${INDEX_NEW}_*.txt
  do
    echo ${file%_*}
  done | sort -u))
  declare -p files_array
  printf "\n"

  n_groups=$((${#files_array[@]}-1)) # number of types of bar groups, using 0-based indexing

  # permutations to assign subsequent groupings to reference pattern type
  if (( i==${RANGE[0]} ));
  then
    # get stamdard ordering in reference case from the 0-bar group
    for ((j=0; j<=$n_groups; j++))
    do
      n_commas=$(($(echo ${files_array[$j]##*_} | tr -cd "," | wc -c)-1)) # trim to left of last remaining underscore _
      echo Found $n_commas comma$( (( $n_commas!=1 )) && echo "s" || echo "") in file name ${files_array[$j]##*_}
      printf "\n"
      # https://www.linuxjournal.com/content/bash-parameter-expansion
      # https://stackoverflow.com/questions/16679369/count-occurrences-of-a-char-in-a-string-using-bash
      # find 0-bar group by maximum number of commas, equalling SIGMA
      if (( $n_commas==$SIGMA ));
      then
        ref_pattern=${files_array[$j]##*_}
        ref_pattern=($(echo $ref_pattern | perl -pne 's/([A-Za-z0-9]+)[,]*/$1 /g'))
        echo Pattern with no bars is the reference pattern:
        declare -p ref_pattern
        printf "\n"
      fi
    done
  fi
  if (( i>=${RANGE[0]} ));
  then
    # get corresponding ordering in current case from the 0-bar group
    for ((j=0; j<=$n_groups; j++))
    do
      n_commas=$(echo ${files_array[$j]##*_} | tr -cd "," | wc -c)
      if (( $n_commas==$SIGMA ));
      then
        pattern=${files_array[$j]##*_}
        pattern=($(echo $pattern | perl -pne 's/([A-Za-z0-9]+)[,]*/$1 /g'))
        declare -p pattern
        declare -p ref_pattern
        printf "\n"
      fi
    done
  fi

  for ((j=0;j<=$n_groups;j++))
  do
    # Remove files of total snps found in each bar group (i.e., not unique)

    # Initiate a file "Pattern_" for each bar grouping if the overlap INDEX is first in the RANGE
    echo Prefix is ${files_array[$j]}
    if (( i==${RANGE[0]} ));
    then
      pattern_array+=(${files_array[$j]##*_})
      echo ${pattern_array[$j]}

      test -f Pattern${OUTPUT}_${pattern_array[$j]}.txt && echo rm Pattern${OUTPUT}_${pattern_array[$j]}.txt && printf "\n"
      test -f Pattern${OUTPUT}_${pattern_array[$j]}.txt && rm Pattern${OUTPUT}_${pattern_array[$j]}.txt && printf "\n"

      echo touch Pattern${OUTPUT}_${pattern_array[$j]}.txt
      touch Pattern${OUTPUT}_${pattern_array[$j]}.txt
      printf "\n"

      echo Pattern array:
      declare -p pattern_array
      printf "\n"
    fi
    printf "\n"
    # Join the unique SNP paetterns in each bar group in as many columns as there are one fewer than the number of bars
    for file in ${files_array[$j]}*.txt
    do
      test ! -f Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt && echo touch Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt && printf "\n"
      test ! -f Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt && touch Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt

      if [ -f $file ];
      then
        echo cat $file \> overlap${OUTPUT}_${pattern_array[$j]}.txt
        head $file
        cat $file > overlap${OUTPUT}_${pattern_array[$j]}.txt
        printf "\n"

        echo Array of files in the bar group ${pattern_array[$j]}:
        array_f=($(for file_f in ${files_array[$j]}*.txt; do echo $file_f; done | sort)) # array of overlaps for different bar groups in a single pattern, e.g., (...2 ...j)
        declare -p array_f
        printf "\n"

        # starting from $file, see which lines the other files have in common
        echo Get chromosome tuples in common among files of the bar group ${pattern_array[$j]}:
        for ((k=0;k<${#array_f[@]};k++));
        do
          # echo ${array_f[$k]}
          echo awk \-F\"\\t\" \'BEGIN{OFS=FS} NR==FNR{group[\$1]=\$1\;next} \$1 in group {print \$0}\' overlap${OUTPUT}_${pattern_array[$j]}.txt ${array_f[$k]} \> overlap${OUTPUT}_${pattern_array[$j]}_temp.txt
          awk -F"\t" 'BEGIN{OFS=FS} NR==FNR{group[$1]=$1;next} $1 in group {print $0}' overlap${OUTPUT}_${pattern_array[$j]}.txt ${array_f[$k]} > overlap${OUTPUT}_${pattern_array[$j]}_temp.txt # first field is the chromosome tuple, second type of sharing (e.g., 00,11), third is the list of shared SNPs
          echo head overlap${OUTPUT}_${pattern_array[$j]}_temp.txt
          head overlap${OUTPUT}_${pattern_array[$j]}_temp.txt
          printf "\n"

          test -f overlap${OUTPUT}_${pattern_array[$j]}_temp.txt && echo mv overlap${OUTPUT}_${pattern_array[$j]}_temp.txt overlap${OUTPUT}_${pattern_array[$j]}.txt && printf "\n"
          test -f overlap${OUTPUT}_${pattern_array[$j]}_temp.txt && mv overlap${OUTPUT}_${pattern_array[$j]}_temp.txt overlap${OUTPUT}_${pattern_array[$j]}.txt
        done
        # only pick lines from file shared by all other files in bar grouping; append column of found SNPs in bar group onto others in grouping
        echo Keep only chromosome tuples \(col 1\) from $file shared by all members of the bar group ${pattern_array[$j]}:
        echo awk \-F\"\\t\" \'BEGIN{OFS=FS} NR==FNR{group[\$1]=\$1\;next} \$1 in group {print \$0}\' overlap${OUTPUT}_${pattern_array[$j]}.txt $file \| cut \-f 3- \| paste Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt \- \> Pattern${OUTPUT}_${pattern_array[$j]}_temp1.txt
        awk -F"\t" 'BEGIN{OFS=FS} NR==FNR{group[$1]=$1;next} $1 in group {print $0}' overlap${OUTPUT}_${pattern_array[$j]}.txt $file | cut -f 3- | paste Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt - > Pattern${OUTPUT}_${pattern_array[$j]}_temp1.txt
        printf "\n"

        echo cat Pattern${OUTPUT}_${pattern_array[$j]}_temp1.txt \> Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt
        cat Pattern${OUTPUT}_${pattern_array[$j]}_temp1.txt > Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt
        echo head Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt
        head Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt
        printf "\n"

        test -f Pattern${OUTPUT}_${pattern_array[$j]}_temp1.txt && echo rm Pattern${OUTPUT}_${pattern_array[$j]}_temp1.txt && printf "\n"
        test -f Pattern${OUTPUT}_${pattern_array[$j]}_temp1.txt && rm Pattern${OUTPUT}_${pattern_array[$j]}_temp1.txt

        test -f overlap${OUTPUT}_${pattern_array[$j]}.txt && echo rm overlap${OUTPUT}_${pattern_array[$j]}.txt && printf "\n"
        test -f overlap${OUTPUT}_${pattern_array[$j]}.txt && rm overlap${OUTPUT}_${pattern_array[$j]}.txt

        unset array_f
      fi
    done
    # Append to the list of found SNP patterns for each bar grouping.  First replace elements of the current pattern with approriate elements from the reference pattern based on the pattern with no bars
    pattern_length=${#pattern[@]}
    let pattern_length=$pattern_length-1
    old=$(echo "(\b${pattern[0]}\b)(?![A-Za-z0-9])([,]*)")
    new=${ref_pattern[0]}$(echo "x\$2")
    # replace pattern integer, being not followed or trailed by integers, or by x, using negative lookahead ?! and word boundaries \b. The first captured group in the pattern integer, andt he second is the possibly trailing comma or space
    # https://www.regular-expressions.info/lookaround.html
    # https://www.rexegg.com/regex-boundaries.html
    str=$(echo "echo ${files_array[$j]##*_} | perl -pne 's/$old/$new/g'")
    # echo $str
    pattern_temp=$(eval "$str")
    for k in `seq 1 $pattern_length`
    do
      old=$(echo "(\b${pattern[$k]}\b)(?![A-Za-z0-9])([,]*)")
      new=${ref_pattern[$k]}$(echo "x\$2")
      str=$(echo "echo $pattern_temp | perl -pne 's/$old/$new/g'")
      echo $str
      eval $str
      pattern_temp=$(eval "$str")
      printf "\n"
    done
    pattern_temp=$(echo $pattern_temp | perl -pne 's/x//g')
    echo Append to Pattern${OUTPUT}_${pattern_temp}.txt
    echo cut \-f1 \--complement Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt \>\> Pattern${OUTPUT}_${pattern_temp}.txt
    cut -f1 --complement Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt >> Pattern${OUTPUT}_${pattern_temp}.txt
    # https://stackoverflow.com/questions/32812916/how-to-delete-the-first-column-which-is-in-fact-row-names-from-a-data-file-in
    printf "\n"

    test -f Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt && echo rm Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt && printf "\n"
    test -f Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt && rm Pattern${OUTPUT}_${pattern_array[$j]}_temp.txt

    unset pattern_temp
  done

  for file in snps*_${INDEX_NEW}_*.txt
  do
    test -f $file && found+=1
    test -f $file && echo rm $file
    test -f $file && rm $file
  done
  test ! -z $found && (( $found>0 )) && found=0 && printf "\n"
done

# replace pattern files with unique rows (columns 2,3,..) and their counts (1)
for file in Pattern${OUTPUT}*.txt
do
  echo Totaling unique lines of $file:
  nf=$(awk 'NR==1{ x=NF } NR>1{ if (NF > x) { x=NF } } END{ print x }' $file)
  printf "\n"
  if [ ! -z $nf ]
  then
    echo $nf field$( (($nf!=1)) && echo "s"|| echo "" ).
    printf "\n"
    # Arrange columns 2 to $nf in increasing order using asort() https://www.gnu.org/software/gawk/manual/html_node/Array-Sorting-Functions.html#Array-Sorting-Functions
    str=$(echo "awk -F\"\t\" 'BEGIN{ OFS=FS } { split(")
    for ((i=1; i<$nf; i++));
    do
      str=${str}$(echo "\$$i\"\t\"")
    done
    str=${str}$(echo "\$$nf,array,\"\t\"); asort(array); { print ")
    for ((i=1; i<$nf; i++));
    do
      str=${str}$(echo "array[$i],")
    done
    str=${str}$(echo "array[$nf] } }' $file > ${file%.txt}_temp.txt ")
    echo $str
    eval "$str"
    printf "\n"

    echo \"Slide over\" fields to the right into empty fields to the left: # https://www.unix.com/shell-programming-and-scripting/221003-how-detect-empty-field-awk.html
    echo awk \-F\"\\t\" \'BEGIN{ OFS=FS } { for \(i=1\; i\<NF\; i++\) { if \(\$i==\"\"\) { \$i=\$\(i+1\)\; \$\(i+1\)=\"\" } } } { print }\' ${file%.txt}_temp.txt \> $file
    awk -F"\t" 'BEGIN{ OFS=FS } { for (i=1; i<NF; i++) { if ($i=="") { $i=$(i+1); $(i+1)="" } } } { print }' ${file%.txt}_temp.txt > $file
    printf "\n"

    echo Sort lines, print counts of unique lines in new field on the left:
    echo cat $file \| sort \| uniq \-c \> ${file%.txt}_temp.txt
    cat $file | sort | uniq -c > ${file%.txt}_temp.txt

    let nf=$nf+1 # get maximum number of fields over all lines

    echo Format output as tab-delimited:
    str=$(echo "awk 'BEGIN{ OFS=\"\t\" } { print ")
    for ((i=1; i<$nf; i++));
    do
      str=${str}$(echo "\$$i,")
    done
    str=${str}$(echo "\$$nf }' ${file%.txt}_temp.txt > $file")
    echo $str
    eval "$str"
    printf "\n"

    # permute columns of file if there are multiple bar groups
    n_bar_groups=$(($(echo ${file##*_} | tr -c -d "+" | wc -c)+1)) #number of bar groups
    if (($n_bar_groups>1))
    then
      echo $n_bar_groups bar group$((($n_bar_groups==1)) && echo "" || echo "s") in $file
      printf "\n"
      # generate all permutations of the set {2,3,...,n_bar_groups+1}
      echo Generate permutations of the bar groups:
      str="2";
      for ((k=3; k<=$((n_bar_groups+1)); k++))
      do
        str=${str}",$k"
      done
      str0=$str # original, unpermuted list of columns
      str=$(perl -e "print '{$str},'x$n_bar_groups") # generate all $n_bar_groups-letter words on $n_bar_groups letters, with repeats
      str="for perm in ${str%,*}; do echo \$perm; done;"
      echo $str
      perms_array=($(eval $str | sort | uniq))
      declare -p perms_array
      printf "\n"

      echo Unpermuted file:
      file0=${file%_*}.unpermuted_${file##*_} # unpermuted file
      echo cat $file \> $file0
      cat $file > $file0
      printf "\n"

      for ((k=0; k<${#perms_array[@]}; k++))
      do
        # evaluate permuation if every element is unique
        if [ ${perms_array[k]} != $str0 ] # skip the un-permuted case
        then
          echo Determine if permutation ${perms_array[k]} is acceptable:
          str=$(echo "awk -F\"\t\" -v perm=${perms_array[k]} '{n=split(perm,cols,\",\"); idx=1; for (i=2;i<=n;i++) {for (j=1; j<i; j++) {if (cols[i]==cols[j]) {idx=0; break} } } } END{if (idx==1) {print \"1\"} else {print \"0\"}  }' $file0") # loop through permuted columns and check that they are distinct
          echo $str
          use_perm=$(eval "$str")
          printf "\n"

          if (($use_perm==1));
          then
            echo Permute columns:
            str=$(echo "awk -F\"\t\" -v perm=${perms_array[k]} 'BEGIN{OFS=FS} {n=split(perm,cols,\",\"); idx=1; for (i=2;i<=n;i++) {for (j=1; j<i; j++) {if (cols[i]==cols[j]) {idx=0; break} } } } {if (idx==1) {")
            for ((j=1; j<=$n_bar_groups; j++))
            do
              str=${str}$(echo " col$j=cols[$j];") # store permuted columns in array as variables
            done
            str=${str}$(echo " print \$1")
            for ((j=1; j<=$n_bar_groups; j++))
            do
              str=${str}$(echo ",\$col$j")
            done
            str=${str}$(echo " } }' $file0 > ${file%.txt}_temp.txt")
            echo $str
            eval "$str"
            printf "\n"

            echo Append permutation \(${perms_array[k]}\) to $file
            echo cat ${file%.txt}_temp.txt \>\> $file
            cat ${file%.txt}_temp.txt >> $file
            printf "\n"
          fi
          if [ -f ${file%.txt}_temp.txt ];
          then
            echo rm ${file%.txt}_temp.txt
            rm ${file%.txt}_temp.txt
          fi
        fi
      done
    fi
  fi
  # move output to output directory if not already in it
  if [ "$PWD" != "$DIRECTORY" ];
  then
    echo Moving list of patterns to output directory:
    test -f $file && echo mv $file ${DIRECTORY}/${file} && printf "\n"
    test -f $file && mv $file ${DIRECTORY}/${file}
  fi

  test -f ${file%.txt}_temp.txt && echo rm ${file%.txt}_temp.txt && printf "\n"
  test -f ${file%.txt}_temp.txt && rm ${file%.txt}_temp.txt
done

test -f samples${OUTPUT}.sigma_${SIGMA}.txt && echo rm samples${OUTPUT}.sigma_${SIGMA}.txt && printf "\n"
test -f samples${OUTPUT}.sigma_${SIGMA}.txt && rm samples${OUTPUT}.sigma_${SIGMA}.txt
