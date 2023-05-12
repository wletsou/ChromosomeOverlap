#! /bin/bash
set -e

# Forms all SIGMA-tuples of the rows of FILE and looks for common alleles.

HAPLOTYPE=$1 # haplotype to be included in every overlap
FILE=$2 # Pattern_combined file with first field counts and second and further fields patterns
TRANSPOSE_FILE=$3 # optional transposed version of FILE; searched for and/or created if empty
RANGE=$4 # Range (comma-separated list, lower,upper) of intersections (of one group of SIGMA rows with all others) to do, from 0 to (n_rows choose sigma) - 1.  OR a period-separated list "step_size.step_number" defining the range "(i-1) * step,i * step - 1"
STEP_SIZE=$5 # number of overlaps to do in one step
SIGMA=$6 # number of rows to overlap
NAME=$7 # optional tag for output files Overlap_tuples.NAME..., by iteration
DIRECTORY=$8 # output directory to store Output.Tuples.Iteration00x.Range files
HOME_DIR=$9 # location of program files

if [ -z $HOME_DIR ]
then
  unset HOME_DIR
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
echo Overlap files stored in $DIRECTORY
printf "\n"

module load R/3.6.1

if [ -z $SIGMA ];
then
  SIGMA=1 # number of rows to overlap with a given row
fi

[ -z $FILE ] && (>&2 echo "Pattern_combined file not supplied."; exit 1)
[ ! -f $FILE ] && (>&2 echo "Pattern_combined file not found."; exit 1)

if [ -z $NAME ]
then
  unset NAME
fi

n_rows=$(awk 'END{print NR}' $FILE) # number of patterns
num0=$n_rows
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
(($n_rows==0)) && n_tuples=1 # Correction for one way of arranging 0 things
echo $n_tuples total combination$((($n_tuples==1)) || echo "s").

if [ ! -z $RANGE ]; # user-provided comma-separated list lower,upper limits on combinations to perform in this round
then
  RANGE=($(echo $RANGE | perl -pne 's/([0-9]+)[.]+/$1 /g'))
  if ((${#RANGE[@]}>1)) # length=1 if TUPLE_RANGE is a comma-separated list
  then
    # For the case that RANGE is of the form STEP_SIZE.STEP
    ll=$(($((RANGE[1]-1))*${RANGE[0]})) # (i - 1) * step_size
    ul=$((RANGE[1]*RANGE[0]-1)) # i * step_size
  elif ((${#RANGE[@]}==2))
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

((${RANGE[0]}>$((n_tuples-1)))) && (>&2 echo "Range exceeds total number of tuples."; exit 1)
((${RANGE[1]}>$((n_tuples-1)))) && RANGE[1]=$((n_tuples-1))
char=$(( $( echo $( (( $n_tuples>${RANGE[1]} )) && echo $n_tuples || echo ${RANGE[1]} ) | wc -m)-1)) # number of characters in the number of tuples
n_fields=$(awk 'BEGIN{x=0} (NF>x){x=NF} END{print x}' $FILE)

if [ -z $STEP_SIZE ]
then
  STEP_SIZE=$((${RANGE[1]}-${RANGE[0]}+1))
fi

if [ -z $TRANSPOSE_FILE ]
then
  TRANSPOSE_FILE=${FILE%.*}.transpose.${FILE##*.}
  exec 3<>$TRANSPOSE_FILE
  # (
  if flock -x -n 3
  then
    echo Obtained lock.
    echo wc \-l $TRANSPOSE_FILE
    if [ ! -s $TRANSPOSE_FILE ]
    then
      echo Transpose input file:
      echo awk \'BEGIN{OFS=\"\\t\"}\; {for\(j=1\;j\<=NF\;j++\) {a[NR,j]=\$j\; n_rows=NR\; n_cols=\(n_cols\<NF?NF:n_cols\)} } END{for \(j=1\;j\<=n_cols\;j++\) {for \(i=1\;i\<=n_rows\;i++\) {printf \"%s%s\",a[i,j],\(i==n_rows?\"\\n\":\"\\t\"\)} } }\' $FILE \> $TRANSPOSE_FILE # transpose to counts in first row, snp patterns by bar groups in subsequent rows
      awk 'BEGIN{OFS="\t"}; {for(j=1;j<=NF;j++) {a[NR,j]=$j; n_rows=NR; n_cols=(n_cols<NF?NF:n_cols)} } END{for (j=1;j<=n_cols;j++) {for (i=1;i<=n_rows;i++) {printf "%s%s",a[i,j],(i==n_rows?"\n":"\t")} } }' $FILE > $TRANSPOSE_FILE && printf "\n"
    fi
  else
    while ! flock -x -n 3
    do
      echo Did not obtain lock.
      sleep 0.1 && printf "\n"
    done
  fi
fi

echo wc \-l $TRANSPOSE_FILE
wc -l $TRANSPOSE_FILE && printf "\n"
echo awk \'{print NF}\' $TRANSPOSE_FILE
awk '{print NF}' $TRANSPOSE_FILE && printf "\n"

file_ID=${FILE%.*} # drop file extension
file_ID=${file_ID##*/} # delete to the left of the last /
echo File identifier is $file_ID
printf "\n"

test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && echo rm Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && printf "\n"
test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && rm Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
echo touch Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
touch Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt

if (($n_rows > 0))
then
  for ((k=1;k<=$(($((${RANGE[1]}-${RANGE[0]}+1+STEP_SIZE-1))/STEP_SIZE));k++)) # https://stackoverflow.com/questions/2395284/round-a-divided-number-in-bash
  do
    range[0]=$((${RANGE[0]}+(k-1)*STEP_SIZE))
    ((${range[0]}>${RANGE[1]})) && range[0]=${RANGE[1]}
    range[1]=$((${RANGE[0]}+k*STEP_SIZE-1))
    ((${range[1]}>${RANGE[1]})) && range[1]=${RANGE[1]}
    declare -p range
    echo Rscript ${HOME_DIR/%/\/}index2combo2.R I=${range[0]} n=$n_rows sigma=$SIGMA
    Rscript ${HOME_DIR/%/\/}index2combo2.R I=${range[0]} n=$n_rows sigma=$SIGMA
    out=($(Rscript ${HOME_DIR/%/\/}index2combo2.R I=${range[0]} n=$n_rows sigma=$SIGMA)) # get the tuple corresponding to combination ${RANGE[0]} in the list of combinations from 0 to (n_rows choose SIGMA) - 1
    printf "\n"

    echo First joining tuple \(of rows\) is:
    declare -p out
    printf "\n"

    str=$(echo "awk 'BEGIN{OFS=\"\\t\"; n0=split(\"$HAPLOTYPE\",array,\",\"); for (j=1;j<=n0;j++) {snps_0[array[j]]}; delete array; i0=1} NR>1{ k=${range[0]}; ")
    for ((i=1;i<${#out[@]};i++))
    do
      str=${str}$(echo "if (i$((i-1))>${out[$((i-2))]}) {ll=i$((i-1))+1} else {ll=${out[$((i-1))]}}; for (i$i=ll;i$i<=$n_rows;i$i++) {n$i=split(\$i$i,array,\",\"); delete snps_$i; for (j=1;j<=n$i;j++) {snps_$i[array[j]]}; delete array; ")
    done
    str=${str}$(echo "if (i$((${#out[@]}-1))>${out[$((${#out[@]}-2))]}) {ll=i$((${#out[@]}-1))+1} else {ll=${out[$((${#out[@]}-1))]}}; for (i${#out[@]}=ll;i${#out[@]}<=$n_rows;i${#out[@]}++) { n${#out[@]}=split(\$i${#out[@]},snps_${#out[@]},\",\"); idx=0; printf \"%s\t\", 1; for (j=1;j<=n${#out[@]};j++) { if (")
    for ((i=1;i<${#out[@]};i++))
    do
      if ((i>1))
      then
        str=$str$(echo " && ")
      fi
      str=$str$(echo "snps_${#out[@]}[j] in snps_$i")
    done
    if [ ! -z $HAPLOTYPE ] && (( $(awk 'BEGIN{n0=split("'$HAPLOTYPE'",array,","); print n0}')>0 )) # if a non-zero-length haplotype was supplied, then add it to the comparrison
    then
      if ((i>1))
      then
        str=$str$(echo " && ")
      fi
      str=$str$(echo "snps_${#out[@]}[j] in snps_0")
    fi
    str=${str}$(echo ") {printf \"%s%s\",(idx==0?\"\":\",\"),snps_${#out[@]}[j]; idx+=1 } ")
    for ((i=${#out[@]};i>1;i--))
    do
      str=${str}$(echo "} ")
    done
    if ((${#out[@]}==1))
    then
      str=${str}$(echo "} ")
    fi
    str=${str}$(echo "if (idx==0) {printf \"0\"}; k+=1; printf \"\n\"; if (k>${range[1]}) {next} ")
    for ((i=${#out[@]};i>=1;i--))
    do
      str=${str}$(echo "} ")
    done
    str=${str}$(echo "}' $TRANSPOSE_FILE >> Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt")
    echo $str
    printf "\n"
    start=$(date +%s.%N)
    eval $str
    finish=$(date +%s.%N)
    duration=$(awk 'BEGIN{printf "%0.2f\n", '$finish'-'$start'}')
    echo Overlap took $duration second$( [ "$duration" == "1.00" ] && echo "" || echo "s").
    printf "\n"
  done
fi

if (($n_fields>2))
then
  # Make a tab-delimined list of columns 2,3,...,n_fields and split into an array. Sort the array and print
  echo Sort fields of Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt:
  str=$(echo "awk 'BEGIN{ OFS=\"\t\" } { split(")
  for ((i=2; i<$n_fields; i++));
  do
    str=${str}$(echo "\$$i\"\t\"")
  done
  str=${str}$(echo "\$$n_fields,array,\"\t\"); asort(array); { print \$1,")
  for ((i=1; i<$(($n_fields-1)); i++));
  do
    str=${str}$(echo "array[$i],")
  done
  str=${str}$(echo "array[$(($n_fields-1))] } }' Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt > Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp ")
  echo $str
  start=$(date +%s.%N)
  printf "\n"
  output=$(eval "$str" 2> /dev/null) # eliminate error message
  finish=$(date +%s.%N)
  duration=$(awk 'BEGIN{printf "%0.2f\n", '$finish'-'$start'}')
  echo Sorting fields took $duration second$( [ "$duration" == "1.00" ] && echo "" || echo "s").
  if (( $(echo $output | wc -w) != 0 )) # check if output is 0-length (i.e., an error)
  then
    echo cat Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt > Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp
    #cat Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt > Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp
  fi
  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && echo mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && printf "\n"
  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
fi

while ! test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
do
  echo Waiting for Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt...
  sleep 10
done

# ! test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && (&>2 echo "Could not find Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt"; exit 1)

# Sort SNPs within a field
sort_snps=0
if (($sort_snps==1))
then
  echo Sort SNPs within each group:
  start=$(date +%s.%N)
  str=$(echo "awk 'BEGIN{ OFS=\"\t\" } { split(\$2,array2,\",\"); n=asort(array2); idx=0; for (j=1;j<=n;j++) {if (idx==0) {\$2=array2[j]; idx++} else {\$2=\$2\",\"array2[j]} } }")
  for ((i=3; i<=$n_fields; i++));
  do
    str=${str}$(echo " { split(\$$i,array$i,\",\"); n=asort(array$i); idx=0; for (j=1;j<=n;j++) {if (idx==0) {\$$i=array$i[j]; idx++} else {\$$i=\$$i\",\"array$i[j]} } }")
  done
  str=${str}$(echo " { print }' Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt > Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp")
  echo $str
  eval $str
  finish=$(date +%s.%N)
  duration=$(awk 'BEGIN{printf "%0.2f\n", '$finish'-'$start'}')
  echo Sorting SNPs took $duration second$( [ "$duration" == "1.00" ] && echo "" || echo "s").
  printf "\n"
  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && echo mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && printf "\n"
  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
fi

while ! test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
do
  echo Waiting for Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt...
  sleep 10
done

if (($n_fields>2))
then
  start=$(date +%s.%N)
  echo \"Slide over\" fields to the right into empty fields to the left:
  echo awk \'BEGIN{ OFS=\"\\t\" } { for \(i=1\; i\<NF\; i++\) { if \(\$i==\"\"\) { \$i=\$\(i+1\)\; \$\(i+1\)=\"\" } } } { print }\' Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt \> Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp
  awk 'BEGIN{ OFS="\t" } { for (i=1; i<NF; i++) { if ($i=="") { $i=$(i+1); $(i+1)="" } } } { print }' Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt > Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp
  printf "\n"
  finish=$(date +%s.%N)
  duration=$(awk 'BEGIN{printf "%0.2f\n", '$finish'-'$start'}')
  echo Sliding over took $duration second$( [ "$duration" == "1.00" ] && echo "" || echo "s").
  printf "\n"
  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && echo mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && printf "\n"
fi

while ! test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
do
  echo Waiting for Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt...
  sleep 10
done

# Sum the first column of rows with matching patterns in all of columns 2,3,...,n_fields
echo Sum counts of matching rows:
start=$(date +%s.%N)
str=$(echo "awk 'BEGIN{OFS=\"\t\"} { array[")
for ((i=2; i<$n_fields; i++));
do
  str=${str}$(echo "\$$i\"\t\"")
done
str=${str}$(echo "\$$n_fields]+=\$1 } END{ for (row in array) { print array[row],row } }' Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt > Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp")
echo $str
eval "$str"
finish=$(date +%s.%N)
duration=$(awk 'BEGIN{printf "%0.2f\n", '$finish'-'$start'}')
echo Summing took $duration second$( [ "$duration" == "1.00" ] && echo "" || echo "s").
printf "\n"

while ! test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp
do
  echo Waiting for Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp...
  sleep 10
done

test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && echo mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && printf "\n"
test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt

while ! test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
do
  echo Waiting for Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt...
  sleep 10
done

# eliminate null rows
check_zeros=1
if (($check_zeros==1))
then
  echo Eliminate rows with all 0\'s:
  start=$(date +%s.%N)
  str=$(echo "awk 'BEGIN{OFS=\"\t\"} { if (")
  for ((i=2; i<$n_fields; i++))
  do
    str=${str}$(echo "\$$i != 0 &&")
  done
  str=${str}$(echo "\$$n_fields != 0) { print \$0 } }' Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt > Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp")
  echo $str
  eval "$str"
  finish=$(date +%s.%N)
  duration=$(awk 'BEGIN{printf "%0.2f\n", '$finish'-'$start'}')
  echo Elimination took $duration second$( [ "$duration" == "1.00" ] && echo "" || echo "s").
  printf "\n"

  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && echo mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && printf "\n"
  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
fi

while ! test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
do
  echo Waiting for Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt...
  sleep 10
done

# Move to output directory
f_name=$(printf "Overlap_tuples${NAME/#/.}.%0.${char}d-%0.${char}d.%s.txt" ${RANGE[0]} ${RANGE[1]} ${file_ID})
# check if error during mv operation
if [[ "${PWD/%/\/}Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt" != "${DIRECTORY/%/\/}${f_name}" ]]
then
  test -f  Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && echo mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt ${DIRECTORY/%/\/}${f_name} && printf "\n"
  test -f  Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt ${DIRECTORY/%/\/}${f_name}
fi

while ! test -f ${DIRECTORY/%/\/}${f_name}
do
  echo Waiting for ${DIRECTORY/%/\/}${f_name}...
  sleep 10
done
