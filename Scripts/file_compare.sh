#! /bin/bash

# Look for lines in TEST_FILE (from column 2 over) not in REF_FILE and write to OUT.INDEX, where INDEX indicates the comparison job index when dividing the comparison job into groups of DELTA lines in the TEST_FILE

TEST_FILE=$1 # file to be compared
REF_FILE=$2 # reference for comparison
INDEX=$3 # comparison group
DELTA=$4 # number of lines per comparison job
OUT=$5 # output prefix for file of unique lines
N_FIELDS=$6 # one more than the number of groups in the pattern, less than or equal to the maximum number of fields in the reference file

if [ -z $INDEX ];
then
  INDEX="1"
fi

if [ -z $OUT ];
then
  OUT="comparison"
fi

N_LINES=$(awk 'END{print NR}' $TEST_FILE) # total number of comparisons to be divided up into blocks of size $DELTA

N_JOBS=$(echo $DELTA | awk -v x=$N_LINES '{print x/$1}') # exact number of comparison jobs
N_JOBS=$(echo $N_JOBS | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }') # round up number of comparison jobs

if [ -z $N_FIELDS ];
then
  N_FIELDS=$(awk 'BEGIN{x=2} NR>=1{if (NF > x) x=NF} END{print x}' $REF_FILE) # maximum number of fields in pattern, including first column of counts
fi

char=$(($(echo $N_JOBS | wc -m)-1)) # the number of characters in the total number of jobs

ll=$((DELTA*(INDEX-1)+1))
ul=$((DELTA*INDEX))
str=$(echo "awk -F\"\t\" 'FILENAME==ARGV[1]{lines[\$2") # print lines in TEST_FILE that are not in REF_FILE
for ((j=3;j<=$N_FIELDS;j++))
do
  str=${str}$(echo "\"\t\"\$$j")
done
str=${str}$(echo "]; next} (!(\$2")
for ((j=3;j<=$N_FIELDS;j++))
do
  str=${str}$(echo "\"\t\"\$$j")
done
str=${str}$(echo " in lines) && FNR>=$ll && FNR<=$ul) {print \$2")
for ((j=3;j<=$N_FIELDS;j++))
do
  str=${str}$(echo "\"\t\"\$$j")
done
str=${str}$(echo "}' $REF_FILE $TEST_FILE > "$(printf " $OUT.%0.${char}d.txt" ${INDEX})"")
#str="awk 'NR==FNR{lines[\$2]=\$2; next} (!(\$2 in lines) && FNR>=$ll && FNR<$ul){print \$2}' $REF_FILE $TEST_FILE > "$(printf " $OUT.%0.${char}d.txt" ${INDEX})""
echo $str
eval "$str"
