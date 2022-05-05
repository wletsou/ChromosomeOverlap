#! /bin/bash
set -e

POPULATION=$1 # comma-separated list of two haplotypes files (subjectx x rsid with header and row names)
PATTERNS=$2 # two-column file of pattern count and patterns in the form column1_value1,column2_value2 in the first population
INDEX=$3 # optional 1-based index in the form step_size.step_number
DIRECTORY=$4
HOME_DIR=$5

module load R/3.6.1
module load python/3.7.0
module load python/conda/3.8.1
module load conda3/5.1.0

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ]
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY
echo Working directory is $DIRECTORY
printf "\n"

POPULATION=($(echo $POPULATION | perl -pne 's/[,]/ /g'))

echo Create combined haplotypes file:
test -f haplotype_estimates_combined.txt && echo rm haplotype_estimates_combined.txt
test -f haplotype_estimates_combined.txt && rm haplotype_estimates_combined.txt
if [ ! -f haplotype_estimates_combined.txt ]
then
  echo touch haplotype_estimates_combined.txt
  touch haplotype_estimates_combined.txt
  for ((j=0;j<${#POPULATION[@]};j++))
  do
    echo awk \'BEGIN{OFS=\"\\t\"} NR\>1{\$1=$((${#POPULATION[@]}-1-j))\; print \$0}\' ${POPULATION[j]} \>\> haplotype_estimates_combined.txt
    awk 'BEGIN{OFS="\t"} NR>1{$1='$((${#POPULATION[@]}-1-j))'; print $0}' ${POPULATION[j]} >> haplotype_estimates_combined.txt # combined cases and controls haplotypes without header, 1st field is case(1)/control(0) status
  done
fi
printf "\n"

# n_patterns=$(cat $PATTERNS | wc -l)
# echo $n_patterns
echo awk \'{print \$2}\' $PATTERNS
patterns=($(awk '{print $2}' $PATTERNS))
n_patterns=${#patterns[@]}
echo Total number of unique patterns is $n_patterns.
printf "\n"

if [ ! -z $INDEX ]
then
  INDEX="${INDEX}."
  INDEX_ARRAY=($(echo $INDEX | perl -pne 's/([0-9]+)[.]*/$1 /g')) # array (STEP_SIZE STEP_NO)
else
  INDEX=""
  INDEX_ARRAY=($n_patterns 1)
fi
ll=$(((${INDEX_ARRAY[1]}-1)*${INDEX_ARRAY[0]}+1))
ul=$((${INDEX_ARRAY[1]}*${INDEX_ARRAY[0]}))
(( $ll<1 )) && ll=1
(( $ul>$n_patterns )) && ul=$n_patterns

n_jobs=$(awk 'BEGIN{printf "%0.25f\n",'$n_patterns'/'${INDEX_ARRAY[0]}'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }') # number of jobs needed to complete all patterns in steps of size ${INDEX_ARRAY[0]}
char=$(($(echo $n_jobs | wc -c)-1)) # how big to make "step" field in filename

filename=${PATTERNS/combined/differences}
filename=$(printf "${filename/%.txt}_step%0.${char}d.txt" ${INDEX_ARRAY[1]})

test -f $filename && echo rm $filename && printf "\n"
test -f $filename && rm $filename
test ! -f $filename && echo touch $filename && printf "\n"
test ! -f $filename && touch $filename

echo Range: \(${ll},${ul}\)
for ((i=$ll;i<=$ul;i++))
do
  start=$(date +%s.%N) # start time in seconds.nanoseconds
  echo Evaluate pattern $i:
  # echo awk \'BEGIN{OFS=\"\\t\"} \(NR==FNR \&\& NR==$i\){ pattern=\$2\; m=split\(\$2,pattern_array,\",\"\)\; for \(i=1\;i\<=m\;i++\) {b=gensub\(\"0*\([1-9][0-9]*\)_\([0-9]*\)\",\"\\\\1\",\"g\",pattern_array[i]\)\; c=gensub\(\"0*\([1-9][0-9]*\)_\([0-9]*\)\",\"\\\\2\",\"g\",pattern_array[i]\)\; array[b+1]=c}\; next} \(FILENAME==ARGV[2]\){ n_chromosomes[\$1]+=1\; count=1\; for \(i=2\;i\<=\(m+1\)\;i++\) {if \(\$i!=array[i]\) {count=0\; break} }\; total[\$1]+=count\; } END{mu[1]=total[1] / n_chromosomes[1]\; mu[0]=total[0] / n_chromosomes[0]\; sigma2[1]=\(total[1] \* \(1 \- mu[1]\) ^ 2 + \(n_chromosomes[1] \- total[1]\) \* mu[1] ^ 2\) / \(n_chromosomes[1] \- 1\)\; sigma2[0]=\(total[0] \* \(1 \- mu[0]\) ^ 2 + \(n_chromosomes[0] \- total[0]\) \* mu[0] ^ 2\) / \(n_chromosomes[0] \- 1\)\; se[1]=sqrt\(sigma2[1]/n_chromosomes[1]\)\; se[0]=sqrt\(sigma2[0]/n_chromosomes[0]\)\; print pattern,mu[1],se[1],n_chromosomes[1],mu[0],se[0],n_chromosomes[0]}\' $PATTERNS haplotype_estimates_combined.txt \>\> $filename
  # awk 'BEGIN{OFS="\t"} (NR==FNR && NR=='$i'){ pattern=$2; m=split($2,pattern_array,","); for (i=1;i<=m;i++) {b=gensub("0*([1-9][0-9]*)_([0-9]*)","\\1","g",pattern_array[i]); c=gensub("0*([1-9][0-9]*)[_]([0-9]*)","\\2","g",pattern_array[i]); array[b+1]=c}; next} (FILENAME==ARGV[2]){ n_chromosomes[$1]+=1; count=1; for (i=2;i<=(m+1);i++) {if ($i!=array[i]) {count=0; break} }; total[$1]+=count; } END{mu[1]=total[1] / n_chromosomes[1]; mu[0]=total[0] / n_chromosomes[0]; sigma2[1]=(total[1] * (1 - mu[1]) ^ 2 + (n_chromosomes[1] - total[1]) * mu[1] ^ 2) / (n_chromosomes[1] - 1); sigma2[0]=(total[0] * (1 - mu[0]) ^ 2 + (n_chromosomes[0] - total[0]) * mu[0] ^ 2) / (n_chromosomes[0] - 1); se[1]=sqrt(sigma2[1]/n_chromosomes[1]); se[0]=sqrt(sigma2[0]/n_chromosomes[0]); print pattern,mu[1],se[1],n_chromosomes[1],mu[0],se[0],n_chromosomes[0]}' $PATTERNS haplotype_estimates_combined.txt >> $filename
  # finish=$(date +%s.%N)
  # awk 'BEGIN{duration='$finish'-'$start';printf "Completed in %0.2f s.\n",duration}'
  # printf "\n"

  # Data for t-test
  # echo awk \'BEGIN{OFS=\"\\t\"\; pattern=\"${patterns[$i-1]}\"\; m=split\(pattern,pattern_array,\",\"\)\; for \(i=1\;i\<=m\;i++\) {b=gensub\(\"0*\([1-9][0-9]*\)_\([0-9]*\)\",\"\\\\1\",\"g\",pattern_array[i]\)\; c=gensub\(\"0*\([1-9][0-9]*\)_\([0-9]*\)\",\"\\\\2\",\"g\",pattern_array[i]\)\; array[b+1]=c} } { n_chromosomes[\$1]+=1\; count=1\; for \(i=2\;i\<=\(m+1\)\;i++\) {if \(\$i!=array[i]\) {count=0\; break} }\; total[\$1]+=count\; } END{mu[1]=total[1] / n_chromosomes[1]\; mu[0]=total[0] / n_chromosomes[0]\; sigma2[1]=\(total[1] \* \(1 \- mu[1]\) ^ 2 + \(n_chromosomes[1] \- total[1]\) \* mu[1] ^ 2\) / \(n_chromosomes[1] \- 1\)\; sigma2[0]=\(total[0] \* \(1 \- mu[0]\) ^ 2 + \(n_chromosomes[0] \- total[0]\) \* mu[0] ^ 2\) / \(n_chromosomes[0] \- 1\)\; se[1]=sqrt\(sigma2[1]/n_chromosomes[1]\)\; se[0]=sqrt\(sigma2[0]/n_chromosomes[0]\)\; print pattern,mu[1],se[1],n_chromosomes[1],mu[0],se[0],n_chromosomes[0]}\' haplotype_estimates_combined.txt \>\> $filename
  # awk 'BEGIN{OFS="\t"; pattern="'${patterns[$i-1]}'"; m=split(pattern,pattern_array,","); for (i=1;i<=m;i++) {b=gensub("0*([1-9][0-9]*)_([0-9]*)","\\1","g",pattern_array[i]); c=gensub("0*([1-9][0-9]*)_([0-9]*)","\\2","g",pattern_array[i]); array[b+1]=c} } { n_chromosomes[$1]+=1; count=1; for (i=2;i<=(m+1);i++) {if ($i!=array[i]) {count=0; break} }; total[$1]+=count; } END{mu[1]=total[1] / n_chromosomes[1]; mu[0]=total[0] / n_chromosomes[0]; sigma2[1]=(total[1] * (1 - mu[1]) ^ 2 + (n_chromosomes[1] - total[1]) * mu[1] ^ 2) / (n_chromosomes[1] - 1); sigma2[0]=(total[0] * (1 - mu[0]) ^ 2 + (n_chromosomes[0] - total[0]) * mu[0] ^ 2) / (n_chromosomes[0] - 1); se[1]=sqrt(sigma2[1]/n_chromosomes[1]); se[0]=sqrt(sigma2[0]/n_chromosomes[0]); print pattern,mu[1],se[1],n_chromosomes[1],mu[0],se[0],n_chromosomes[0]}' haplotype_estimates_combined.txt >> $filename
  # finish=$(date +%s.%N)
  # awk 'BEGIN{duration='$finish'-'$start';printf "Completed in %0.2f s.\n",duration}'
  # printf "\n"

  # Data for Fisher's exact test
  echo awk \'BEGIN{OFS=\"\\t\"\; pattern=\"${patterns[$i-1]}\"\; m=split\(pattern,pattern_array,\",\"\)\; for \(i=1\;i\<=m\;i++\) {b=gensub\(\"0*\([1-9][0-9]*\)_\([0-9]*\)\",\"\\\\1\",\"g\",pattern_array[i]\)\; c=gensub\(\"0*\([1-9][0-9]*\)_\([0-9]*\)\",\"\\\\2\",\"g\",pattern_array[i]\)\; array[b+1]=c} } { n_chromosomes[\$1]+=1\; count=1\; for \(i in array\) {if \(\$i!=array[i]\) {count=0\; break} }\; total[\$1]+=count\; } END{a=total[1]\; b=n_chromosomes[1]\-a\; c=total[0]\; d=n_chromosomes[0]\-c\; print pattern,a,b,c,d}\' haplotype_estimates_combined.txt \>\> $filename
  awk 'BEGIN{OFS="\t"; pattern="'${patterns[$i-1]}'"; m=split(pattern,pattern_array,","); for (i=1;i<=m;i++) {b=gensub("0*([1-9][0-9]*)_([0-9]*)","\\1","g",pattern_array[i]); c=gensub("0*([1-9][0-9]*)_([0-9]*)","\\2","g",pattern_array[i]); array[b+1]=c} } { n_chromosomes[$1]+=1; count=1; for (i in array) {if ($i!=array[i]) {count=0; break} }; total[$1]+=count; } END{a=total[1]; b=n_chromosomes[1]-a; c=total[0]; d=n_chromosomes[0]-c; print pattern,a,b,c,d}' haplotype_estimates_combined.txt >> $filename
  finish=$(date +%s.%N)
  awk 'BEGIN{duration='$finish'-'$start';printf "Completed in %0.2f s.\n",duration}'
  printf "\n"
done

# test -f haplotype_estimates_combined.txt && echo rm haplotype_estimates_combined.txt && printf "\n"
# test -f haplotype_estimates_combined.txt && rm haplotype_estimates_combined.txt

# echo Get 1-sided p_values for pattern differences:
# echo Rscript ${HOME_DIR}/ukbb_pattern_difference.R file=$filename
# Rscript ${HOME_DIR}/ukbb_pattern_difference.R file=$filename

echo Get Fisher\'s exact test p_values for pattern partitioning:
echo python ${HOME_DIR}/ukbb_fisher_exact.py -f $filename -d $DIRECTORY
python ${HOME_DIR}/ukbb_fisher_exact.py -f $filename -d $DIRECTORY
