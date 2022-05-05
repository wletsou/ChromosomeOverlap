#! /bin/bash
set -e

# sh /home/wletsou/scripts/ukbb_haplotype_partition.sh haplotype_estimates.ukbb_bca_20200116_cases.chr17.5222961-5436196.txt,haplotype_estimates.ukbb_bca_20200116_controls.chr17.5222961-5436196.txt Core_patterns_names.ukbb_bca_20200116_cases.subset_9346.chr17.5222961-5436196_2,3,j.txt 17 5222961,5436196 5.1

# sh /home/wletsou/scripts/ukbb_haplotype_partition.sh haplotype_estimates.ukbb_bca_20200116_cases.chr11.68881416-69455872.txt,haplotype_estimates.ukbb_bca_20200116_controls.chr11.68881416-69455872.txt "rs2046494_C=0,rs4980661_A=0,rs657315_T=0,rs637185_T=0,rs498931_G=0,rs79241527_T=0,rs1122316_A=0" 11 68881416,69455872

POPULATION=$1 # comma-separated list of case,control haplotype_estimates files
HAPLOTYPES=$2 # three-column file with iteration number in first column, pattern length in second, pattern name in third
CHR=$3 # chromosome number
BP_RANGE=$4 # comma-separated list from-bp,to-bp (separate with semicolons or commas for multiple ranges; use quotes)
INDEX=$5 # optional index of the form STEP_SIZE.STEP_NO to indicate the range (STEP_NO-1)*STEP_SIZE to STEP_NO*STEP_SIZE-1 of haplotypes to analyze
ITERATION=$6 # optional overlap iteration at which pattern was found
DIRECTORY=$7 # location to store output/where .indiv files are found
HOME_DIR=$8 # location of program files

module load R/3.6.1
module load python/3.7.0
module load python/conda/3.8.1
module load conda3/5.1.0

printf "\n"
if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi
echo Home directory is $HOME_DIR

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY
echo Directory path is $DIRECTORY
printf "\n"

BP_RANGE=($(echo $BP_RANGE | perl -pne 's/([0-9]+)[,]*/$1 /g'))

POPULATION=($(echo $POPULATION | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))

if [ -f $HAPLOTYPES ]
then
  # if HAPLOTYPES supplied as the third column of a file
  if [ ! -z $ITERATION ] # if ITERATION is not the empty string ""
  then
    HAPLOTYPES_LIST=($(awk '($1 == '$ITERATION' ){printf "%s ",$3}' $HAPLOTYPES)) # select only haplotypes found at iteration ITERATION
  else
    HAPLOTYPES_LIST=($(awk '{printf "%s ",$3}' $HAPLOTYPES)) # select haplotypes at all iterations
  fi
else
  # if HAPLOTYPES supplied as a colon- (:) separated list
  HAPLOTYPES_LIST=($(echo $HAPLOTYPES | perl -pne 's/([0-9A-Za-z_.,=-]+)[:;]*/$1 /g')) # separate snp sets by (:) or (;), but not (,)
fi

if [ ! -z $INDEX ]
then
  INDEX="${INDEX}."
  INDEX_ARRAY=($(echo $INDEX | perl -pne 's/([0-9]+)[.]*/$1 /g')) # array (STEP_SIZE STEP_NO)
else
  INDEX=""
  INDEX_ARRAY=(${#HAPLOTYPES_LIST[@]} 1)
fi
ll=$(((${INDEX_ARRAY[1]}-1)*${INDEX_ARRAY[0]}))
ul=$((${INDEX_ARRAY[1]}*${INDEX_ARRAY[0]}-1))
(( $ul>$((${#HAPLOTYPES_LIST[@]}-1)) )) && ul=$((${#HAPLOTYPES_LIST[@]}-1))
echo Range: \(${ll},${ul}\)
printf "\n"

test -f haplotype_segregation.patterns_${ll}-${ul}.txt && echo rm haplotype_segregation.patterns_${ll}-${ul}.txt && printf "\n"
test -f haplotype_segregation.patterns_${ll}-${ul}.txt && rm haplotype_segregation.patterns_${ll}-${ul}.txt
touch haplotype_segregation.patterns_${ll}-${ul}.txt

for ((i=$ll;i<=$ul;i++))
do
  test -f haplotype_counts.pattern_${i}.txt && echo rm haplotype_counts.pattern_${i}.txt && printf "\n"
  test -f haplotype_counts.pattern_${i}.txt && rm haplotype_counts.pattern_${i}.txt
  touch haplotype_counts.pattern_${i}.txt

  for ((j=0;j<${#POPULATION[@]};j++))
  do
    # if haplotypes are in the form rsid1_[ATCG]=[01],rsid2_[ATCG]=[01]
    # echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${HAPLOTYPES_LIST[i]}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; snps[a]=b} } NR==1{ for \(i = 2\;i\<=NF\;i++\) {if \(\$i in snps\) {col[i]=\$i} } print \"sjlife\",\"${HAPLOTYPES_LIST[i]}\",\"affected\" } NR\>1{count=1\; for \(j in col\) {if \(\$j!=snps[col[j]]\) {count=0\; break} }\; print \$1,count,\'$((${#POPULATION[@]}-j-1))\' }\' ${POPULATION[j]} \>\> haplotype_counts.pattern_${i}.txt
    # awk 'BEGIN{OFS="\t"; n=split("'${HAPLOTYPES_LIST[i]}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); snps[a]=b} } NR==1{ for (i = 2;i<=NF;i++) {if ($i in snps) {col[i]=$i} } print "sjlife","'${HAPLOTYPES_LIST[i]}'","affected" } NR>1{count=1; for (j in col) {if ($j!=snps[col[j]]) {count=0; break} }; print $1,count,'$((${#POPULATION[@]}-j-1))' }' ${POPULATION[j]} >> haplotype_counts.pattern_${i}.txt

    # if haplotypes are in the form 01_[01],02_[01],...
    echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${HAPLOTYPES_LIST[i]}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"0*\([1-9][0-9]*\)_\([0-9]*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"0*\([1-9][0-9]*\)_\([0-9]*\)\",\"\\\\2\",\"g\",array[i]\)\; snps[a+1]=b} } NR==1{ print \"sjlife\",\"${HAPLOTYPES_LIST[i]}\",\"affected\" } NR\>1{count=1\; for \(j in snps\) {if \(\$j!=snps[j]\) {count=0\; break} }\; print \$1,count,$((${#POPULATION[@]}-j-1)) }\' ${POPULATION[j]} \>\> haplotype_counts.pattern_${i}.txt
    awk 'BEGIN{OFS="\t"; n=split("'${HAPLOTYPES_LIST[i]}'",array,","); for (i=1;i<=n;i++) {a=gensub("0*([1-9][0-9]*)_([0-9]*)","\\1","g",array[i]); b=gensub("0*([1-9][0-9]*)_([0-9]*)","\\2","g",array[i]); snps[a+1]=b} } NR==1{ print "sjlife","'${HAPLOTYPES_LIST[i]}'","affected" } NR>1{count=1; for (j in snps) {if ($j!=snps[j]) {count=0; break} }; print $1,count,'$((${#POPULATION[@]}-j-1))' }' ${POPULATION[j]} >> haplotype_counts.pattern_${i}.txt
    printf "\n"
  done
  echo Get segregation of pattern $i between populations
  echo awk \'BEGIN{OFS=\"\\t\"} NR==1{hap=\$2} \(NR\>1 \&\& \$3==1\){a+=\$2\; n+=\$3} \(NR\>1 \&\& \$3==0\){c+=\$2\; m+=\(1-\$3\)} END{b=n-a\; d=m-c\; print hap,a,b,c,d}\' haplotype_counts.pattern_${i}.txt \>\> haplotype_segregation.patterns_${ll}-${ul}.txt
  awk 'BEGIN{OFS="\t"} NR==1{hap=$2} (NR>1 && $3==1){a+=$2; n+=$3} (NR>1 && $3==0){c+=$2; m+=(1-$3)} END{b=n-a; d=m-c; print hap,a,b,c,d}' haplotype_counts.pattern_${i}.txt >> haplotype_segregation.patterns_${ll}-${ul}.txt # a = cases counts, b = controls counts, c = cases non-counts, d = controls non-counts of haplotype hap
  printf "\n"

  test -f haplotype_counts.pattern_${i}.txt && echo rm haplotype_counts.pattern_${i}.txt && printf "\n"
  test -f haplotype_counts.pattern_${i}.txt && rm haplotype_counts.pattern_${i}.txt
done

echo python ${HOME_DIR}/ukbb_fisher_exact.py -f haplotype_segregation.patterns_${ll}-${ul}.txt
python ${HOME_DIR}/ukbb_fisher_exact.py -f haplotype_segregation.patterns_${ll}-${ul}.txt

if [ -f $HAPLOTYPES ] && [ -f fisher_exact.patterns_${ll}-${ul}.txt ]
then
  echo Get pattern iteration:
  echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{iteration[\$3]=\$1\; next} \(\$1 in iteration\){\$1=sprintf\(\"%s\\t%s\"\,\$1\,iteration[\$1]\)\; print \$0}\' $HAPLOTYPES fisher_exact.patterns_${ll}-${ul}.txt \> fisher_exact.patterns_${ll}-${ul}.tmp
  awk 'BEGIN{OFS="\t"} NR==FNR{iteration[$3]=$1; next} ($1 in iteration){$1=sprintf("%s\t%s",$1,iteration[$1]); print $0}' $HAPLOTYPES fisher_exact.patterns_${ll}-${ul}.txt > fisher_exact.patterns_${ll}-${ul}.tmp # inserts iteration number at second field
  printf "\n"

else
  echo Insert NA as iteration:
  echo awk \'BEGIN{OFS=\"\\t\"} {\$1=sprintf\(\"%s\\t%s\"\,\$1\,\"NA\"\)\; print \$0}\' fisher_exact.patterns_${ll}-${ul}.txt \> fisher_exact.patterns_${ll}-${ul}.tmp
  awk 'BEGIN{OFS="\t"} {$1=sprintf("%s\t%s",$1,"NA"); print $0}' fisher_exact.patterns_${ll}-${ul}.txt > fisher_exact.patterns_${ll}-${ul}.tmp # inserts NA at second field
fi

test -f fisher_exact.patterns_${ll}-${ul}.tmp && echo mv fisher_exact.patterns_${ll}-${ul}.tmp fisher_exact.patterns_${ll}-${ul}.txt && printf "\n"
test -f fisher_exact.patterns_${ll}-${ul}.tmp && mv fisher_exact.patterns_${ll}-${ul}.tmp fisher_exact.patterns_${ll}-${ul}.txt
