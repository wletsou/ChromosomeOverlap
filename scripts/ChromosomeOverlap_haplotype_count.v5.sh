#! /bin/bash
set -e

# For each haplotype, prints a = cases counts, b = controls counts, c = cases non-counts, d = controls non-counts

# sh HOME_DIR/ChromosomeOverlap_haplotype_count.v2.sh haplotype_estimates.ukbb_bca_cases.chr11.69231642-69431642.txt,haplotype_estimates.ukbb_bca_controls.chr11.69231642-69431642.txt Pattern_combined.Iteration000.chr11.69231642-69431642_2,j.txt "" "chr11.69231642-69431642"

# Rscript HOME_DIR/ChromosomeOverlap_fisher_exact.R -f haplotype_segregation.chr11.69231642-69431642.patterns_0000-XXXX.txt -o chr11.69231642-69431642.patterns_0000-XXXX

POPULATION=$1 # comma-separated list of case,control haplotype_estimates files
HAPLOTYPES=$2 # multi-column file with (optional) pattern name in the last column
INDEX=$3 # optional index of the form STEP_SIZE.STEP_NO to indicate the range (STEP_NO-1)*STEP_SIZE to STEP_NO*STEP_SIZE-1 of haplotypes to analyze
NAME=$4 # optional name for output files
DIRECTORY=$5 # location to store output

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY
echo Directory path is $DIRECTORY
printf "\n"

POPULATION=($(echo $POPULATION | perl -pne 's/([^,]+)[,]*/$1 /g'))
(( ${#POPULATION[@]}== 2)) || (>&2 echo "Two POPULATION files must be supplied."; exit 1)

if [ ! -z $HAPLOTYPES ] && [ -f $HAPLOTYPES ]
then
  n_haplotypes=$(cat $HAPLOTYPES | wc -l)
  char=$(($( echo $n_haplotypes | wc -m)-1))
else
  test -z $HAPLOTYPES && (>&2 echo "HAPLOTYPES file must be supplied."; exit 1)
  test -f $HAPLOTYPES || (>&2 echo "HAPLOTYPES file not found."; exit 1)
fi

if [ -z $NAME ]
then
  unset NAME
fi

if [ ! -z $INDEX ]
then
  INDEX_ARRAY=($(echo $INDEX | perl -pne 's/([0-9]+)[.]*/$1 /g')) # array (STEP_SIZE STEP_NO)
else
  INDEX_ARRAY=(${n_haplotypes} 1) # do all haplotypes in one step
fi
ll=$(((${INDEX_ARRAY[1]}-1)*${INDEX_ARRAY[0]}+1))
ul=$((${INDEX_ARRAY[1]}*${INDEX_ARRAY[0]}))
echo Range: \(${ll},${ul}\)
printf "\n"

test -f ${DIRECTORY/%/\/}haplotype_segregation${NAME/#/.}.txt && echo rm ${DIRECTORY/%/\/}haplotype_segregation${NAME/#/.}.txt || printf ""
test -f ${DIRECTORY/%/\/}haplotype_segregation${NAME/#/.}.txt && rm ${DIRECTORY/%/\/}haplotype_segregation${NAME/#/.}.txt || printf ""
touch ${DIRECTORY/%/\/}haplotype_segregation${NAME/#/.}.txt && printf "\n" || printf ""

declare -a hap_array

i=1
j=$ll
awk '{print $NF}' $HAPLOTYPES | while read pattern
do
  start=$(date +%s.%N)
  if [[ ${pattern} =~ rs[0-9]+_[A-Z]+=[0-1]+ ]]
  then
    # if haplotypes are in the form rsid1_[ATCG]=[01],rsid2_[ATCG]=[01]
    echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${pattern}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; snps[a]=b} } NR==1{ for \(i = 2\;i\<=NF\;i++\) {if \(\$i in snps\) {col[i]=\$i} } } FNR\>1{count=1\; for \(j in col\) {if \(\$j !~ /[0-1]/\) {next}\; if \(\$j!=snps[col[j]]\) {count=0\; break} }\; hap[NR==FNR]+=count\; chrom[NR==FNR]+=1} END{a=hap[1]+0\; b=chrom[1]\-a+0\; c=hap[0]+0\; d=chrom[0]\-c+0\; print \"${pattern}\",a+0,b+0,c+0,d+0}\' ${POPULATION[*]}
    hap_array+=($(awk 'BEGIN{OFS="\t"; n=split("'${pattern}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); snps[a]=b} } NR==1{ for (i = 2;i<=NF;i++) {if ($i in snps) {col[i]=$i} } } FNR>1{count=1; for (j in col) {if ($j !~ /[0-1]/) {next}; if ($j!=snps[col[j]]) {count=0; break} }; hap[NR==FNR]+=count; chrom[NR==FNR]+=1} END{a=hap[1]+0; b=chrom[1]-a+0; c=hap[0]+0; d=chrom[0]-c+0; print "'${pattern}'",a+0,b+0,c+0,d+0}' ${POPULATION[*]})) # a = cases counts, b = controls counts, c = cases non-counts, d = controls non-counts of haplotype hap
  elif [[ ${pattern} =~ 0*[1-9][0-9]*_[0-9]+ ]]
  then
    # if haplotypes are in the form 01_[01],02_[01],...
    echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${pattern}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"0*\([1-9][0-9]*\)_\([0-9]*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"0*\([1-9][0-9]*\)_\([0-9]*\)\",\"\\\\2\",\"g\",array[i]\)\; snps[a+1]=b} } FNR\>1{count=1\; for \(j in snps\) {if \(\$j !~ /[0-1]/\) {next}\; if \(\$j!=snps[j]\) {count=0\; break} }\; hap[NR==FNR]+=count\; chrom[NR==FNR]+=1} END{a=hap[1]\; b=chrom[1]\-a\; c=hap[0]\; d=chrom[0]\-c\; print \"${pattern}\",a+0,b+0,c+0,d+0}\' ${POPULATION[*]}
    hap_array+=($(awk 'BEGIN{OFS="\t"; n=split("'${pattern}'",array,","); for (i=1;i<=n;i++) {a=gensub("0*([1-9][0-9]*)_([0-9]*)","\\1","g",array[i]); b=gensub("0*([1-9][0-9]*)_([0-9]*)","\\2","g",array[i]); snps[a+1]=b} } FNR>1{count=1; for (j in snps) {if ($j !~ /[0-1]/) {next}; if ($j!=snps[j]) {count=0; break} }; hap[NR==FNR]+=count; chrom[NR==FNR]+=1} END{a=hap[1]+0; b=chrom[1]-a+0; c=hap[0]+0; d=chrom[0]-c+0; print "'${pattern}'",a+0,b+0,c+0,d+0}' ${POPULATION[*]})) # a = cases counts, b = controls counts, c = cases non-counts, d = controls non-counts of haplotype hap
  else
    echo Pattern not recognized.
  fi
  finish=$(date +%s.%N)
  duration=$(awk 'BEGIN{printf "%0.2f\n", '$finish'-'$start'}')
  echo Haplotype $i took $duration second$( [ "$duration" == "1.00" ] && echo "" || echo "s").
  i=$((i+1))
  j=$((j+1))
  printf "\n"
  if (($i > $ul)) || (( $i > $n_haplotypes ))
  then
    echo Write haplotype counts:
    echo printf \"%s\\t%s\\t%s\\t%s\\t%s\\n\" \"\${hap_array[@]}\" \>\> ${DIRECTORY/%/\/}haplotype_segregation${NAME/#/.}.txt
    printf "%s\t%s\t%s\t%s\t%s\n" "${hap_array[@]}" >> ${DIRECTORY/%/\/}haplotype_segregation${NAME/#/.}.txt
  fi
done
