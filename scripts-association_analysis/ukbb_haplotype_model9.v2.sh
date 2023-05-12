#! /bin/bash
set -e

# For counting number of copies of INCLUDED_HAPLOTYPE (h1) + TEST_HAPLOTYPE (= h2) for each TEST_HAPLOTYPE (in supplied RANGE) for each subject

POPULATION=$1 # combined cases/controls haplotype_estimates file
TEST_HAPLOTYPES=$2 # file with translated haplotypes in first column, or a colon-separated list of comma-separated lists rsid1_allele=[0,1] to be joined to INCLUDED_HAPLOTYPES
INCLUDED_HAPLOTYPE=$3 # one haplotype to be joined to each TEST_HAPLOTYPE, comma-separated list of rsid1_allele=[0,1]
RANGE=$4 # Range (ll,ul) of haplotypes to do, from 0 to n_haplotypes - 1.  OR a period-separated list "STEP_SIZE.STEP_NO" defining the range "(i-1) * step,i * step - 1"
NAME=$5 # optional name prepended to output file
DIRECTORY=$6
HOME_DIR=$7

if [ -z $HOME_DIR ];
then
  unset HOME_DIR
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY

POPULATION=($(echo $POPULATION | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))

SNP_LIST=$(awk 'NR==1{for (j=2;j<=NF;j++) {printf "%s%s",$j,(j<NF?",":"")} }' ${POPULATION[0]}) # get SNPs in order from first haplotype_estimates file
char_s=$(($(echo $(($(echo $SNP_LIST | tr -cd "," | wc -c)+1)) | wc -c)-1)) # number of characters in the total number snps

if [ -z $NAME ]
then
  unset NAME
fi

if [ -f $TEST_HAPLOTYPES ] && [ ! -z $TEST_HAPLOTYPES ] # list can be supplied from a file
then
  TEST_HAPLOTYPES=($(cut -f1 $TEST_HAPLOTYPES))
else
  TEST_HAPLOTYPES=($(echo $TEST_HAPLOTYPES | perl -pne 's/[:]/ /g'))
fi
n_haplotypes=${#TEST_HAPLOTYPES[@]}
char_h=$(( $(echo $n_haplotypes | wc -c)-1)) # number of characters in the total number haplotypes
printf "\n"

RANGE=($(echo $RANGE | perl -pne 's/([0-9]+)[.]+/$1 /g'))
if ((${#RANGE[@]}>1)) # length=1 if RANGE is a comma-separated list, =2 if in the form step_size.step_no
then
  ll=$(($((RANGE[1]-1))*$((RANGE[0])))) # (i - 1) * step_size
  ul=$(($((RANGE[1]*RANGE[0]))-1)) # i * step_size
  RANGE=$(echo ${ll},${ul})
fi

# two-element array of the lower and upper ranges of sigma0-tuples to sample
RANGE=($(echo $RANGE | perl -pne 's/([0-9]+)[,]+/$1 /g'))
echo Supplied range is \(${RANGE[0]},${RANGE[1]}\).
printf "\n"
step_size=$((${RANGE[1]}-${RANGE[0]}+1))
STEP=$(awk 'BEGIN{print int('"${RANGE[0]}"'/'"$step_size"' + 1)}')

((${RANGE[0]}>$((n_haplotypes-1)))) && (>&2 echo "Range exceeds total number of haplotypes."; exit 1)
((${RANGE[1]}>$((n_haplotypes-1)))) && RANGE[1]=$((n_haplotypes-1))
STEP_SIZE=$((${RANGE[1]}-${RANGE[0]}+1)) # corrected step size (if going over the maximum)
# STEP=$(awk 'BEGIN{print int('"${RANGE[0]}"'/'"$STEP_SIZE"' + 1)}')
echo Corrected range is \(${RANGE[0]},${RANGE[1]}\).
printf "\n"

# new file for counts of haplotypes in RANGE[0] to RANGE[1]
test -f ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt && echo rm ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt
test -f ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt && rm ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt && printf "\n"
echo touch ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt
touch ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt && printf "\n"

for ((i=${RANGE[0]};i<=${RANGE[1]};i++))
do
  start=$(date +%s.%N) # start time in seconds.nanoseconds

  pattern=${INCLUDED_HAPLOTYPE}$( ((${#INCLUDED_HAPLOTYPE[@]}>0)) && echo "," || echo "")${TEST_HAPLOTYPES[i]} # join test haplotype to select included haplotypes
  pattern=$(echo $pattern | awk 'BEGIN{m=split("'$SNP_LIST'",snp_list,",")} {n=split($0,alleles,","); for (i=1;i<=n;i++) {for (j=1;j<=m;j++) {if (alleles[i] ~ snp_list[j]) {out[sprintf("%0'${char_s}'d",j)]=alleles[i]} } } } END{n=asorti(out,dest); for (i=1;i<=n;i++) {printf "%s%s",out[dest[i]],(i<n?",":"")} }') # sort pattern in haplotype order
  printf "\n"

  echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${pattern}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; haplotype_snps[a]=b} } NR==1{ for \(i=2\;i\<=NF\;i++\) {if \(\$i in haplotype_snps\) {hap_col[i]=\$i} }\; if \(length\(hap_col\)!=n\) {exit 1}\; print \"sid\",\"${pattern}\"} NR\>1{seen[\$1]+=1\; hap_count=1\; for \(j in hap_col\) {if \(\$j!=haplotype_snps[hap_col[j]]\) {hap_count=0}\; if \(\$j !~ /[0-9]/\) {hap_count=\"NA\"} }\; printf \"%s_%s\\t%s\\n\",\$1,seen[\$1],hap_count }\' ${POPULATION[0]} \> ${NAME/%/.}new_allele_counts.haplotype_${i}.txt
  awk 'BEGIN{OFS="\t"; n=split("'${pattern}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } NR==1{ for (i = 2;i<=NF;i++) {if ($i in haplotype_snps) {hap_col[i]=$i} }; if (length(hap_col)!=n) {exit 1}; print "sid","'${pattern}'"} NR>1{seen[$1]+=1; hap_count=1; for (j in hap_col) {if ($j!=haplotype_snps[hap_col[j]]) {hap_count=0}; if ($j !~ /[0-9]/) {hap_count="NA"} }; printf "%s_%s\t%s\n",$1,seen[$1],hap_count }' ${POPULATION[0]} > ${NAME/%/.}new_allele_counts.haplotype_${i}.txt
  printf "\n"

  if ((i==${RANGE[0]}))
  then
    echo mv ${NAME/%/.}new_allele_counts.haplotype_${i}.txt ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt && printf "\n"
    mv ${NAME/%/.}new_allele_counts.haplotype_${i}.txt ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt
  else
    # concatenate new column to haplotypes count table
    echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{seen[\$1]=\$2\; next} \(\$1 in seen\){print \$0,seen[\$1]}\' ${NAME/%/.}new_allele_counts.haplotype_${i}.txt ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt \> ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").tmp
    awk 'BEGIN{OFS="\t"} NR==FNR{seen[$1]=$2; next} ($1 in seen){print $0,seen[$1]}' ${NAME/%/.}new_allele_counts.haplotype_${i}.txt ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt > ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").tmp && printf "\n"

    test -f ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").tmp && echo mv ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").tmp ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt
    test -f ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").tmp && mv ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").tmp ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt && printf "\n"
  fi
  test -f ${NAME/%/.}new_allele_counts.haplotype_${i}.txt && echo rm ${NAME/%/.}new_allele_counts.haplotype_${i}.txt
  test -f ${NAME/%/.}new_allele_counts.haplotype_${i}.txt && rm ${NAME/%/.}new_allele_counts.haplotype_${i}.txt && printf "\n"

  finish=$(date +%s.%N) # end time in seconds.nanoseconds
  awk 'BEGIN{duration='$finish'-'$start';printf "Count for haplotype '$((i-${RANGE[0]}+1))' of '$STEP_SIZE' completed in %0.2f s.\n",duration}'
done

# if no output is produced, then print chromosomes 1 and 2 for all subjects
if [ ! -f ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt ]
then
  echo No output produced for \(${RANGE[0]},${RANGE[1]}\).
  echo awk \'NR==1{print \"sid\"} NR\>1{seen[\$1]+=1\; printf \"%s_%s\\n\",\$1,seen[\$1]}\' ${POPULATION} \> ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt
  awk 'NR==1{print "sid"} NR>1{seen[$1]+=1; printf "%s_%s\n",$1,seen[$1]}' ${POPULATION} > ${NAME/%/.}new_allele_counts.$(eval "printf 'job%d' $STEP").txt && printf "\n"
fi
