#! /bin/bash
set -e

# for reducing common haplotypes using rpart in R

# sh haplotype_rpart.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.69231642-69431642.txt "rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs2298764_C=0,rs117222887_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "" 1

POPULATION=$1 # haplotype estiamtes file
TEST_HAPLOTYPE=$2 # comma-separated list of alleles corresponding to columns of POPULATION to be used for scoring haplotype and for partitioning
INCLUDED_HAPLOTYPE=$3 # comma-separated list of alleles corresponding to columns of POPULATION to be included in haplotype definition but not for partitioning
MINSPLIT=$4 # rpart minsplit parameter, minimum number of items in split-off bin

if [ -z $HOME_DIR ];
then
  unset HOME_DIR
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY

module load R/3.6.1

if [ -z $MINSPLIT ]
then
  MINSPLIT=10
fi

POPULATION=($(echo $POPULATION | perl -pne 's/([^,]+)[,]*/$1 /g'))

TEST_HAPLOTYPE=($(echo $TEST_HAPLOTYPE | perl -pne 's/[:]/ /g'))

for ((i=0;i<${#POPULATION[@]};i++))
do
  for ((j=0;j<${#TEST_HAPLOTYPE[@]};j++))
  do
    echo Haplotype $j \(${TEST_HAPLOTYPE[j]}\):
    printf "\n"
    # get counts of TEST_HAPLOTYPE only if INCLUDED_HAPLOTYPE is included in line; print counts and all snp alleles
    echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${TEST_HAPLOTYPE[j]}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; test_snps[a]=b}\; delete array\; m=split\(\"${INCLUDED_HAPLOTYPE}\",array,\",\"\)\; for \(i=1\;i\<=m\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; included_snps[a]=b} } NR==1{ delete test_col\; delete included_col\; for \(i=2\;i\<=NF\;i++\) {if \(\$i in test_snps\) {test_col[i]=\$i}\; if \(\$i in included_snps\) {included_col[i]=\$i} }\; if \(length\(test_col\)!=n \|\| length\(included_col\)!=m\) {print \"Wrong number of alleles found.\"\; exit 1}\; printf \"id\\ttest_haplotype\\tincluded_haplotype\\t\" } NR\>1{ included_count=1\; test_count=0\; for \(j in included_col\) {if \(\$j!=included_snps[included_col[j]]\) {included_count=0}\; if \(\$j !~ /[0-9]/\) {included_count=\"NA\"} }\; if \(included_count==1\) {test_count=1\; for \(j in test_col\) {if \(\$j!=test_snps[test_col[j]]\) {test_count=0}\; if \(\$j !~ /[0-9]/\) {test_count=\"NA\"} } } printf \"%s\\t%s\\t%s\\t\",\$1,test_count,included_count } NR\>=1{for \(j=2\;j\<=NF\;j++\) {printf \"%s%s\",\$j,\(j\<NF?\"\\t\":\"\\n\"\)} }\' ${POPULATION[i]} \> haplotype_counts.${POPULATION[i]%.*}.txt
    printf "\n"
    awk 'BEGIN{OFS="\t"; n=split("'${TEST_HAPLOTYPE[j]}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); test_snps[a]=b}; delete array; m=split("'${INCLUDED_HAPLOTYPE}'",array,","); for (i=1;i<=m;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); included_snps[a]=b} } NR==1{ delete test_col; delete included_col; for (i=2;i<=NF;i++) {if ($i in test_snps) {test_col[i]=$i}; if ($i in included_snps) {included_col[i]=$i} }; if (length(test_col)!=n || length(included_col)!=m) {print "Wrong number of alleles found."; exit 1}; printf "id\ttest_haplotype\tincluded_haplotype\t" } NR>1{ included_count=1; test_count=0; for (j in included_col) {if ($j!=included_snps[included_col[j]]) {included_count=0}; if ($j !~ /[0-9]/) {included_count="NA"} }; if (included_count==1) {test_count=1; for (j in test_col) {if ($j!=test_snps[test_col[j]]) {test_count=0}; if ($j !~ /[0-9]/) {test_count="NA"} } }; printf "%s\t%s\t%s\t",$1,test_count,included_count } NR>=1{for (j=2;j<=NF;j++) {printf "%s%s",$j,(j<NF?"\t":"\n")} }' ${POPULATION[i]} > haplotype_counts.${POPULATION[i]%.*}.txt
    printf "\n"

    echo Rscript ${HOME_DIR/%/\/}haplotype_rpart.R file=haplotype_counts.${POPULATION[i]%.*}.txt haplotype=\"included_haplotype,${TEST_HAPLOTYPE[j]}\" min_split=$MINSPLIT
    Rscript ${HOME_DIR/%/\/}haplotype_rpart.R file=haplotype_counts.${POPULATION[i]%.*}.txt haplotype="included_haplotype,${TEST_HAPLOTYPE[j]}" min_split=$MINSPLIT
    printf "\n"

    test -f haplotype_counts.${POPULATION[i]%.*}.txt && echo rm haplotype_counts.${POPULATION[i]%.*}.txt && printf "\n"
    test -f haplotype_counts.${POPULATION[i]%.*}.txt && rm haplotype_counts.${POPULATION[i]%.*}.txt
  done
done
