#! /bin/bash
set -e

# bsub -P SJLIFE -J ukbb_haplotype_model9_sub -oo ukbb_haplotype_model9_sub.out -eo ukbb_haplotype_model9_sub.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.sh bca.BR.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69231641.txt fisher_exact.ukbb_bca_cases.Results.translated.txt 1 \"rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs76987532_T=0,rs1122316_A=0,rs2298764_C=0,rs117222887_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0\" 1 11 68850000,69231641 50 \"\" \"\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.33/phase2 /home/wletsou/scripts"

PHENOTYPES=$1 # e.g. bca.BR.txt
POPULATION=$2 # combined cases/controls haplotype_estimates file
TEST_HAPLOTYPES=$3 # file with translated haplotypes in first column, or a colon-separated list of comma-separated lists rsid1_allele=[0,1] to be joined to INCLUDED_HAPLOTYPES
JOIN=$4 # comma-separated list of which INCLUDED_HAPLOTYPES to join to (start with "1", correction will be done for 0-based indexing); use all if value is empty
INCLUDED_HAPLOTYPES=$5 # colon-separated list or file of comma-separated lists rsid1_allele=[0,1]
GROUPING=$6 # comma-separated list 1,2-4,5,6-8,... of how the INCLUDED_HAPLOTYPES should be grouped in the model; single-grouping if empty
CHR=$7 # chromosome number
BP_RANGE=$8 # from_bp,to_bp
DELTA=$9 # number of combinations to do in one job
OUTPUT=${10} # optional name of output file
NAME=${11}
DIRECTORY=${12}
HOME_DIR=${13}

module load R/3.6.1

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
echo cd $DIRECTORY
cd $DIRECTORY

MAX_JOBS=$(cat "/hpcf/lsf/lsf_prod/conf/lsbatch/hpcf_research_cluster/configdir/lsb.users" | grep "#wletsou" | awk '{print $2}')
if [ -z $MAX_JOBS ]
then
  MAX_JOBS=100 # revised limit
fi

test -f $PHENOTYPES || (>&2 echo "Phenotypes file does not exist in ${DIRECTORY}"; exit 1)

BP_RANGE=($(echo $BP_RANGE | perl -pne 's/([0-9]+)[,]*/$1 /g'))

POPULATION=($(echo $POPULATION | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))

SNP_LIST=$(awk 'NR==1{for (j=2;j<=NF;j++) {printf "%s%s",$j,(j<NF?",":"")} }' ${POPULATION[0]}) # get SNPs in order from first haplotype_estimates file
char=$(($(echo $(($(echo $SNP_LIST | tr -cd "," | wc -c)+1)) | wc -c)-1)) # number of characters in the total number snps

if [ -f $TEST_HAPLOTYPES ] && [ ! -z $TEST_HAPLOTYPES ] # list can be supplied from a file, one haplotype per line
then
  n_haplotypes=$(cat $TEST_HAPLOTYPES | wc -l)
else
  array=($(echo $TEST_HAPLOTYPES | perl -pne 's/[:]/ /g'))
  n_haplotypes=${#array[@]}
fi

char_h=$(( $(echo $n_haplotypes | wc -c)-1)) # number of characters in the total number haplotypes
printf "\n"

if [ -f $INCLUDED_HAPLOTYPES ] && [ ! -z $INCLUDED_HAPLOTYPES ] # list can be supplied from a file
then
  INCLUDED_HAPLOTYPES=$(cat $INCLUDED_HAPLOTYPES) # supplied as a colon-separated list in a file
fi
INCLUDED_HAPLOTYPES=($(echo $INCLUDED_HAPLOTYPES | perl -pne 's/[:]/ /g')) # colon-separated list into an array
declare -p INCLUDED_HAPLOTYPES
printf "\n"

if [ ! -z $JOIN ]
then
  JOIN=($(echo $JOIN | perl -pne 's/[,]/ /g')) # array of haplotypes to join test haplotype to
else
  JOIN=($(eval "echo {1..$( ((${#INCLUDED_HAPLOTYPES[@]}>0)) && echo ${#INCLUDED_HAPLOTYPES[@]} || echo 1 )}")) # indicates joining of test haplotype to all included haplotypes
fi
# correction for 0-based indexing
for ((i=0;i<${#JOIN[@]};i++))
do
  JOIN[i]=$((JOIN[i]-1))
done
declare -p JOIN

group_array=($(echo $GROUPING | perl -pne 's/[,]/ /g'))
declare -p group_array
printf "\n"

for ((k=0;k<${#group_array[@]};k++))
do
  if ((${#JOIN[@]} > 1))
  then
    VAR="$((${JOIN[0]}+1))-$((${JOIN[$((${#JOIN[@]}-1))]}+1))" # express joining list as a range
  else
    VAR="$((${JOIN[0]}+1))" # express joining list as a single number
  fi
  if [[ "$VAR" == "${group_array[k]}" ]] # determine to which element of GROUPING the list of joining haplotypes belongs
  then
    group_no=$((k+1)) # correction for 0-based indexing
    break
  else
    group_no=0
  fi
done
test -z $group_no && group_no=0 || (echo Joining test haplotypes to group_h${group_no} of included haplotypes. && printf "\n")

if [ ! -z $NAME ]
then
  NAME=${NAME}"."
fi

if ((${#INCLUDED_HAPLOTYPES[@]}==0))
then
  # if there are no INCLUDED_HAPLOTYPES, print just list of chromosomes 1 and 2 for each subject
  echo awk \'NR==1{print \"sid\"} NR\>1{seen[\$1]+=1\; printf \"%s_%s\\n\",\$1,seen[\$1]}\' ${POPULATION} \> ${NAME}allele_counts.txt
  awk 'NR==1{print "sid"} NR>1{seen[$1]+=1; printf "%s_%s\n",$1,seen[$1]}' ${POPULATION} > ${NAME}allele_counts.txt && printf "\n"
else
  for ((i=0;i<${#INCLUDED_HAPLOTYPES[@]};i++))
  do
    # get haplotype counts (0 or 1) for each chromosome
    echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${INCLUDED_HAPLOTYPES[i]}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; haplotype_snps[a]=b} } NR==1{ for \(i = 2\;i\<=NF\;i++\) {if \(\$i in haplotype_snps\) {hap_col[i]=\$i} }\; if \(length\(hap_col\)!=n\) {exit 1}\; print \"sid\",\"${INCLUDED_HAPLOTYPES[i]}\"} NR\>1{hap_count=1\; for \(j in hap_col\) {if \(\$j!=haplotype_snps[hap_col[j]]\) {hap_count=0}\; if \(\$j !~ /[0-9]/\) {hap_count=\"NA\"} }\; print \$1,hap_count }\' ${POPULATION[0]} \> ${NAME}allele_counts.haplotype_${i}.txt
    awk 'BEGIN{OFS="\t"; n=split("'${INCLUDED_HAPLOTYPES[i]}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } NR==1{ for (i = 2;i<=NF;i++) {if ($i in haplotype_snps) {hap_col[i]=$i} }; if (length(hap_col)!=n) {exit 1}; print "sid","'${INCLUDED_HAPLOTYPES[i]}'"} NR>1{seen[$1]+=1;hap_count=1; for (j in hap_col) {if ($j!=haplotype_snps[hap_col[j]]) {hap_count=0}; if ($j !~ /[0-9]/) {hap_count="NA"} }; printf "%s_%s\t%s\n",$1,seen[$1],hap_count }' ${POPULATION[0]} > ${NAME}allele_counts.haplotype_${i}.txt
    printf "\n"

    if ((i==0))
    then
      echo mv ${NAME}allele_counts.haplotype_${i}.txt ${NAME}allele_counts.txt && printf "\n"
      mv ${NAME}allele_counts.haplotype_${i}.txt ${NAME}allele_counts.txt
    else
      # concatenate new column to haplotypes count table
      echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{seen[\$1]=\$2\; next} \(\$1 in seen\){print \$0,seen[\$1]}\' ${NAME}allele_counts.haplotype_${i}.txt ${NAME}allele_counts.txt \> ${NAME}allele_counts.tmp
      awk 'BEGIN{OFS="\t"} NR==FNR{seen[$1]=$2; next} ($1 in seen){print $0,seen[$1]}' ${NAME}allele_counts.haplotype_${i}.txt ${NAME}allele_counts.txt > ${NAME}allele_counts.tmp && printf "\n"

      test -f ${NAME}allele_counts.tmp && echo mv ${NAME}allele_counts.tmp ${NAME}allele_counts.txt
      test -f ${NAME}allele_counts.tmp && mv ${NAME}allele_counts.tmp ${NAME}allele_counts.txt && printf "\n"
    fi
  done
fi

if [ -z $DELTA ]
then
  DELTA=50 # initial number of haplotypes to count per job
fi

# calculate number of jobs for the TEST_HAPLOTYPES
n_jobs=$(awk 'BEGIN{printf "%0.25f\n",'$n_haplotypes'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
echo Number of jobs for $n_haplotypes haplotype$( (($n_haplotypes>1)) && echo "s" || echo "" ) is $n_jobs with DELTA = $DELTA haplotype$( (($DELTA>1)) && echo "s" || echo "" ) per job.
while (( $n_jobs > $MAX_JOBS ))
do
  DELTA=$((2*DELTA))
  n_jobs=$(awk 'BEGIN{printf "%0.25f\n",'$n_haplotypes'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
  echo Revised number of jobs is $n_jobs with DELTA = $DELTA haplotype$( (($DELTA>1)) && echo "s" || echo "" ) per job.
done
printf "\n"

for ((i=0;i<${#JOIN[@]};i++))
do
  echo bsub \-P SJLIFE \-J \"myJob[1-$n_jobs]\" \-oo ukbb_haplotype_model9.%I.out \-eo ukbb_haplotype_model9.%I.err \-R \"rusage[mem=1000]\" \-R \"select[ut \< 0.8]\" \-R \"order[!ut]\" \-q standard \"${HOME_DIR}/ukbb_haplotype_model9.sh ${POPULATION} $TEST_HAPLOTYPES \\\"${INCLUDED_HAPLOTYPES[${JOIN[i]}]}\\\" $DELTA.\\\$LSB_JOBINDEX \\\"$( ((${#JOIN[@]}>0)) && echo "included_haplotype_${JOIN[i]}" || echo "")\\\" $DIRECTORY $HOME_DIR\"
  bsub -P SJLIFE -J "myJob[1-$n_jobs]" -oo ukbb_haplotype_model9.%I.out -eo ukbb_haplotype_model9.%I.err -R "rusage[mem=1000]" -R "select[ut < 0.8]" -R "order[!ut]" -q standard "${HOME_DIR}/ukbb_haplotype_model9.sh ${POPULATION} $TEST_HAPLOTYPES \"${INCLUDED_HAPLOTYPES[${JOIN[i]}]}\" $DELTA.\$LSB_JOBINDEX \"$( ((${#JOIN[@]}>0)) && echo "included_haplotype_${JOIN[i]}" || echo "")\" $DIRECTORY $HOME_DIR"
  printf "\n"

  # wait until all jobs complete
  job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}')) # get job name from 7th field (or 6th if no exectution host yet) in all non-header rows of bjobs
  declare -p job_array
  while (( ${#job_array[@]}>0 ))
  do
    sleep 60
    job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}'))
    # declare -p job_array
  done
  printf "\n"
done

# concatenate results horizontally by subject's chromosome 1 or 2 (first column)
# for file in *new_allele_counts.*.txt
# do
#   if [ -f $file ]
#   then
#     echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{str=\"\"\; for \(i=2\;i\<=NF\;i++\) {str=sprintf\(\"%s%s\",\(i\>2?\"\\t\":\"\"\),str\)} seen[\$1]=str} \(\$1 in seen\){print \$0,seen[\$1]}\' ${file} ${NAME}allele_counts.txt \> ${NAME}allele_counts.tmp
#     awk 'BEGIN{OFS="\t"} NR==FNR{str=""; for (i=2;i<=NF;i++) {str=sprintf("%s%s",(i>2?"\t":""),str)} seen[$1]=str} ($1 in seen){print $0,seen[$1]}' ${file} ${NAME}allele_counts.txt > ${NAME}allele_counts.tmp # join all but first field of current haplotype counts file to master list
#
#     test -f ${NAME}allele_counts.tmp && echo mv ${NAME}allele_counts.tmp ${NAME}allele_counts.txt
#     test -f ${NAME}allele_counts.tmp && mv ${NAME}allele_counts.tmp ${NAME}allele_counts.txt && printf "\n"
#
#     rm $file && echo rm $file && printf "\n"
#   fi
# done
