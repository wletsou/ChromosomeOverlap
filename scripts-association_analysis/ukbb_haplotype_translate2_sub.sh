#! /bin/bash
set -e

# bsub -P SJLIFE -J ukbb_haplotype_translate2_sub -oo ukbb_haplotype_translate2_sub.out -eo ukbb_haplotype_translate2_sub.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr11.68850000-69231641.subset.txt fisher_exact.ukbb_bca_cases.Results.txt 50"

POPULATION=$1 # haplotype_estimates file for population with subject id in first column and rsids in columns 2 and beyond
HAPLOTYPES=$2 # colon-separated list of haplotypes, or a file with haplotypes in the first column.  Can be of the form column_allele,... or rsid_allele=[0/1]
DELTA=$3 # number of haplotypes to do in one job
DIRECTORY=$4 # current working directory
HOME_DIR=$5 # location of scripts to be run

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY

# MAX_JOBS=$(cat "/hpcf/lsf/lsf_prod/conf/lsbatch/hpcf_research_cluster/configdir/lsb.users" | grep "#wletsou" | awk '{print $2}')
if [ -z $MAX_JOBS ]
then
  MAX_JOBS=1500 # revised limit
fi

if [ -f $HAPLOTYPES ]
then
  n_haplotypes=$(cat $HAPLOTYPES | wc -l) # first column of file, assuming no header
else
  array=($(echo $HAPLOTYPES | perl -pne 's/[:]/ /g')) # colon-separated list into an array
  n_haplotype=${#array[@]}
fi

n_jobs=$(awk 'BEGIN{printf "%0.25f\n",'$n_haplotypes'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
echo Number of jobs for $n_haplotypes haplotype$( (($n_haplotypes>1)) && echo "s" || echo "" ) is $n_jobs with DELTA = $DELTA haplotype$( (($DELTA>1)) && echo "s" || echo "" ) per job.
while (( $n_jobs > $MAX_JOBS ))
do
  DELTA=$((2*DELTA))
  n_jobs=$(awk 'BEGIN{printf "%0.25f\n",'$n_haplotypes'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
  echo Revised number of jobs is $n_jobs with DELTA = $DELTA haplotype$( (($DELTA>1)) && echo "s" || echo "" ) per job.
done
printf "\n"

echo bsub \-P SJLIFE \-J \"myJob[1-$n_jobs]\" \-oo ukbb_haplotype_translate2.%I.out \-eo ukbb_haplotype_translate2.%I.err \-R \"rusage[mem=2048]\" \-q standard \"${HOME_DIR}/ukbb_haplotype_translate2.sh ${POPULATION} $HAPLOTYPES $DELTA.\\\$LSB_JOBINDEX $DIRECTORY $HOME_DIR\"
bsub -P SJLIFE -J "myJob[1-$n_jobs]" -oo ukbb_haplotype_translate2.%I.out -eo ukbb_haplotype_translate2.%I.err -R "rusage[mem=2048]" -q standard "${HOME_DIR}/ukbb_haplotype_translate2.sh ${POPULATION} $HAPLOTYPES $DELTA.\$LSB_JOBINDEX $DIRECTORY $HOME_DIR"
printf "\n"

job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}')) # get job name from 7th field (or 6th if no exectution host yet) in all non-header rows of bjobs
declare -p job_array
while (( ${#job_array[@]}>0 ))
do
  sleep 20
  job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}'))
  declare -p job_array
done
printf "\n"

test -f ${HAPLOTYPES%.*}.translated.txt && echo rm ${HAPLOTYPES%.*}.translated.txt
test -f ${HAPLOTYPES%.*}.translated.txt && rm ${HAPLOTYPES%.*}.translated.txt && printf "\n"

echo touch ${HAPLOTYPES%.*}.translated.txt
touch ${HAPLOTYPES%.*}.translated.txt && printf "\n"

for file in ${HAPLOTYPES%.*}.*.translated.txt
do
  if [ -f $file ]
  then
    cat $file >> ${DIRECTORY}/${HAPLOTYPES%.*}.translated.txt && echo cat $file \>\> ${DIRECTORY}/${HAPLOTYPES%.*}.translated.txt
    rm $file && echo rm $file
  fi
done
