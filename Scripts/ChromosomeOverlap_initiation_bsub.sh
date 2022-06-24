#! /bin/bash
set -e

# loops through all SIGMA-tuples of chromosomes in the columns of HAPLOTYPES_FILE and look for SNPs that are the same across the group.  Output the unique patterns and their counts for each group type.  Relies on subject_overlap.sh to do overlaps and comparisons

# HOME_DIR/ChromosomeOverlap_initiation_bsub.sh haplotype_estimates_transpose.ukbb_bca_cases.chr11.69231642-69431642.txt 1 chr11.69231642-69431642 50 DIRECTORY HOME_DIR

HAPLOTYPES_FILE=$1 # transposed haplotypes file of snps x chromosomes
SIGMA=$2 # (One fewer than) the number of chromosomes to overlap
OUTPUT=$3 # Base file name for output Pattern.${OUTPUT}_*.txt
DELTA=$4 # number of overlaps to do per job
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
echo cd $DIRECTORY
cd $DIRECTORY && printf "\n"

[ ! -f $HAPLOTYPES_FILE ] && (>&2 echo "Haplotypes file not found."; exit 1)
[ -z $HAPLOTYPES_FILE ] && (>&2 echo "Haplotypes file not supplied."; exit 1)
n_samples=$(awk 'BEGIN{nf=0} {if (NF>nf) {nf=NF} } END{print nf-1}' $HAPLOTYPES_FILE) # columns (does not include rsid column)
n_snps=$(awk 'END{print NR}' $HAPLOTYPES_FILE)

echo Determine number of $SIGMA\-combinations of fixed chromosomes for $HAPLOTYPES_FILE:
num0=$n_samples
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
let n_tuples=($num/$den)
echo Number of operations to be performed is $n_tuples.
printf "\n"

if [ -z $MAX_JOBS ]
then
  MAX_JOBS=100 # revised limit
fi

N_JOBS=$(awk 'BEGIN{printf "%0.25f\n",'$n_tuples'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
echo Requested number of jobs for $n_tuples combination$( (($n_tuples>1)) && echo "s" || echo "" ) is $N_JOBS with DELTA = $DELTA combination$( (($DELTA>1)) && echo "s" || echo "" ) per job.
while (( $N_JOBS > $MAX_JOBS ))
do
DELTA=$((2*DELTA))
N_JOBS=$(awk 'BEGIN{printf "%0.25f\n",'$n_tuples'/'$DELTA'}' | awk '{if ($1 != int($1)) {$1=int(($1))+1} {printf "%.0f\n",$1}  }')
echo Revised number of jobs is $N_JOBS with DELTA = $DELTA combination$( (($DELTA>1)) && echo "s" || echo "" ) per job.
done
printf "\n"

echo Submit haplotype jobs:
echo bsub \-P SJLIFE \-J \"myJob[1-$N_JOBS]\" \-oo ${DIRECTORY}/${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_sub.%I.out \-eo ${DIRECTORY}/${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_sub.%I.err \-R \"rusage[mem=256]\" \-R \"select[ut \< 0.8]\" \-R \"order[!ut]\" \"sh ${HOME_DIR}ChromosomeOverlap_initiation_sub.sh ${HAPLOTYPES_FILE} $SIGMA \\\"${OUTPUT}\\\" ${DELTA}.\\\$LSB_JOBINDEX $DIRECTORY $HOME_DIR\"
job_id=$(bsub -P SJLIFE -J "myJob[1-$N_JOBS]" -oo ${DIRECTORY}/${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_sub.%I.out -eo ${DIRECTORY}/${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_sub.%I.err -R "rusage[mem=256]" -R "select[ut < 0.8]" -R "order[!ut]" "sh ${HOME_DIR}ChromosomeOverlap_initiation_sub.sh ${HAPLOTYPES_FILE} $SIGMA \"${OUTPUT}\" ${DELTA}.\$LSB_JOBINDEX $DIRECTORY $HOME_DIR")
printf "\n"
echo $job_id
job_id=$(echo $job_id | awk 'b=gensub(/.*<([0-9]*)>.*/,"\\1","g",$0) {print b}') #extract job_id (number) from output

echo Wait until all jobs start:
echo bsub \-P SJLIFE \-J sleep.${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_sub \-w \"numrun\($job_id,\*\) \|\| numended\($job_id,\*\)\" \-R \"rusage[mem=32]\" \-oo ${DIRECTORY}/sleep.${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_sub.out \-eo ${DIRECTORY}/sleep.${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_sub.err \-K \"sleep 10\"
bsub -P SJLIFE -J sleep.${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_sub -w "numrun($job_id,*) || numended($job_id,*)" -R "rusage[mem=32]" -oo ${DIRECTORY}/sleep.${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_sub.out -eo ${DIRECTORY}/sleep.${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_sub.err -K "sleep 10"
printf "\n"

echo Wait until all jobs done:
job_array=($(bjobs -J myJob 2> /dev/null | awk 'NR>1{print $7}')) # get job name from 7th field in all non-header rows of bjobs
declare -p job_array
while (( ${#job_array[@]}>0 ))
do
  sleep 60
  job_array=($(bjobs -J myJob 2> /dev/null | awk 'NR>1{print $7}'))
done
printf "\n"

echo Combine results from Iteration000:
echo bsub \-P SJLIFE \-J ${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_combine.Iteration000 \-oo ${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_combine.Iteration000.out \-eo ${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_combine.Iteration000.err \-R "rusage[mem=1024]" \-K "sh ${HOME_DIR}ChromosomeOverlap_initiation_combine.sh ${OUTPUT} $SIGMA \"Iteration000\" ${DIRECTORY}"
bsub -P SJLIFE -J ${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_combine.Iteration000 -oo ${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_combine.Iteration000.out -eo ${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_combine.Iteration000.err -R "rusage[mem=1024]" -K "sh ${HOME_DIR}ChromosomeOverlap_initiation_combine.sh ${OUTPUT} $SIGMA \"Iteration000\" ${DIRECTORY}"
printf "\n"

echo Wait until combinations complete:
echo bsub \-P SJLIFE \-J sleep.${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_combine.Iteration000 \-oo ${DIRECTORY}/sleep.${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_combine.Iteration000.out \-eo ${DIRECTORY}/sleep.${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_combine.Iteration000.err \-w \"done\(${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_combine.Iteration000\)\" \-R \"rusage[mem=32]\" \-K \"sleep 10\"
bsub -P SJLIFE -J sleep.${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_combine.Iteration000 -oo ${DIRECTORY}/sleep.${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_combine.Iteration000.out -eo ${DIRECTORY}/sleep.${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_combine.Iteration000.err -w "done(${OUTPUT}$( [ ! -z $OUTPUT ] && echo "." || echo "" )ChromosomeOverlap_initiation_combine.Iteration000)" -R "rusage[mem=32]" -K "sleep 10"
echo Done.
printf "\n"
