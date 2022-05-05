#! /bin/bash

POPULATION=$1 # comma-separated list of haplotype_estimates_transpose files
CHR=$2 # chromosome number
BP_RANGE=$3 # comma-separated list "from,to" where from < to on chromosome
SIGMA=$4 # (one-fewer than the) initial number of chromosomes to overlap
TUPLE_RANGE=$5 # comma-separated list "from,to" of tuples of chromosomes to be overlapped, out of a total (2n choose sigma_0); OR a period-separated list "step_size.step_number" defining the range "(i-1) * step + 1,i * step"
NAMES=$6 # optional comma-separated list of population names in output file
DIRECTORY=$7 # Folder to store output Pattern file, e.g. ./snp_files; use $PWD if none supplied
HOME_DIR=$8 # location of program files

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
# Do NOT cd $DIRECTORY, because POPULATION files are in the submission directory, not the output directory

POPULATION=($(echo $POPULATION | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))

if [ ! -z $NAMES ]
then
  NAMES=($(echo $NAMES | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))
else
  for ((i=0;i<${#POPULATION[@]};i++))
  do
    NAMES+=("population_${i}")
  done
fi

module load vcftools/0.1.13
module load R/3.6.1

BP_RANGE=($(echo $BP_RANGE | perl -pne 's/([0-9]+)[,]*/$1 /g'))
if [ ! -z $TUPLE_RANGE ]
then
  TUPLE_RANGE_IN=$TUPLE_RANGE
fi

if [ -z $SIGMA ]
then
  SIGMA=2 #default is to join 2 chromosomes to every chromosome, creating triples
fi

for ((i=0;i<${#POPULATION[@]};i++))
do
  echo Population ${POPULATION[i]}
  printf "\n"
  # evaluate binomial coefficient (2n choose sigma) for n cases; total number of haplotypes might not be an even number due to rounding error in Rscript
  #https://stackoverflow.com/questions/25586346/bash-script-for-finding-binomical-coefficient
  let num0=$(awk 'BEGIN{nf=0} {if (NF>nf) {nf=NF} } END{print nf}' ${POPULATION[i]})
  let num0=num0-1 # total number of haplotypes; first column is rsid
  if (($num0>=$SIGMA))
  then
    let num=1;
    let den=1
    for ((j=1; j<=$SIGMA; j++));
    do
      let num=$num*$((num0-j+1)) # error if negative
      let den=$den*$j
    done
    let n_tuples=($num/$den)
  else
    n_tuples=0
  fi

  if [ -z $TUPLE_RANGE_IN ]
  then
    ll=0 # lower limit, first sigma-tuple of chromosomes, indexed from 0
    let ul=$n_tuples-1 # index starting from 0
    TUPLE_RANGE=$(echo ${ll},${ul})

    echo Evaluate tuples in range {$TUPLE_RANGE}.
    TUPLE_RANGE=($(echo $TUPLE_RANGE | perl -pne 's/([0-9]+)[.]+/$1 /g'))
  else
    # for the case that TUPLE_RANGE is in the format STEP_SIZE.STEP_NUMBER
    TUPLE_RANGE=($(echo $TUPLE_RANGE_IN | perl -pne 's/([0-9]+)[.]+/$1 /g'))
  fi

  echo ${TUPLE_RANGE[@]}
  if ((${#TUPLE_RANGE[@]}>1)) # length=1 if TUPLE_RANGE is a comma-separated list
  then
    let ll=$((TUPLE_RANGE[1]-1))*$((TUPLE_RANGE[0])); # (i - 1) * step_size
    let ll=$(($ll<=$((n_tuples-1))?$ll:$((n_tuples-1)))) # in case index i exceeds n_tuples
    let ul=$((TUPLE_RANGE[1]*TUPLE_RANGE[0]))-1 # i * step_size
    let ul=$(($ul<=$((n_tuples-1))?$ul:$((n_tuples-1)))) # in case index i exceeds n_tuples
    # https://stackoverflow.com/questions/10415064/how-to-calculate-the-minimum-of-two-variables-simply-in-bash
    TUPLE_RANGE=$(echo ${ll},${ul})
    echo $TUPLE_RANGE
    printf "\n"
  fi
  # two-element array of the lower and upper ranges of sigma0-tuples to sample
  TUPLE_RANGE=($(echo $TUPLE_RANGE | perl -pne 's/([0-9]+)[,]+/$1 /g'))
  echo ${TUPLE_RANGE[@]}
  let ll=$((${TUPLE_RANGE[0]}))
  let ll=$(($ll<=$((n_tuples-1))?$ll:$((n_tuples-1))))
  let ul=$((${TUPLE_RANGE[1]}))
  let ul=$(($ul<=$((n_tuples-1))?$ul:$((n_tuples-1))))
  TUPLE_RANGE=$(echo ${ll},${ul})
  TUPLE_RANGE=($(echo $TUPLE_RANGE | perl -pne 's/([0-9]+)[,]+/$1 /g'))
  echo Tuples:
  declare -p TUPLE_RANGE
  printf "\n"
  # Do first found of overlaps
  if ((${TUPLE_RANGE[1]}>${TUPLE_RANGE[0]})); #will not be true if both $ll and $ul exceed n_tuples
  then
    # fixed character length for tuples field in output file name
    printf -v tl "%0.${#n_tuples}d" ${TUPLE_RANGE[0]}
    printf -v tu "%0.${#n_tuples}d" ${TUPLE_RANGE[1]}
    echo Now running sh ${HOME_DIR}/subject_overlap_loop2.sh ${POPULATION[i]} $SIGMA ${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.tuples.${tl}-${tu} ${TUPLE_RANGE[0]},${TUPLE_RANGE[1]} $DIRECTORY $HOME_DIR
    sh ${HOME_DIR}/subject_overlap_loop2.sh ${POPULATION[i]} $SIGMA ${NAMES[i]}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}.tuples.${tl}-${tu} ${TUPLE_RANGE[0]},${TUPLE_RANGE[1]} $DIRECTORY $HOME_DIR

  else
    # echo Tuple index \(${TUPLE_RANGE[0]},${TUPLE_RANGE[1]}\) out of range
    echo Tuple index ${TUPLE_RANGE[0]} out of range.
    printf "\n"
  fi
  test -f ${POPULATION[i]} && echo rm ${POPULATION[i]} && printf "\n"
  test -f ${POPULATION[i]} && rm ${POPULATION[i]}
done
