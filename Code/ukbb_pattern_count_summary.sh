#! /bin/bash
set -e

# sh /home/wletsou/scripts/ukbb_pattern_count_summary.sh ukbb_bca_cases 11 69231642,69431642 2,j

# sh /home/wletsou/scripts/ukbb_pattern_count_summary.sh ukbb_bca_cases 11 68850000,69231641 2,j

# sh /home/wletsou/scripts/ukbb_pattern_count_summary.sh bca_pre_cases 11 69231642,69431642 2,j

# sh /home/wletsou/scripts/ukbb_pattern_count_summary.sh ukbb_bca_cases_complete.row_sample_4409 11 69231642,69431642 2,j

POPULATION=$1 # sinlge population name as it appears in Pattern_combined files
CHR=$2 # single chromosome number
BP_RANGE=$3 # comma-separated list from_bp,to_bp of regions on chromosome CHR
PATTERN=$4 # bar group type, e.g. 2,j or 2,3,j or 2,3+j
DIRECTORY=$5 # current working directory, specify replicate folder if REPLICATE_NO > 0
HOME_DIR=$6 # location of scripts to be run

module load R/3.6.1

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi

BP_RANGE=($(echo $BP_RANGE | perl -pne 's/([0-9]+)[,]*/$1 /g'))

test -f Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.txt && rm Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.txt
touch Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.txt
# for file in Pattern_combined_old.Iteration*.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.txt
for file in Pattern_combined.Iteration*.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.txt
do
  if [ -f $file ]
  then
    Iteration=${file#*.Iteration}
    Iteration=${Iteration%%.*}
    printf "%s\t%s\n" $Iteration $(cat $file | wc -l) >> Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.txt
  fi
done

for file in Pattern_combined.Iteration*.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.txt
do
  if [ -f $file ]
  then
    Iteration=${file#*.Iteration}
    Iteration=${Iteration%%.*}
    awk 'BEGIN{OFS="\t"} ($1 ~ /'$Iteration'/){print $0,'$(cat $file | wc -l)'} !($1 ~ /'$Iteration'/){print $0}' Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.txt > Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.temp.txt
    test -f Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.temp.txt && mv Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.temp.txt Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.txt
  fi
done

for file in fisher_exact.${POPULATION}.Iteration*.txt
do
  if [ -f $file ]
  then
    Iteration=${file#*.Iteration}
    Iteration=${Iteration%%.*}
    awk 'BEGIN{OFS="\t"} ($1 ~ /'$Iteration'/){print $0,'$(cat $file | wc -l)'} !($1 ~ /'$Iteration'/){print $0}' Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.txt > Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.temp.txt
    test -f Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.temp.txt && mv Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.temp.txt Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.txt
  fi
done


Rscript ${HOME_DIR}/ukbb_pattern_count_summary_plot.R file=Pattern_counts.${POPULATION}.chr${CHR}.${BP_RANGE[0]}-${BP_RANGE[1]}_${PATTERN}.txt directory=$DIRECTORY
