#! /bin/bash
set -e

# wrapper script for doing forward-selection in DRIVE on each haplotype in counts files produced by ukbb_haplotype_model9_sub.v2.sh

# bsub -P SJLIFE -J dbgap28544_haplotype_model9_iterate.v2.replication -oo dbgap28544_haplotype_model9_iterate.v2.replication.out -eo dbgap28544_haplotype_model9_iterate.v2.replication.err -R "rusage[mem=256]" "sh dbgap28544_haplotype_model9_iterate.v2.sh drive_replication.allele_counts.txt dbgap28544.pheno.txt \"\" \"\" 1e-5 drive_replication.Significant_patterns.txt \"drive_replication\""

INCLUDED_HAPLOTYPES=$1 # counts of included_haplotypes for each chromosome (can be just a list of chromosomes if no included_haplotypes in model)
PHENOTYPES=$2 # phenotypes file for getting age ("ageonset"), principal components of ancestry ("pc1" to "pc10") and affected status ("BCa") from each subject ("sid")
JOIN=$3 # comma-separated list of which INCLUDED_HAPLOTYPES to join to (start with "1", correction will be done for 0-based indexing); use all if value is empty
GROUPING=$4 # comma-separated list 1,2-4,5,6-8,... of how the INCLUDED_HAPLOTYPES should be grouped in the model; single-grouping if empty
CUTOFF=$5 # LRT p-value cutoff; stop when all p exceed
OUTPUT=$6 # optional name of output file
NAME=$7 # optional name to prepend to included_haplotypes files
DIRECTORY=$8
HOME_DIR=$9

if [ -z $HOME_DIR ];
then
  unset HOME_DIR
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY

if [ -z $OUTPUT ]
then
  OUTPUT="Significant_patterns.txt"
fi

awk 'BEGIN{printf "variable\thaplotype\tHR\tll\tul\tp\tfreq\tcount\tcases_freq\tcontrols_freq\tcases_count\tcontrols_count\n"}' > ${OUTPUT}

n_included_haplotypes=$(awk 'NR==1{print NF-1}' ${INCLUDED_HAPLOTYPES}) # number of new haplotypes to be tested

if [ ! -z $JOIN ]
then
  JOIN=($(echo $JOIN | perl -pne 's/[,]/ /g')) # array of haplotypes to join test haplotype to
else
  JOIN=() # indicates joining to the empty haplotype
fi
max_haplotype=$(awk 'BEGIN{max=0; n=split("'$(echo ${GROUPING} | sed 's/[-]/,/g')'",array,","); for (i=1;i<=n;i++) {a=gensub("([0-9]*)-[0-9]*","\\1","g",array[i]); if (a>max) {max=a}; b=gensub("[0-9]*-([0-9]*)","\\1","g",array[i]); if (b>max) {max=b} }; print max}' ) # maximum variable name in GROUPING array
if [ ! -z $n_included_haplotypes ]
then
  GROUPING=($(echo $GROUPING | sed 's/[,]/ /g')$(test ! -z $GROUPING && echo " " || echo "")$((($max_haplotype+1<=$n_included_haplotypes)) && eval "echo {$(($max_haplotype+1))..$( (($n_included_haplotypes>$max_haplotype)) && echo $n_included_haplotypes || echo $(($max_haplotype+1)) )}" || echo "")) # grouping into single-haplotype groups of included haplotypes which have not already been grouped
fi
declare -p GROUPING
groups=$(echo "${GROUPING[*]}" | sed 's/\s\+/,/g') # string of groups in baseline model
group_array=($(echo "${GROUPING[*]}")) # array of groups in baseline model
max_new_haplotypes=$(awk 'BEGIN{n=split("'$(echo ${GROUPING[@]} | sed 's/\s\+/,/g')'",array,","); for (i=1;i<=n;i++) {GROUPING[array[i]]}; delete array; m=split("'$(echo ${JOIN[@]} | sed 's/\s\+/,/g')'",array,","); for (i=1;i<=m;i++) {JOIN[array[i]]}; if (m>0) {for (i in JOIN) {for (j in GROUPING) {a=gensub("([0-9]*)-[0-9]*","\\1","g",j); b=gensub("[0-9]*-([0-9]*)","\\1","g",j); if (i>=a && i<=b) {print j} } } } else {print 0} }' | awk '{seen[$0]+=1} END{print length(seen)}') # find the number of unique elements of GROUPING to which elements of JOIN belong; this is the maximum number of extra groups in the model
n_new_haplotypes=$(awk 'BEGIN{if ('$max_new_haplotypes'-'$(echo $(((${#JOIN[@]}>0)) && echo ${#JOIN[@]} || echo 1))'>=0) {print '$max_new_haplotypes'-('$max_new_haplotypes'-'$(echo $(((${#JOIN[@]}>0)) && echo ${#JOIN[@]} || echo 1))')} else {print '$max_new_haplotypes'} }') # total number new haplotypes to add to model (by joining to haplotypes in JOIN)
(($n_new_haplotypes <= 0)) && (>&2 echo "No new haplotypes added."; exit 1)
if ((${#group_array[@]}==0))
then
  group_array=() # for defining the null group but not counting it
fi
declare -p group_array

if [ -z $CUTOFF ]
then
  CUTOFF=0.01 # LRT p-value cutoff; stop when all p exceed
fi

k=0 # initiate rounds variable
pval=0 # initiate p value variable

test -f $INCLUDED_HAPLOTYPES && echo cp $INCLUDED_HAPLOTYPES ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.txt || (>&2 echo "Included haplotypes file not found"; exit 1)
test -f $INCLUDED_HAPLOTYPES && cp $INCLUDED_HAPLOTYPES ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.txt && printf "\n"

while (( $(awk 'BEGIN{ if ('$pval'+0<'$CUTOFF') {print 1} else {print 0} }')==1 ))
do
  k=$((k+1))
  test -f ${NAME/%/.}Conditional_haplotype_effects.$(eval "for i in {$((${#group_array[@]}+1))..$((${#group_array[@]}+$n_new_haplotypes))}; do printf \"%s,\" h\$i; done | sed 's/,\$//g'").txt && echo rm ${NAME/%/.}Conditional_haplotype_effects.$(eval "for i in {$((${#group_array[@]}+1))..$((${#group_array[@]}+$n_new_haplotypes))}; do printf \"%s,\" h\$i; done | sed 's/,\$//g'").txt
  test -f ${NAME/%/.}Conditional_haplotype_effects.$(eval "for i in {$((${#group_array[@]}+1))..$((${#group_array[@]}+$n_new_haplotypes))}; do printf \"%s,\" h\$i; done | sed 's/,\$//g'").txt && rm ${NAME/%/.}Conditional_haplotype_effects.$(eval "for i in {$((${#group_array[@]}+1))..$((${#group_array[@]}+$n_new_haplotypes))}; do printf \"%s,\" h\$i; done | sed 's/,\$//g'").txt
  echo touch ${NAME/%/.}Conditional_haplotype_effects.$(eval "for i in {$((${#group_array[@]}+1))..$((${#group_array[@]}+$n_new_haplotypes))}; do printf \"%s,\" h\$i; done | sed 's/,\$//g'").txt
  touch ${NAME/%/.}Conditional_haplotype_effects.$(eval "for i in {$((${#group_array[@]}+1))..$((${#group_array[@]}+$n_new_haplotypes))}; do printf \"%s,\" h\$i; done | sed 's/,\$//g'").txt && printf "\n"

  N_JOBS=0
  for file in ${NAME/%/.}included_haplotype_0.*
  do
    N_JOBS=$((N_JOBS+1)) # file count
  done
  x=$( ((${#JOIN[@]}>0)) && echo $((${#JOIN[@]}-1)) || echo ${#JOIN[@]}) # total number of included haplotypes, starting with 0

  echo bsub \-P SJLIFE \-J \"myJob[1-$N_JOBS]\" \-oo dbgap28544_haplotype_model9_fit.v2.${NAME/%/.}round_${k}.%I.out \-eo dbgap28544_haplotype_model9_fit.v2.${NAME/%/.}round_${k}.%I.err \-R \"rusage[mem=4000]\" \-R \"order[!ut]\" \-R \"select[ut \< 0.9]\" \-q standard \"sh ${HOME_DIR/%/\/}dbgap28544_haplotype_model9_fit.v2.sh \\\"${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.txt\\\" \\\"$(eval "echo ${NAME/%/.}included_haplotype_{0..${x}}.new_allele_counts.job\\\$LSB_JOBINDEX.txt | tr \" \" \",\"")\\\" \\\"$PHENOTYPES\\\" \\\"$3\\\" \\\"$4\\\" job\$LSB_JOBINDEX ${DIRECTORY/%/\/}${NAME/%/.}Conditional_haplotype_effects.$(eval "for i in {$((${#group_array[@]}+1))..$((${#group_array[@]}+$n_new_haplotypes))}; do printf \"%s,\" h\$i; done | sed 's/,\$//g'").txt $DIRECTORY $HOME_DIR\"
  bsub -P SJLIFE -J "myJob[1-$N_JOBS]" -oo dbgap28544_haplotype_model9_fit.v2.${NAME/%/.}round_${k}.%I.out -eo dbgap28544_haplotype_model9_fit.v2.${NAME/%/.}round_${k}.%I.err -R "rusage[mem=4000]" -R "order[!ut]" -R "select[ut < 0.9]" -q standard "sh ${HOME_DIR/%/\/}dbgap28544_haplotype_model9_fit.v2.sh \"${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.txt\" \"$(eval "echo ${NAME/%/.}included_haplotype_{0..${x}}.new_allele_counts.job\\\$LSB_JOBINDEX.txt | tr \" \" \",\"")\" \"$PHENOTYPES\" \"$3\" \"$4\" job\$LSB_JOBINDEX ${DIRECTORY/%/\/}${NAME/%/.}Conditional_haplotype_effects.$(eval "for i in {$((${#group_array[@]}+1))..$((${#group_array[@]}+$n_new_haplotypes))}; do printf \"%s,\" h\$i; done | sed 's/,\$//g'").txt $DIRECTORY $HOME_DIR" && printf "\n"


  echo Wait until all jobs done:
  job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}')) # get job name from 7th field (or 6th if no exectution host yet) in all non-header rows of bjobs
  declare -p job_array
  while (( ${#job_array[@]}>0 ))
  do
    sleep 60
    job_array=($(bjobs 2> /dev/null | awk '($7 ~ /myJob/){print $7} ($6 ~ /myJob/){print $6}'))
    # declare -p job_array
  done
  printf "\n"

  pval=$(awk 'BEGIN{pval=1} ($5 ~ /[0-9]+/ && $5+0<pval){pval=$5} END{print pval}' ${NAME/%/.}Conditional_haplotype_effects.$(eval "for i in {$((${#group_array[@]}+1))..$((${#group_array[@]}+$n_new_haplotypes))}; do printf \"%s,\" h\$i; done | sed 's/,\$//g'").txt) # LRT p-value from kth round
  test_haplotype=$(awk 'BEGIN{pval=1} ($5 ~ /[0-9]+/ && $5+0<pval){pval=$5; haplotype=$1} END{print haplotype}' ${NAME/%/.}Conditional_haplotype_effects.$(eval "for i in {$((${#group_array[@]}+1))..$((${#group_array[@]}+$n_new_haplotypes))}; do printf \"%s,\" h\$i; done | sed 's/,\$//g'").txt) # haplotype corresponding to minimum p-value
  awk 'BEGIN{pval=1} ($5 ~ /[0-9]+/ && $5+0<pval){pval=$5; line=$0} END{printf "%s\t%s\n","'$(eval "for i in {$((${#group_array[@]}+1))..$((${#group_array[@]}+$n_new_haplotypes))}; do printf \"%s,\" h\$i; done | sed 's/,\$//g'")'",line}' ${NAME/%/.}Conditional_haplotype_effects.$(eval "for i in {$((${#group_array[@]}+1))..$((${#group_array[@]}+$n_new_haplotypes))}; do printf \"%s,\" h\$i; done | sed 's/,\$//g'").txt >> $OUTPUT # add line corresponding to minimum p-value to output file with a label
  if ((${#JOIN[@]}>0))
  then
    for ((i=0;i<${#JOIN[@]};i++)) # loop through included haplotypes which test haplotypes were joined to
    do
      included_haplotype=$(awk 'NR==1{print $'$((i+2))'}' ${INCLUDED_HAPLOTYPES})
      for file in ${NAME/%/.}*new_allele_counts*
      do
        col=$(awk 'BEGIN{n=split("'${test_haplotype}$(test ! -z $included_haplotype && echo ",")${included_haplotype}'",array,","); for (i=1;i<=n;i++) {hap[array[i]]} } NR==1{for (i=2;i<=NF;i++) {delete array;  m=split($i,array,","); if (m==n) {for (j=1;j<=m;j++) {if (array[j] in hap==0) {break} else {if (j==m) {print i} } } } } }' $file) # column corresponding to the smallest p-value in $file; will be empty if haplotype is not found
        if [ ! -z $col ]
        then
          echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{seen[\$1]=\$${col}\; next} \(\$1 in seen\){print \$0,seen[\$1]}\' $file ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.txt \> ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.tmp
          awk 'BEGIN{OFS="\t"} NR==FNR{seen[$1]=$'$col'; next} ($1 in seen){print $0,seen[$1]}' $file ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.txt > ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.tmp && printf "\n" # concatenate the column corresponding to the new haplotype with the lowest p-value, repeat for all included haplotypes

          # updated file of included haplotypes
          test -f ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.tmp && echo mv ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.tmp ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.txt
          test -f ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.tmp && mv ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.tmp ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.txt && printf "\n"

          unset col

          break # go to next included haplotype
        fi
      done
    done
  else # if the test haplotypes are not joined to any included haplotypes
    for file in ${NAME/%/.}*new_allele_counts*
    do
      col=$(awk 'BEGIN{n=split("'${test_haplotype}'",array,","); for (i=1;i<=n;i++) {hap[array[i]]} } NR==1{for (i=2;i<=NF;i++) {delete array; m=split($i,array,","); if (m==n) {for (j=1;j<=m;j++) {if (array[j] in hap==0) {break} else {if (j==m) {print i} } } } } }' $file) # column corresponding to the smallest p-value in $file; will be empty if haplotype is not found
      if [ ! -z $col ]
      then
        echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{seen[\$1]=\$$col\; next} \(\$1 in seen\){print \$0,seen[\$1]}\' $file ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.txt \> ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.tmp
        awk 'BEGIN{OFS="\t"} NR==FNR{seen[$1]=$'$col'; next} ($1 in seen){print $0,seen[$1]}' $file ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.txt > ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.tmp && printf "\n" # concatenate the column corresponding to the new haplotype with the lowest p-value, repeat for all included haplotypes

        # updated file of included haplotypes
        test -f ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.tmp && echo mv ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.tmp ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.txt
        test -f ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.tmp && mv ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.tmp ${INCLUDED_HAPLOTYPES%.*}.included_haplotypes.txt && printf "\n"

        unset col

        break # done with this round
      fi
    done
  fi

  # for ((i=1;i<=$n_new_haplotypes;i++)) # add one new group for each included_haplotype group to which the test_haplotype is joined
  # do
  #   if ((${#group_array[@]}>0))
  #   then
  #     group_array+=($((${group_array[-1]}+1)))
  #   else
  #     group_array+=(1)
  #   fi
  # done


  ll=$((${#group_array[@]}+1)) # lower limit of joined haplotypes in list
  ul=$ll # initial upper limit
  join_group=$(awk -v var=${JOIN[0]} 'BEGIN{n=split("'$(echo ${GROUPING[@]} | sed 's/\s\+/,/g')'",array,","); if (n>0) {for (i=1;i<=n;i++) {a=gensub("([0-9]*)-[0-9]*","\\1","g",array[i]); b=gensub("[0-9]*-([0-9]*)","\\1","g",array[i]); if (var+0>=a && var+0<=b) {print array[i]} } } else {print 0} }') # find position of first JOIN index in GROUPING; returns 0 if GROUPING is empty, which can only happen if JOIN is empty
  if ((${#JOIN[@]}>1)) # JOIN to more than one haplotype in GROUPING
  then
    for ((i=1;i<${#JOIN[@]};i++))
    do
      if [ $(awk -v var=${JOIN[i]} 'BEGIN{n=split("'$(echo ${GROUPING[@]} | sed 's/\s\+/,/g')'",array,","); for (i=1;i<=n;i++) {a=gensub("([0-9]*)-[0-9]*","\\1","g",array[i]); b=gensub("[0-9]*-([0-9]*)","\\1","g",array[i]); if (var>=a && var<=b) {print array[i]} } }') == "$join_group" ] # check if current JOIN is in the same position of GROUPING as the last
      then
        ul=$((ul+1))
      else
        group_array+=("${ll}-${ul}") # complete new group
        ll=$((ul+1)) # start of next group
      fi
      join_group=$(awk -v var=${JOIN[i]} 'BEGIN{n=split("'$(echo ${GROUPING[@]} | sed 's/\s\+/,/g')'",array,","); for (i=1;i<=n;i++) {a=gensub("([0-9]*)-[0-9]*","\\1","g",array[i]); b=gensub("[0-9]*-([0-9]*)","\\1","g",array[i]); if (var>=a && var<=b) {print array[i]} } }')
    done
  else
    group_array+=("${ll}") # only new group is the joinging to the empty group
  fi
  declare -p group_array
  printf "\n"
done
