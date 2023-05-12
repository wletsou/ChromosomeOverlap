#! /bin/bash
set -e

# for evaluating the newest TEST_HAPLOTYPE in a logistic model model for BCa

INCLUDED_HAPLOTYPES=$1 # counts of included_haplotypes for each chromosome (can be just a list of chromosomes if no included_haplotypes in model)
TEST_HAPLOTYPES=$2 # counts of each new haplotype to be joined to each included_haplotype (comma-separated list of files)
PHENOTYPES=$3 # phenotypes file for getting age ("ageonset"), principal components of ancestry ("pc1" to "pc10") and affected status ("BCa") from each subject ("sid")
JOIN=$4 # comma-separated list of which INCLUDED_HAPLOTYPES to join to (start with "1", correction will be done for 0-based indexing); use all if value is empty
GROUPING=$5 # comma-separated list 1,2-4,5,6-8,... of how the INCLUDED_HAPLOTYPES should be grouped in the model; single-grouping if empty
NAME=$6 # optional name appended to temporary files
OUTPUT=$7 # optional name of output file
DIRECTORY=$8
HOME_DIR=$9

module load R/3.6.1

if [ -z $HOME_DIR ];
then
  unset HOME_DIR
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY

if [ -z $NAME ]
then
  unset NAME
fi

test -z $PHENOTYPES && (>&2 echo "Phenotypes file not supplied"; exit 1)
test -f $PHENOTYPES || (>&2 echo "Phenotypes file $PHENOTYPES not found"; exit 1)

test_haplotypes_files=($(echo $TEST_HAPLOTYPES | perl -pne 's/[,]/ /g'))
n_test_haplotypes=$(awk 'NR==1{print NF-1}' ${test_haplotypes_files[0]}) # number of new haplotypes to be tested
declare -p test_haplotypes_files
for ((i=0;i<${#test_haplotypes_files[@]};i++))
do
  test -z ${test_haplotypes_files[i]} && (>&2 echo "Test haplotypes file ${test_haplotypes_files[i]} file not supplied"; exit 1)
  test -f ${test_haplotypes_files[i]} || (>&2 echo "Test haplotypes file ${test_haplotypes_files[i]} file not found"; exit 1)
done
printf "\n"

test -z ${INCLUDED_HAPLOTYPES} && (>&2 echo "Included haplotypes file ${INCLUDED_HAPLOTYPES} file not supplied"; exit 1)
test -f ${INCLUDED_HAPLOTYPES} || (>&2 echo "Included haplotypes file ${INCLUDED_HAPLOTYPES} file not found"; exit 1)
n_included_haplotypes=$(awk 'NR==1{print NF-1}' ${INCLUDED_HAPLOTYPES})

if [ ! -z $JOIN ]
then
  JOIN=($(echo $JOIN | perl -pne 's/[,]/ /g')) # array of haplotypes to join test haplotype to
else
  # JOIN=($(eval "echo {1..$( (($n_included_haplotypes>0)) && echo $n_included_haplotypes || echo 1 )}")) # indicates joining of test haplotype to all included haplotypes (one test_haplotype file per included haplotype)
  JOIN=() # indicates joining to the empty haplotype
fi
max_haplotype=$(awk 'BEGIN{max=0; n=split("'$(echo ${GROUPING} | sed 's/[-]/,/g')'",array,","); for (i=1;i<=n;i++) {a=gensub("([0-9]*)-[0-9]*","\\1","g",array[i]); if (a>max) {max=a}; b=gensub("[0-9]*-([0-9]*)","\\1","g",array[i]); if (b>max) {max=b} }; print max}' ) # maximum variable name in GROUPING array
if [ ! -z $n_included_haplotypes ]
then
  GROUPING=($(echo $GROUPING | sed 's/[,]/ /g')$(test ! -z $GROUPING && echo " " || echo "")$((($max_haplotype+1<=$n_included_haplotypes)) && eval "echo {$(($max_haplotype+1))..$( (($n_included_haplotypes>$max_haplotype)) && echo $n_included_haplotypes || echo $(($max_haplotype+1)) )}" || echo "")) # grouping into single-haplotype groups of included haplotypes which have not already been grouped
fi
declare -p GROUPING
groups=$(echo "${GROUPING[*]}" | sed 's/\s\+/,/g')
group_array=($(echo "${GROUPING[*]}"))
max_new_haplotypes=$(awk 'BEGIN{n=split("'$(echo ${GROUPING[@]} | sed 's/\s\+/,/g')'",array,","); for (i=1;i<=n;i++) {GROUPING[array[i]]}; delete array; m=split("'$(echo ${JOIN[@]} | sed 's/\s\+/,/g')'",array,","); for (i=1;i<=m;i++) {JOIN[array[i]]}; if (m>0) {for (i in JOIN) {for (j in GROUPING) {a=gensub("([0-9]*)-[0-9]*","\\1","g",j); b=gensub("[0-9]*-([0-9]*)","\\1","g",j); if (i>=a && i<=b) {print j} } } } else {print 0} }' | awk '{seen[$0]+=1} END{print length(seen)}') # find the number of unique elements of GROUPING to which elements of JOIN belong; this is the maximum number of extra groups in the model
n_new_haplotypes=$(awk 'BEGIN{if ('$max_new_haplotypes'-'$(echo $(((${#JOIN[@]}>0)) && echo ${#JOIN[@]} || echo 1))'>=0) {print '$max_new_haplotypes'-('$max_new_haplotypes'-'$(echo $(((${#JOIN[@]}>0)) && echo ${#JOIN[@]} || echo 1))')} else {print '$max_new_haplotypes'} }') # total number new haplotypes to add to model (by joining to haplotypes in JOIN)

# if ((${#group_array[@]}==0))
# then
#   group_array=(0) # for defining the null group but not counting it
# fi
# declare -p group_array
# printf "\n"

# for ((i=1;i<=$n_new_haplotypes;i++)) # add one new group for each included_haplotype group to which the test_haplotype is joined
# do
#   groups=${groups}$( test ! -z $groups && echo "," || echo "")$(($n_included_haplotypes+i)) # https://stackoverflow.com/questions/36977855/get-last-element-in-bash-array
#   group_array=($(echo $groups | perl -pne 's/[,]/ /g'))
# done
# declare -p group_array
# printf "\n"

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
      groups="${groups}$(((${#groups}>0)) && echo "," || echo "")${ll}-${ul}" # complete new group
      group_array+=("${ll}-${ul}")
      ll=$((ul+1)) # start of next group
    fi
    join_group=$(awk -v var=${JOIN[i]} 'BEGIN{n=split("'$(echo ${GROUPING[@]} | sed 's/\s\+/,/g')'",array,","); for (i=1;i<=n;i++) {a=gensub("([0-9]*)-[0-9]*","\\1","g",array[i]); b=gensub("[0-9]*-([0-9]*)","\\1","g",array[i]); if (var>=a && var<=b) {print array[i]} } }')
  done
else
  groups="${groups}$(((${#groups}>0)) && echo "," || echo "")${ll}" # only new group is the joinging to the empty group
  group_array+=("${ll}")
fi

for ((i=0;i<$n_test_haplotypes;i++))
do
  start=$(date +%s.%N)
  # echo awk \'BEGIN{OFS=\"\\t\"} {print \$1,\$$((i+2))}\' $INCLUDED_HAPLOTYPES \> ${INCLUDED_HAPLOTYPES%.*}.${NAME/%/.}haplotype_$((i+1)).txt
  # awk 'BEGIN{OFS="\t"} {print $1,$'$((i+2))'}' $INCLUDED_HAPLOTYPES > ${INCLUDED_HAPLOTYPES%.*}.${NAME/%/.}haplotype_$((i+1)).txt && printf "\n" # initialize haplotype count file

  for ((j=0;j<${#test_haplotypes_files[@]};j++))
  do
    if ((${#JOIN[@]}>0)) # empty columns of INCLUDED_HAPLOTYPES will be filled in as other haplotypes are included in the model
    then
      included_haplotype_name=$(awk 'NR==1{print $'$((j+2))'}' ${INCLUDED_HAPLOTYPES}) # only get included_haplotype names from the first ${#JOIN[@]} haplotypes to which test_haplotype are joined
    fi
    test_haplotype_name=$(awk 'NR==1{print $'$((i+2))'}' ${test_haplotypes_files[j]})
    hap_name_new=$(awk 'BEGIN{n=split("'$included_haplotype_name'",array1,","); for (i=1;i<=n;i++) {included[array1[i]]}; m=split("'$test_haplotype_name'",array2,","); for (i=1;i<=m;i++) {if (array2[i] in included==0) {printf "%s%s",(found>0?",":""),array2[i]; found+=1} } }') # alleles in test haplotype not part of included haplotype
    if ((j>0))
    then
      if [ "$hap_name_new" != "$hap_name" ] # see if column names are the same in each file
      then
        >&2 echo "test haplotype names do not match"; exit 1
      fi
    fi
    hap_name=$hap_name_new

  done

  echo Joining haplotype $hap_name:
  echo Rscript ${HOME_DIR/%/\/}dbgap_haplotype_model.v2.R -i ${INCLUDED_HAPLOTYPES} -t ${TEST_HAPLOTYPES} -p $PHENOTYPES -k $((i+1)) -g "$groups" -v 1 -n $((${#group_array[@]}-${#JOIN[@]}))
  if [ -z $OUTPUT ]
  then
    Rscript ${HOME_DIR/%/\/}dbgap_haplotype_model.v2.R -i ${INCLUDED_HAPLOTYPES} -t ${TEST_HAPLOTYPES} -p $PHENOTYPES -k $((i+1)) -g "$groups" -v 1 -n $((${#group_array[@]}-${#JOIN[@]}))
    printf "\n"
  fi
  if [ ! -z $OUTPUT ] # write conditional effect of TEST_HAPLOTYPE[0] (HR, p value) in supplied OUTPUT file
  then
    Rscript ${HOME_DIR/%/\/}dbgap_haplotype_model.v2.R -i ${INCLUDED_HAPLOTYPES} -t ${TEST_HAPLOTYPES} -p $PHENOTYPES -k $((i+1)) -g "$groups" -v 1 -n $((${#group_array[@]}-${#JOIN[@]})) > routput.${NAME/%/.}txt
    test -f routput.${NAME/%/.}txt && cat routput.${NAME/%/.}txt # print R output file
    printf "\n"
    cat routput.${NAME/%/.}txt | grep LRT | awk 'BEGIN{OFS="\t"} {a=gensub(".*OR = ([NA0-9e.+-]*.*),.*","\\1","g",$0); b=gensub(".*LRT p value = ([NA0-9e.+-]*).*","\\1","g",$0); c=gensub(".*frequency = ([NA0-9e.+-]*.*)$","\\1","g",$0); print a,b,c}' | sed 's/[^NA0-9e+-.]/ /g' | sed 's/\s\+/ /g' | awk '{n=split($0,array,/\s+/); printf "'$hap_name'\t"; for (i=1;i<=n;i++) {printf "%s%s",array[i],(i<n?"\t":"\n")} }' >> ${OUTPUT}

    test -f routput.${NAME/%/.}txt && echo rm routput.${NAME/%/.}txt
    test -f routput.${NAME/%/.}txt && rm routput.${NAME/%/.}txt && printf "\n"
  fi
  finish=$(date +%s.%N)
  duration=$(awk 'BEGIN{printf "%0.2f\n", '$finish'-'$start'}')
  echo Haplotype $((i+1)) completed in $duration second$( [ "$duration" == "1.00" ] && echo "" || echo "s").
  printf "\n"
done
