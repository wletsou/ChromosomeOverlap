#! /bin/bash
set -e

PHENOTYPES=$1 # phenotypes file for getting age and affected status
POPULATION=$2 # haplotype_estimates file for combined population
TEST_HAPLOTYPES=$3 # colon-separated list of comma-separated lists rsid1_allele=[0,1] to be joined to INCLUDED_HAPLOTYPES
JOIN=$4 # comma-separated list of which INCLUDED_HAPLOTYPES to join to (start with "1", correction will be done for 0-based indexing);
INCLUDED_HAPLOTYPES=$5 # colon-separated list or file of comma-separated lists rsid1_allele=[0,1]
GROUPING=$6 # comma-separated list 1,2-4,5,6-8,... of how the INCLUDED_HAPLOTYPES should be grouped in the model; single-grouping if empty
OUTPUT=$7 # optional name of output file
NAME=$8
DIRECTORY=$9
HOME_DIR=${10}

module load R/3.6.1

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY

test -f ${PHENOTYPES} || (>&2 echo "Phenotypes file not found."; exit 1)

BP_RANGE=($(echo $BP_RANGE | perl -pne 's/([0-9]+)[,]*/$1 /g'))

POPULATION=($(echo $POPULATION | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))

SNP_LIST=$(awk 'NR==1{for (j=2;j<=NF;j++) {printf "%s%s",$j,(j<NF?",":"")} }' ${POPULATION[0]}) # get SNPs in order from first haplotype_estimates file
char=$(($(echo $(($(echo $SNP_LIST | tr -cd "," | wc -c)+1)) | wc -c)-1)) # number of characters in the total number snps

if [ -f $TEST_HAPLOTYPES ] && [ ! -z $TEST_HAPLOTYPES ] # list can be supplied from a file
then
  TEST_HAPLOTYPES=$(cat $TEST_HAPLOTYPES)
fi
TEST_HAPLOTYPES=($(echo $TEST_HAPLOTYPES | perl -pne 's/[:]/ /g'))
declare -p TEST_HAPLOTYPES
printf "\n"

if [ -f $INCLUDED_HAPLOTYPES ] && [ ! -z $INCLUDED_HAPLOTYPES ] # list can be supplied from a file
then
  INCLUDED_HAPLOTYPES=$(cat $INCLUDED_HAPLOTYPES)
fi
INCLUDED_HAPLOTYPES=($(echo $INCLUDED_HAPLOTYPES | perl -pne 's/[:]/ /g'))
declare -p INCLUDED_HAPLOTYPES
n_included_haplotypes=${#INCLUDED_HAPLOTYPES[@]}
printf "\n"

if [ ! -z $JOIN ]
then
  JOIN=($(echo $JOIN | perl -pne 's/[,]/ /g')) # array of haplotypes to join test haplotype to
else
  # JOIN=($(eval "echo {1..$( ((${#INCLUDED_HAPLOTYPES[@]}>0)) && echo ${#INCLUDED_HAPLOTYPES[@]} || echo 1 )}")) # indicates joining of test haplotype to all included haplotypes
  JOIN=() # indicates joining to the empty haplotype
fi


max_haplotype=$(awk 'BEGIN{max=0; n=split("'$(echo ${GROUPING} | sed 's/[-]/,/g')'",array,","); for (i=1;i<=n;i++) {a=gensub("([0-9]*)-[0-9]*","\\1","g",array[i]); if (a>max) {max=a}; b=gensub("[0-9]*-([0-9]*)","\\1","g",array[i]); if (b>max) {max=b} }; print max}' ) # maximum variable name in GROUPING array
if [ ! -z $n_included_haplotypes ] && [ -z $GROUPING ]
then
  GROUPING=($(echo $GROUPING | sed 's/[,]/ /g')$(test ! -z $GROUPING && echo " " || echo "")$((($max_haplotype+1<=$n_included_haplotypes)) && eval "echo {$(($max_haplotype+1))..$( (($n_included_haplotypes>$max_haplotype)) && echo $n_included_haplotypes || echo $(($max_haplotype+1)) )}" || echo "")) # grouping into single-haplotype groups of included haplotypes which have not already been grouped
fi
declare -p GROUPING
groups=$(echo "${GROUPING[*]}" | sed 's/\s\+/,/g')

max_new_haplotypes=$(awk 'BEGIN{n=split("'$(echo ${GROUPING[@]} | sed 's/\s\+/,/g')'",array,","); for (i=1;i<=n;i++) {GROUPING[array[i]]}; delete array; m=split("'$(echo ${JOIN[@]} | sed 's/\s\+/,/g')'",array,","); for (i=1;i<=m;i++) {JOIN[array[i]]}; if (m>0) {for (i in JOIN) {for (j in GROUPING) {a=gensub("([0-9]*)-[0-9]*","\\1","g",j); b=gensub("[0-9]*-([0-9]*)","\\1","g",j); if (i>=a && i<=b) {print j} } } } else {print 0} }' | awk '{seen[$0]+=1} END{print length(seen)}') # find the number of unique elements of GROUPING to which elements of JOIN belong; this is the maximum number of extra groups in the model
n_new_haplotypes=$(awk 'BEGIN{if ('$max_new_haplotypes'-'$(echo $(((${#JOIN[@]}>0)) && echo ${#JOIN[@]} || echo 1))'>=0) {print '$max_new_haplotypes'-('$max_new_haplotypes'-'$(echo $(((${#JOIN[@]}>0)) && echo ${#JOIN[@]} || echo 1))')} else {print '$max_new_haplotypes'} }') # number of new variables to add to model (by joining to haplotypes in JOIN) for each TEST_HAPLOTYPE
(($n_new_haplotypes <= 0)) && (>&2 echo "No new haplotypes added."; exit 1)

if [ ! -z $NAME ]
then
  NAME=${NAME}"."
fi

# correction for 0-based indexing
for ((i=0;i<${#JOIN[@]};i++))
do
  JOIN[i]=$((JOIN[i]-1))
done
declare -p JOIN

for ((i=0;i<${#INCLUDED_HAPLOTYPES[@]};i++))
do
  # get haplotype counts (0 or 1) for each chromosome in POPULATION[i]
  echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${INCLUDED_HAPLOTYPES[i]}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; haplotype_snps[a]=b} } NR==1{ for \(i = 2\;i\<=NF\;i++\) {if \(\$i in haplotype_snps\) {hap_col[i]=\$i} }\; print \"sid\",\"${INCLUDED_HAPLOTYPES[i]}\"} NR\>1{hap_count=1\; for \(j in hap_col\) {if \(\$j!=haplotype_snps[hap_col[j]]\) {hap_count=0}\; if \(\$j !~ /[0-9]/\) {hap_count=\"NA\"} }\; print \$1,hap_count }\' ${POPULATION} \> ${NAME}allele_counts.haplotype_${i}.txt
  awk 'BEGIN{OFS="\t"; n=split("'${INCLUDED_HAPLOTYPES[i]}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } NR==1{ for (i = 2;i<=NF;i++) {if ($i in haplotype_snps) {hap_col[i]=$i} }; print "sid","'${INCLUDED_HAPLOTYPES[i]}'"} NR>1{seen[$1]+=1;hap_count=1; for (j in hap_col) {if ($j!=haplotype_snps[hap_col[j]]) {hap_count=0}; if ($j !~ /[0-9]/) {hap_count="NA"} }; printf "%s_%s\t%s\n",$1,seen[$1],hap_count }' ${POPULATION} > ${NAME}allele_counts.haplotype_${i}.txt
  printf "\n"

  if ((i==0))
  then
    echo mv ${NAME}allele_counts.haplotype_${i}.txt ${NAME}allele_counts.txt && printf "\n"
    mv ${NAME}allele_counts.haplotype_${i}.txt ${NAME}allele_counts.txt
  else
    # string of found haplotypes for each chromosome
    echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{hap[\$1]=\$2\; next} \(\$1 in hap\){\$2=sprintf\(\"%s%s\"\,\$2\,\(FNR==1?\"+\":\"\"\),hap[\$1]\)\; print \$0}\' ${NAME}allele_counts.haplotype_${i}.txt ${NAME}allele_counts.txt \> ${NAME}allele_counts.tmp
    awk 'BEGIN{OFS="\t"} NR==FNR{hap[$1]=$2; next} ($1 in hap){ $2=sprintf("%s%s%s",$2,(FNR==1?"+":""),hap[$1]); print $0 }' ${NAME}allele_counts.haplotype_${i}.txt ${NAME}allele_counts.txt > ${NAME}allele_counts.tmp
    printf "\n"

    test -f ${NAME}allele_counts.tmp && echo mv ${NAME}allele_counts.tmp ${NAME}allele_counts.txt
    test -f ${NAME}allele_counts.tmp && mv ${NAME}allele_counts.tmp ${NAME}allele_counts.txt && printf "\n"

    test -f ${NAME}allele_counts.haplotype_${i}.txt && echo rm ${NAME}allele_counts.haplotype_${i}.txt
    test -f ${NAME}allele_counts.haplotype_${i}.txt && rm ${NAME}allele_counts.haplotype_${i}.txt && printf "\n"
  fi
done

included_array=("${INCLUDED_HAPLOTYPES[@]}") # baseline list of included haplotypes already counted https://stackoverflow.com/questions/19417015/how-to-copy-an-array-in-bash
for ((i=0;i<${#TEST_HAPLOTYPES[@]};i++)) # join multiple test haplotypes to included halotypes in the indicated join group
do
  ll=$((${#included_array[@]}+1)) # lower limit of joined haplotypes in list
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
        ll=$((ul+1)) # start of next group
      fi
      join_group=$(awk -v var=${JOIN[i]} 'BEGIN{n=split("'$(echo ${GROUPING[@]} | sed 's/\s\+/,/g')'",array,","); for (i=1;i<=n;i++) {a=gensub("([0-9]*)-[0-9]*","\\1","g",array[i]); b=gensub("[0-9]*-([0-9]*)","\\1","g",array[i]); if (var>=a && var<=b) {print array[i]} } }')
    done
  else
    groups="${groups}$(((${#groups}>0)) && echo "," || echo "")${ll}" # only new group is the joinging to the empty group
  fi

  if ((${#JOIN[@]}>0))
  then
    for ((j=0;j<${#JOIN[@]};j++))
    do
      pattern=${INCLUDED_HAPLOTYPES[${JOIN[j]}]}$( ((${#INCLUDED_HAPLOTYPES[${JOIN[j]}]}>0)) && echo "," || echo "")${TEST_HAPLOTYPES[i]} # join test haplotype to select included haplotypes
      pattern=$(echo $pattern | awk 'BEGIN{m=split("'$SNP_LIST'",snp_list,",")} {n=split($0,alleles,","); for (i=1;i<=n;i++) {for (j=1;j<=m;j++) {if (alleles[i] ~ snp_list[j]) {out[sprintf("%0'${char}'d",j)]=alleles[i]} } } } END{n=asorti(out,dest); for (i=1;i<=n;i++) {printf "%s%s",out[dest[i]],(i<n?",":"")} }') # sort pattern in haplotype order

      included_array+=(${pattern}) # add joined haplotype to list of included haplotypes
    done
  else
    pattern=${TEST_HAPLOTYPES[i]} # no joining needed
    pattern=$(echo $pattern | awk 'BEGIN{m=split("'$SNP_LIST'",snp_list,",")} {n=split($0,alleles,","); for (i=1;i<=n;i++) {for (j=1;j<=m;j++) {if (alleles[i] ~ snp_list[j]) {out[sprintf("%0'${char}'d",j)]=alleles[i]} } } } END{n=asorti(out,dest); for (i=1;i<=n;i++) {printf "%s%s",out[dest[i]],(i<n?",":"")} }') # sort pattern in haplotype order

    included_array+=(${pattern}) # add haplotype to list of included haplotypes
  fi
  included=$(echo ${included_array[*]} | perl -pne 's/ /:/g') # print array as a colon-separated list
done

for ((i=${#INCLUDED_HAPLOTYPES[@]};i<${#included_array[@]};i++)) # loop over newly added haplotypes
do
  # get haplotype counts (0 or 1) for each chromosome in POPULATION[i]
  if ((i==${#INCLUDED_HAPLOTYPES[@]}))
  then
    (test -f ${NAME}allele_counts.txt && echo cp ${NAME}allele_counts.txt ${NAME}new_allele_counts.txt || echo awk \'BEGIN{OFS=\"\\t\"} NR==1{ print \"sid\" } NR\>1{seen[\$1]+=1\; printf \"%s_%s\\n\",\$1,seen[\$1] }\' ${POPULATION} \> ${NAME}new_allele_counts.txt) && printf "\n"
    test -f ${NAME}allele_counts.txt && cp ${NAME}allele_counts.txt ${NAME}new_allele_counts.txt || awk 'BEGIN{OFS="\t"} NR==1{ print "sid" } NR>1{seen[$1]+=1; printf "%s_%s\n",$1,seen[$1] }' ${POPULATION} > ${NAME}new_allele_counts.txt && printf "\n"
  fi

  echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${included_array[i]}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; haplotype_snps[a]=b} } NR==1{ for \(i = 2\;i\<=NF\;i++\) {if \(\$i in haplotype_snps\) {hap_col[i]=\$i} }\; print \"sid\",\"${included_array[j]}\"} NR\>1{hap_count=1\; for \(j in hap_col\) {if \(\$j!=haplotype_snps[hap_col[j]]\) {hap_count=0}\; if \(\$j !~ /[0-9]/\) {hap_count=\"NA\"} }\; print \$1,hap_count }\' ${POPULATION} \> ${NAME}new_allele_counts.haplotype_${i}.txt
  awk 'BEGIN{OFS="\t"; n=split("'${included_array[i]}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } NR==1{ for (i = 2;i<=NF;i++) {if ($i in haplotype_snps) {hap_col[i]=$i} }; print "sid","'${included_array[j]}'"} NR>1{seen[$1]+=1;hap_count=1; for (j in hap_col) {if ($j!=haplotype_snps[hap_col[j]]) {hap_count=0}; if ($j !~ /[0-9]/) {hap_count="NA"} }; printf "%s_%s\t%s\n",$1,seen[$1],hap_count }' ${POPULATION} > ${NAME}new_allele_counts.haplotype_${i}.txt && printf "\n"

  # string of found haplotypes for each chromosome
  echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{hap[\$1]=\$2\; next} \(\$1 in hap\){\$2=sprintf\(\"%s%s%s\"\,\$2\,\(FNR==1?\"+\":\"\"\),hap[\$1]\)\; print \$0}\' ${NAME}new_allele_counts.haplotype_${i}.txt ${NAME}new_allele_counts.txt \> ${NAME}new_allele_counts.tmp
  awk 'BEGIN{OFS="\t"} NR==FNR{hap[$1]=$2; next} ($1 in hap){ $2=sprintf("%s%s%s",$2,(FNR==1?"+":""),hap[$1]); print $0 }' ${NAME}new_allele_counts.haplotype_${i}.txt ${NAME}new_allele_counts.txt > ${NAME}new_allele_counts.tmp && printf "\n"

  test -f ${NAME}new_allele_counts.tmp && echo mv ${NAME}new_allele_counts.tmp ${NAME}new_allele_counts.txt
  test -f ${NAME}new_allele_counts.tmp && mv ${NAME}new_allele_counts.tmp ${NAME}new_allele_counts.txt && printf "\n"

  test -f ${NAME}new_allele_counts.haplotype_${i}.txt && echo rm ${NAME}new_allele_counts.haplotype_${i}.txt
  test -f ${NAME}new_allele_counts.haplotype_${i}.txt && rm ${NAME}new_allele_counts.haplotype_${i}.txt && printf "\n"

done
# add haplotype counts for each chromosome, including enough columns for all observed haplotypes
echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{ if \(NR\>1\) {seen[\$2]}\; next } \(NR!=FNR\){ if \(FNR==1\) {printf \"%s\\t\",\"sid\"\; m=asorti\(seen,haps\)\; for \(j=1\;j\<=m\;j++\) {printf \"%s%s\",haps[j],\(j\<m?\"\\t\":\"\\n\"\)\; col[haps[j]]=j+1} } else {b=gensub\(\"\([0-9]*\)_.*\",\"\\\\1\",\"g\",\$1\)\; chr[b]+=1\; array[sprintf\(\"%s\,%s\"\,b\,\$2\)]+=1\; if \(chr[b]==2\) {\$1=b\; for \(i=2\;i\<=m+1\;i++\) {\$i=0}\; for \(j in array\) {if \(j \~ b\) {c=gensub\(\".*[,]\([0-9]*\)\",\"\\\\1\",\"g\",j\)\; \$col[c]=array[j]\; delete array[j]} }\; print \$0 } } }\' ${NAME}new_allele_counts.txt ${NAME}new_allele_counts.txt \> ${NAME}new_allele_counts.tmp
awk 'BEGIN{OFS="\t"} NR==FNR{ if (NR>1) {seen[$2]}; next } (NR!=FNR){ if (FNR==1) {printf "%s\t","sid"; m=asorti(seen,haps); for (j=1;j<=m;j++) {printf "%s%s",haps[j],(j<m?"\t":"\n"); col[haps[j]]=j+1} } else {b=gensub("([0-9]*)_.*","\\1","g",$1); chr[b]+=1; array[sprintf("%s,%s",b,$2)]+=1; if (chr[b]==2) {$1=b; for (i=2;i<=m+1;i++) {$i=0}; for (j in array) {if (j ~ b) {c=gensub(".*[,]([0-9]*)","\\1","g",j); $col[c]=array[j]; delete array[j]} }; print $0 } } }' ${NAME}new_allele_counts.txt ${NAME}new_allele_counts.txt > ${NAME}new_allele_counts.tmp && printf "\n"

test -f ${NAME}new_allele_counts.tmp && echo mv ${NAME}new_allele_counts.tmp ${NAME}new_allele_counts.txt
test -f ${NAME}new_allele_counts.tmp && mv ${NAME}new_allele_counts.tmp ${NAME}new_allele_counts.txt && printf "\n"

# get age and affected status
# echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{age[\$1]=\$4\; next} \(\$1 in age\){if \(FNR==1\) {printf \"%s\\t%s\\t%s\\n\",\$0,\"age\",\"affected\"} if \(FNR \> 1\) {print \$0,age[\$1],$((${#POPULATION[@]}-i-1))} }\' ${PHENOTYPES} ${NAME}new_allele_counts.txt \> ${NAME}new_allele_counts.tmp
# awk 'BEGIN{OFS="\t"} NR==FNR{age[$1]=$4; next} ($1 in age){if (FNR==1) {printf "%s\t%s\t%s\n",$0,"age","affected"} if (FNR > 1) {print $0,age[$1],'$((${#POPULATION[@]}-i-1))'} }' ${PHENOTYPES} ${NAME}new_allele_counts.txt > ${NAME}new_allele_counts.tmp
# printf "\n" # total haplotype counts for each individual; get age and affected status

max_pcs=10 # maximum number of pcs to include in model
# max_pcs=0
# get age and affected status and first 10 pcs
echo awk \'BEGIN{OFS=\"\\t\"} NR==1{for \(i=1\;i\<=NF\;i++\) {if \(\$i~\"pc\"\) {b=gensub\(\"pc\([0-9]*\)\",\"\\\\1\",\"g\",\$i\)\; pc[b,\$1]=\$i\; col_pc[b]=i}\; if \(\$i==\"ageonset\"\) {col_age=i\; age[\$1]=\$i}\; if \(\$i==\"BCa\"\) {col_diag=i\; diag[\$1]=\$i} }\; next } \(NR==FNR \&\& NR\>1\){age[\$1]=\$col_age\; diag[\$1]=\$col_diag\; for \(i=1\;i\<=$max_pcs\;i++\) \{ if \(i in col_pc\) {pc[i\,\$1]=\$col_pc[i]\} \}\; next} \(\$1 in age\){if \(FNR==1\) {printf \"%s\",\$0\; for \(i=1\;i\<=$max_pcs\;i++\) {printf \"\\tpc%02d\",i}\; printf \"\\t%s\\t%s\\n\",\"age\",\"affected\" } if \(FNR \> 1\) {printf \"%s\",\$0\; for \(i=1\;i\<=$max_pcs\;i++\) {printf \"\\t%0.5f\",pc[i,\$1]}\; printf \"\\t%s\\t%s\\n\",age[\$1],diag[\$1]} }\' $PHENOTYPES ${NAME}new_allele_counts.txt \> ${NAME}new_allele_counts.tmp
awk 'BEGIN{OFS="\t"} NR==1{for (i=1;i<=NF;i++) {if ($i~"pc") {b=gensub("pc([0-9]*)","\\1","g",$i); pc[b,$1]=$i; col_pc[b]=i}; if ($i=="ageonset") {col_age=i; age[$1]=$i}; if ($i=="BCa") {col_diag=i; diag[$1]=$i} }; next } (NR==FNR && NR>1){age[$1]=$col_age; diag[$1]=$col_diag; for (i=1;i<='$max_pcs';i++) { if (i in col_pc) {pc[i,$1]=$col_pc[i]} }; next} ($1 in age){if (FNR==1) {printf "%s",$0; for (i=1;i<='$max_pcs';i++) {printf "\tpc%02d",i}; printf "\t%s\t%s\n","age","affected" } if (FNR > 1) {printf "%s",$0; for (i=1;i<='$max_pcs';i++) {printf "\t%0.5f",pc[i,$1]}; printf "\t%s\t%s\n",age[$1],diag[$1]} }' $PHENOTYPES ${NAME}new_allele_counts.txt > ${NAME}new_allele_counts.tmp && printf "\n" # total haplotype counts for each individual; get age and affected status

test -f ${NAME}new_allele_counts.tmp && echo mv ${NAME}new_allele_counts.tmp ${NAME}new_allele_counts.txt
test -f ${NAME}new_allele_counts.tmp && mv ${NAME}new_allele_counts.tmp ${NAME}new_allele_counts.txt && printf "\n"

echo Testing haplotype$(((${#TEST_HAPLOTYPES[@]}>1)) && echo "s") ${TEST_HAPLOTYPES[*]}
echo Rscript ${HOME_DIR}/dbgap_haplotype_model3.R file=${DIRECTORY}/${NAME}new_allele_counts.txt groups="$groups" verbose=1 n=$((${#TEST_HAPLOTYPES[@]}*$n_new_haplotypes)) n_covs=$max_pcs
Rscript ${HOME_DIR}/dbgap_haplotype_model3.R file=${DIRECTORY}/${NAME}new_allele_counts.txt groups="$groups" verbose=1 n=$((${#TEST_HAPLOTYPES[@]}*$n_new_haplotypes)) n_covs=$max_pcs
if [ ! -z $OUTPUT ] # write conditional effect of TEST_HAPLOTYPE[0] (HR, p value) in supplied OUTPUT file
then
  output=$(Rscript ${HOME_DIR}/dbgap_haplotype_model3.R file=${DIRECTORY}/${NAME}new_allele_counts.txt groups="$groups" verbose=0 n=$((${#TEST_HAPLOTYPES[@]}*$n_new_haplotypes)) n_covs=$max_pcs)
  printf "\n"
  echo $output | awk 'BEGIN{OFS="\t"; n=split("'$(echo ${TEST_HAPLOTYPES[*]} | sed 's/\s\+/:/g')'",haps,":")}; {a=gensub(".*HR = ([NA0-9e.+-]*).*","\\1","g",$0); b=gensub(".*LRT p value = ([NA0-9e.+-]*).*","\\1","g",$0); printf "%s\t",'${group_no}'; for (i=n;i>0;i--) {printf "%s\t",haps[i]}; printf "%s\t%s\n",a,b}' >> $OUTPUT
fi
printf "\n"