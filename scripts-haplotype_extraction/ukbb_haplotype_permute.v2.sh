#! /bin/bash
set -e

# for randomly permuting cases and controls such that the number of HAPLOTYPE hetero- and homozygotes is preserved

# sh /home/wletsou/scripts/ukbb_haplotype_permute.sh ukb.bca.pheno.unpermuted.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.69231642-69431642.txt rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0 20200116 ukbb_bca_cases,ukbb_bca_controls ukb.bca.pheno.txt

PHENOTYPES=$1 # case/controls status
POPULATION=$2 # combined haplotype_estimates file for population
HAPLOTYPE=$3 # comma-separated list of rsid_allele=[0/1],... to condition on
SEED=$4 # random seed
NAMES=$5 # comma-separated list of base names of output files for cases and controls
OUTPUT=$6 # phenotype output file (in working directory)
DIRECTORY=$7
HOME_DIR=$8

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

test -z ${PHENOTYPES} && (>&2 echo "Phenotypes file not supplied."; exit 1)
test -f ${PHENOTYPES} || (>&2 echo "Phenotypes file not found."; exit 1)
[ "$DIRECTORY" != "${PHENOTYPES%/*}" ] || (>&2 echo "Phenotypes file in working directory."; exit 1)

test -z ${POPULATION} && (>&2 echo "Haplotypes file not supplied."; exit 1)
test -f ${POPULATION} || (>&2 echo "Haplotypes file not found."; exit 1)

test -z ${HAPLOTYPE} && (>&2 echo "No conditioning haplotype supplied."; exit 1)

if [ -z $SEED ]
then
  SEED=20200116
fi

NAMES=($(echo $NAMES | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))
if ((${#NAMES[@]}!=2))
then
  NAMES=(cases controls)
fi

if [ -z $OUTPUT ]
then
  OUTPUT=phenotypes.txt
fi

test -f ${NAMES[0]}.indiv && echo rm ${NAMES[0]}.indiv
test -f ${NAMES[0]}.indiv && rm ${NAMES[0]}.indiv && printf "\n"
test -f ${NAMES[1]}.indiv && echo rm ${NAMES[1]}.indiv
test -f ${NAMES[1]}.indiv && rm ${NAMES[1]}.indiv && printf "\n"

echo touch ${NAMES[0]}.indiv
touch ${NAMES[0]}.indiv && printf "\n"
echo touch ${NAMES[1]}.indiv
touch ${NAMES[1]}.indiv && printf "\n"

# unpermuted cases and controls
echo awk \'\$3==1{print \$1,\$1}\' $PHENOTYPES \> cases.tmp
awk '$3==1{print $1,$1}' $PHENOTYPES > cases.tmp && printf "\n"
echo awk \'\$3==0{print \$1,\$1}\' $PHENOTYPES \> controls.tmp
awk '$3==0{print $1,$1}' $PHENOTYPES > controls.tmp && printf "\n"

for ((i=1;i<=2;i++))
do

  file=$(echo $( (($i==1)) && echo "heterozygotes" || ( (($i==2)) && echo "homozygotes") ).txt)

  echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${HAPLOTYPE}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; haplotype_snps[a]=b} } NR==1{ for \(i = 2\;i\<=NF\;i++\) {if \(\$i in haplotype_snps\) {hap_col[i]=\$i} } if \(length\(hap_col\)!=n\) {print \"Wrong number of SNPs found.\"\; exit 1} } NR\>1{hap_count=1\; for \(j in hap_col\) {if \(\$j!=haplotype_snps[hap_col[j]]\) {hap_count=0}\; if \(\$j !~ /[0-9]/\) {hap_count=\"NA\"} }\; seen[\$1]+=hap_count+0 } END{for \(i in seen\) {if \(seen[i]==${i}\) {print i,i} } }\' ${POPULATION[*]} \> $file
  awk 'BEGIN{OFS="\t"; n=split("'${HAPLOTYPE}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } NR==1{ for (i = 2;i<=NF;i++) {if ($i in haplotype_snps) {hap_col[i]=$i} } if (length(hap_col)!=n) {print "Wrong number of SNPs found."; exit 1} } NR>1{hap_count=1; for (j in hap_col) {if ($j!=haplotype_snps[hap_col[j]]) {hap_count=0}; if ($j !~ /[0-9]/) {hap_count="NA"} }; seen[$1]+=hap_count+0 } END{for (i in seen) {if (seen[i]=='$i') {print i,i} } }' ${POPULATION[*]} > $file && printf "\n" # find carriers of HAPLOTYPE and print sid,sid

  echo awk \'NR==FNR{seen[\$1]=\$3\; next} \(\$1 in seen \&\& seen[\$1]==1\){cases+=1} END{print cases}\' $PHENOTYPES $file
  n_samples=$(awk 'NR==FNR{seen[$1]=$3; next} ($1 in seen && seen[$1]==1){cases+=1} END{print cases}' $PHENOTYPES $file) && printf "\n" # number of hetero- or homozygote cases

  echo Rscript ${HOME_DIR/%//}UKBB_row_sample.R \-f $file \--no_header \-n $n_samples \-s $SEED \-o ${file%.*}.row_sample_${n_samples}.txt \-d $DIRECTORY
  Rscript ${HOME_DIR/%//}UKBB_row_sample.R -f $file --no_header -n $n_samples -s $SEED -o ${file%.*}.row_sample_${n_samples}.txt -d $DIRECTORY && printf "\n"
  # sample n_samples lines to be cases

  test -f ${DIRECTORY/%//}${file%.*}.row_sample_${n_samples}.txt && echo cat ${DIRECTORY/%//}${file%.*}.row_sample_${n_samples}.txt \>\> ${NAMES[0]}.indiv
  test -f ${DIRECTORY/%//}${file%.*}.row_sample_${n_samples}.txt && cat ${DIRECTORY/%//}${file%.*}.row_sample_${n_samples}.txt >> ${NAMES[0]}.indiv && printf "\n" # permuted cases

  echo awk \'NR==FNR{seen[\$1]\; next} \(\$1 in seen==0\){print \$0}\' ${NAMES[0]}.indiv $file \>\> ${NAMES[1]}.indiv
  awk 'NR==FNR{seen[$1]; next} ($1 in seen==0){print $0}' ${NAMES[0]}.indiv $file >> ${NAMES[1]}.indiv && printf "\n" # hetero- or homozygotes not in cases

  echo awk \'NR==FNR{seen[\$1]\; next} \(\$1 in seen==0\){print \$0}\' $file cases.tmp \> cases.1.tmp
  awk 'NR==FNR{seen[$1]; next} ($1 in seen==0){print $0}' $file cases.tmp > cases.1.tmp && printf "\n" # unpermuted cases who are not hetero- or homozygotes
  test -f cases.1.tmp && echo mv cases.1.tmp cases.tmp
  test -f cases.1.tmp && mv cases.1.tmp cases.tmp && printf "\n"

  echo awk \'NR==FNR{seen[\$1]\; next} \(\$1 in seen==0\){print \$0}\' $file controls.tmp \> controls.1.tmp
  awk 'NR==FNR{seen[$1]; next} ($1 in seen==0){print $0}' $file controls.tmp > controls.1.tmp && printf "\n" # unpermuted controls who are not hetero- or homozygotes
  test -f controls.1.tmp && echo mv controls.1.tmp controls.tmp
  test -f controls.1.tmp && mv controls.1.tmp controls.tmp && printf "\n"

  test -f $file && echo rm $file
  test -f $file && rm $file && printf "\n"
  test -f ${DIRECTORY/%//}${file%.*}.row_sample_${n_samples}.txt && echo rm ${DIRECTORY/%//}${file%.*}.row_sample_${n_samples}.txt
  test -f ${DIRECTORY/%//}${file%.*}.row_sample_${n_samples}.txt && rm ${DIRECTORY/%//}${file%.*}.row_sample_${n_samples}.txt && printf "\n"
done

test -f cases.tmp && echo cat cases.tmp \>\> ${NAMES[0]}.indiv
test -f cases.tmp && cat cases.tmp >> ${NAMES[0]}.indiv && printf "\n" # add unpermuted cases into cases file
test -f cases.tmp && echo rm cases.tmp
test -f cases.tmp && rm cases.tmp && printf "\n"

test -f controls.tmp && echo cat controls.tmp \>\> ${NAMES[1]}.indiv
test -f controls.tmp && cat controls.tmp >> ${NAMES[1]}.indiv && printf "\n" # add unpermuted cases into cases file
test -f controls.tmp && echo rm controls.tmp
test -f controls.tmp && rm controls.tmp && printf "\n"

test -f $OUTPUT && echo rm $OUTPUT
test -f $OUTPUT && rm $OUTPUT && printf "\n"

[ "$OUTPUT" != "$PHENOTYPES" ] || (>&2 echo "Output is the original phenotypes file."; exit 1)

echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{seen[\$1]\; next} \(\$1 in seen\){\$3=1}1\' ${NAMES[0]}.indiv $PHENOTYPES \> ${OUTPUT}.tmp
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$1]; next} ($1 in seen){$3=1}1' ${NAMES[0]}.indiv $PHENOTYPES > ${OUTPUT}.tmp && printf "\n" # replace observed cases with permuted cases

echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{seen[\$1]\; next} \(\$1 in seen\){\$3=0}1\' ${NAMES[1]}.indiv ${OUTPUT}.tmp \>\> ${DIRECTORY/%//}${OUTPUT}
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$1]; next} ($1 in seen){$3=0}1' ${NAMES[1]}.indiv ${OUTPUT}.tmp >> ${DIRECTORY/%//}${OUTPUT} && printf "\n" # replace observed controls with permuted controls

test -f ${OUTPUT}.tmp && echo rm ${OUTPUT}.tmp
test -f ${OUTPUT}.tmp && rm ${OUTPUT}.tmp && printf "\n"
