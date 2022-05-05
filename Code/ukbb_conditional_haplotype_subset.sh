#! /bin/bash
set -e

# sh /home/wletsou/scripts/ukbb_conditional_haplotype_subset.sh haplotype_estimates.ukbb_bca_20200116_cases.chr11.68700000-69700000.txt,haplotype_estimates.ukbb_bca_20200116_controls.chr11.68700000-69700000.txt rs79373485_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs79442425_T=0,rs111929748_T=0,rs637185_T=0,rs4980661_A=0,rs2046494_C=0 rs61881030_A,rs149949063_G,rs1123665_G,rs3829241_A,rs74674490_G,rs72930631_C,rs11228502_C,rs12577213_A,rs73520512_T,rs3892895_G,rs10896431_T,rs7116054_T,rs74897859_T,rs17308715_T,rs67039008_A,rs55703863_C,rs17149775_G,rs117236867_G,rs11228530_C,rs17309046_C,rs56134379_A,rs72932540_G,rs12418948_T,rs78033785_T,rs74342245_A,rs35435383_A,rs4495899_G,rs117694794_T,rs10792029_G,rs10896445_C,rs78208050_G,rs61881109_A,rs75590802_A,rs11228565_A,rs7937094_T,rs10896449_G,rs7130881_G,rs7119988_A,rs111762835_A,rs4930672_A,rs11228599_T,rs116951063_T,rs118004919_A,rs35637432_A,rs34737133_T,rs11825015_T,rs78656034_T,rs74551015_A,rs117821797_T,rs12275849_C,rs56103266_C,rs117131950_A,rs7103126_C,rs116926312_T,rs11539762_A,rs142581206_T,rs12274095_A,rs73512137_T,rs11601693_A,rs10896461_G,rs78919175_T,rs7124547_T,rs11603814_G,rs72932105_T,rs12789955_G,rs7481709_C,rs115223540_T,rs28378931_G,rs7102705_G,rs117186144_T,rs79487139_T,rs11600497_A,rs72932198_A,rs12802601_C,rs12796465_G,rs4275647_T,rs61882193_G,rs12790064_G,rs11263638_A,rs11263641_C,rs76855139_T,rs74471298_A,rs10160464_G,rs4980785_T,rs7105934_A,rs6606643_A,rs7927237_C

# sh /home/wletsou/scripts/ukbb_conditional_haplotype_subset.sh haplotype_estimates.ukbb_bca_cases.chr11.68700000-69700000.txt,haplotype_estimates.ukbb_bca_controls.chr11.68700000-69700000.txt s79373485_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs79442425_T=0,rs111929748_T=0,rs637185_T=0,rs4980661_A=0,rs2046494_C=0 "rs2376558_C,rs896973_T,rs61881030_A,rs149949063_G,rs144905179_A,rs1123665_G,rs3829241_A,rs74674490_G,rs72930627_A,rs72930631_C,rs11228502_C,rs12577213_A,rs73520512_T,rs3892895_G,rs10896431_T,rs7116054_T,rs74897859_T,rs17308715_T,rs67039008_A,rs55703863_C,rs17149775_G,rs117236867_G,rs11228530_C,rs17309046_C,rs56134379_A,rs72932540_G,rs12418948_T,rs78033785_T,rs74342245_A,rs35435383_A,rs4495899_G,rs117694794_T,rs10792029_G,rs10896445_C,rs78208050_G,rs61881109_A,rs75590802_A,rs11228565_A,rs7937094_T,rs7931342_G,rs10896449_G,rs7130881_G,rs7119988_A,rs111762835_A,rs4930672_A,rs11228599_T,rs7940107_A,rs72930267_T,rs7946255_C,rs116951063_T,rs118004919_A,rs35637432_A,rs34737133_T,rs11825015_T,rs78656034_T,rs74551015_A,rs117821797_T,rs12275849_C,rs56103266_C,rs117131950_A,rs7103126_C,rs116926312_T,rs11539762_A,rs142581206_T,rs11228610_C,rs12274095_A,rs73512137_T,rs11601693_A,rs10896461_G,rs12224376_A,rs78919175_T,rs7124547_T,rs11603814_G,rs72932105_T,rs12789955_G,rs7481709_C,rs115223540_T,rs28378931_G,rs7102705_G,rs117186144_T,rs79487139_T,rs11600497_A,rs72932198_A,rs12802601_C,rs12796465_G,rs4275647_T,rs61882193_G,rs12790064_G,rs11263638_A,rs11263641_C,rs76855139_T,rs74471298_A,rs10160464_G,rs4980785_T,rs7105934_A,rs6606643_A,rs7927237_C,rs79373485_T" h1

POPULATION=$1 # haplotype_estimates file for population
HAPLOTYPE=$2 # comma-separated list of rsid_allele,... to condition on
SNPS=$3 # SNPs to be included in haplotype
NAME=$4
DIRECTORY=$5 # current working directory, specify replicate folder if REPLICATE_NO > 0
HOME_DIR=$6 # location of scripts to be run

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY

POPULATION=($(echo $POPULATION | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))

if [ ! -z $NAME ]
then
  NAME="${NAME}_"
fi

for ((i=0;i<${#POPULATION[@]};i++))
do
  char=$(($(awk '(NR==1){print NF}' ${POPULATION[j]} | wc -c)-1)) # number of characters in the value of the number of columns in the haplotypes file
  awk 'BEGIN{OFS="\t"; n=split("'${HAPLOTYPE}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b}; } NR==1{ delete array; m=split($0,snps,"\t"); for (i = 2;i<=NF;i++) {if ($i in haplotype_snps) {hap_col[i]=$i} }; print $0 } NR>1{hap_count=1; found=0; for (j in hap_col) {found+=1; if ($j!=haplotype_snps[hap_col[j]]) {hap_count=0; break} }; if (hap_count==1 && found==n) {print $0} }' ${POPULATION[i]} > ${POPULATION[i]%.*}.${NAME}subset.tmp

  if [ ! -z $SNPS ] # only print columns in SNPS variable
  then
    awk 'BEGIN{OFS="\t"} NR==1{printf "%s",$1; for (j = 1; j <= NF; j++) {if ("'$SNPS'" ~ $j) {col[sprintf("%0'${char}'d",j)]=j} }; m = asorti(col,ind); for (j = 1; j <= m; j++) {printf "\t%s%s",$col[sprintf("%0'${char}'d",ind[j])],(j<m?"":"\n") } } NR>1{printf "%s",$1; for (j = 1; j <= m; j++) {printf "\t%s%s",$col[sprintf("%0'${char}'d",ind[j])],(j<m?"":"\n")} }' ${POPULATION[i]%.*}.${NAME}subset.tmp > ${POPULATION[i]%.*}.${NAME}subset.txt

    test -f ${POPULATION[i]%.*}.${NAME}subset.tmp && rm ${POPULATION[i]%.*}.${NAME}subset.tmp
  else
    test -f ${POPULATION[i]%.*}.${NAME}subset.tmp && mv ${POPULATION[i]%.*}.${NAME}subset.tmp ${POPULATION[i]%.*}.${NAME}subset.txt
  fi
done
