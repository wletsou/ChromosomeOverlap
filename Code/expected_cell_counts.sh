#! /bin/bash
set -e

# get expected counts in a 2 x 2 table of cases/controls x haplotype/no haplotype (read a, b, c, d across rows)

# sh /home/wletsou/scripts/expected_cell_counts.sh haplotype_estimates.dbgap28544_cases.chr11.68850000-69231641.txt,haplotype_estimates.dbgap28544_controls.chr11.68850000-69231641.txt "rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs74471298_A=0,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" h2

POPULATION=$1 # two haplotype_estimates files, one for each population (case/controls) with subject id in first column and rsids in columns 2 and beyond
HAPLOTYPE=$2 # rsid_allele=[0/1],...
DIRECTORY=$3 # current working directory
HOME_DIR=$4 # location of scripts to be run

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

POPULATION=($(echo $POPULATION | perl -pne 's/([A-Za-z0-9_.+-]+)[,]*/$1 /g'))

awk 'BEGIN{OFS="\t"; n=split("'${HAPLOTYPE}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } (NR==FNR && NR==1){ for (i = 2;i<=NF;i++) {if ($i in haplotype_snps) {hap_col[i]=$i} }; next } FNR>1{seen[$1]+=1; hap_count=1; for (j in hap_col) {if ($j!=haplotype_snps[hap_col[j]]) {hap_count=0}; if ($j !~ /[0-9]/) {hap_count="NA"; next} }; hap_total+=hap_count; if (NR==FNR) {cases+=1} else {controls+=1} } END{n=cases+controls; a=cases*hap_total/n; b=controls*hap_total/n; c=cases*(n-hap_total)/n; d=controls*(n-hap_total)/n; printf "Expected cell counts: a = %0.6f, b = %0.6f, c = %0.6f, d = %0.6f\n",a,b,c,d}' ${POPULATION[*]}
printf "\n"
