#! /bin/bash
set -e

# computes the LD R2 and D' values between two haplotypes or SNPs

# sh /home/wletsou/scripts/ukbb_haplotype_haplotype_LD2.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0" "rs62237573_T=1"

POPULATION=$1 # comma-separated list of haplotype_estimates files
HAPLOTYPE_1=$2 # colon-separated list of comma-separated lists rsid1_allele=[0,1]
HAPLOTYPE_2=$3 # colon-separated list of comma-separated lists rsid1_allele=[0,1]
DIRECTORY=$4
HOME_DIR=$5

if [ -z $HOME_DIR ];
then
  unset HOME_DIR
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY

POPULATION=($(echo $POPULATION | perl -pne 's/([^,]+)[,]*/$1 /g'))

if [ -z $HAPLOTYPE_2 ]
then
  HAPLOTYPE_2=$HAPLOTYPE_1
fi

HAPLOTYPE_1=($(echo $HAPLOTYPE_1 | perl -pne 's/[:]/ /g'))
HAPLOTYPE_2=($(echo $HAPLOTYPE_2 | perl -pne 's/[:]/ /g'))

for ((i=0;i<${#HAPLOTYPE_1[@]};i++))
do
  for ((j=0;j<${#HAPLOTYPE_2[@]};j++))
  do
    echo Haplotype $i vs. haplotype $j:
    # if haplotypes are in the form rsid1_[ATCG]=[01],rsid2_[ATCG]=[01]

    awk 'BEGIN{OFS="\t"; n=split("'${HAPLOTYPE_1[i]}'",array1,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array1[i]); b=gensub("(.*)=(.*)","\\2","g",array1[i]); haplotype1_snps[a]=b}; m=split("'${HAPLOTYPE_2[j]}'",array2,","); for (i=1;i<=m;i++) {c=gensub("(.*)=(.*)","\\1","g",array2[i]); d=gensub("(.*)=(.*)","\\2","g",array2[i]); haplotype2_snps[c]=d} } NR==1{ for (i = 2;i<=NF;i++) {if ($i in haplotype1_snps) {hap1_col[i]=$i}; if ($i in haplotype2_snps) {hap2_col[i]=$i} }; if (length(hap1_col)!=n || length(hap2_col)!=m) {print "'${HAPLOTYPE_1[i]}'","'${HAPLOTYPE_2[j]}'","NA","NA"; exit 1} } NR>1{chromosome_count=1; hap1_count=1; for (j in hap1_col) {if ($j!=haplotype1_snps[hap1_col[j]]) {hap1_count=0}; if ($j !~ /[0-9]/) {chromosome_count=0;} }; hap2_count=1; for (j in hap2_col) {if ($j!=haplotype2_snps[hap2_col[j]]) {hap2_count=0}; if ($j !~ /[0-9]/) {chromosome_count=0;} }; hap1_total+=hap1_count; hap2_total+=hap2_count; combined_total+=hap1_count*hap2_count; chromosome_total+=chromosome_count } END{p_AB=combined_total/chromosome_total; p_A=hap1_total/chromosome_total; p_B=hap2_total/chromosome_total; print p_A,p_B,p_AB; if (p_A*(1-p_A)*p_B*(1-p_B)!=0) {R2 = (p_AB-p_A*p_B)^2/(p_A*(1-p_A)*p_B*(1-p_B))} else {R2="NA"}; D = p_AB-p_A*p_B; if (D < 0) {Dmax = -p_A*p_B; if ( -(1-p_A)*(1-p_B) > Dmax ) {Dmax = -(1-p_A)*(1-p_B)} } else if (D > 0) {Dmax = p_A*(1-p_B); if (p_B*(1-p_A) < Dmax) {Dmax = p_B*(1-p_A)} } else {Dmax=0}; if (Dmax!=0) {Dprime=D/Dmax} else {Dprime="NA"}; print "'${HAPLOTYPE_1[i]}'","'${HAPLOTYPE_2[j]}'",R2,Dprime}' ${POPULATION[0]}
    printf "\n"
  done
done
