#! /bin/bash
set -e

# /home/wletsou/scripts/ukbb_LD_block_haplotype.sh /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.34/phase2/haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69231641.txt rs74674490_G=0,rs72930631_C=0,rs11228502_C=0,rs10896431_T=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs78033785_T=0,rs117694794_T=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs116951063_T=0,rs11825015_T=0,rs78656034_T=0,rs7103126_C=1,rs12274095_A=0,rs12789955_G=0,rs72932198_A=0,rs61882193_G=0,rs12790064_G=1,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0 ~/Chromosome_Overlap_results/UKBB_chr11.34/phase2/1000GP/1000gp.hg19_11.68850000-69500000.hg38.res.txt /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.34/phase2/ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.haplotypes.hg38.vcf

POPULATION=$1 # comma-separated list of haplotype_estimates files
HAPLOTYPE=$2 # colon-separated list of comma-separated lists rsid1_allele=[0,1]
BIG_LD=$3 # BIG_LD result in region, delimiting LD blocks by chromosome position
VCF=$4 # vcf file for POPULATION (not g-zipped)
OUTPUT=$5 # optional name of output file
DIRECTORY=$6
HOME_DIR=$7

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY

if [ -z $OUTPUT ]
then
  OUTPUT=haplotype.LD_segments.txt
fi

# split into sub-haplotypes by LD block
hap_array=($(awk 'BEGIN{n=split("'${HAPLOTYPE}'",array,",")} (NR==FNR && NR>1){start_bp[NR-1]=$6; end_bp[NR-1]=$7; next} NR!=FNR{ for (j=1;j<=n;j++) {if (array[j]~$3) {total++; for (i=1;i<=length(start_bp);i++) {if ($2>=start_bp[i] && $2<=end_bp[i]) {found[i]+=1; if (found[i]==1 && total>1) {printf " "}; printf "%s%s",(found[i]>1?",":""),array[j]; } } } } } END{if (total>0) {printf "\n"} }' $BIG_LD $VCF))

# chromosome ranges of LD blocks
pos_array=($(awk 'BEGIN{n=split("'${HAPLOTYPE}'",array,",")} (NR==FNR && NR>1){chr[NR-1]=$1; start_bp[NR-1]=$6; end_bp[NR-1]=$7; next} NR!=FNR{ for (j=1;j<=n;j++) {if (array[j]~$3) {total++; for (i=1;i<=length(start_bp);i++) {if ($2>=start_bp[i] && $2<=end_bp[i]) {found[i]+=1; if (found[i]==1 && total>1) {printf " "}; if (found[i]==1) {printf "chr%s:%s-%s",chr[i],start_bp[i],end_bp[i]}; } } } } } END{if (total>0) {printf "\n"} }' $BIG_LD $VCF))

awk 'BEGIN{OFS="\t"; print "range","segment","frequency"}' > ${DIRECTORY}/${OUTPUT%.txt}.tmp

for ((i=0;i<${#hap_array[@]};i++))
do
  awk 'BEGIN{OFS="\t"; n=split("'${hap_array[i]}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } NR==1{ for (i=2;i<=NF;i++) {if ($i in haplotype_snps) {hap_col[i]=$i} } } NR>1{hap_count=1; for (j in hap_col) {if ($j!=haplotype_snps[hap_col[j]]) {hap_count=0}; if ($j !~ /[0-9]/) {hap_count="NA"; next} }; hap_total+=hap_count; chrom_total+=1 } END{print "'${pos_array[i]}'","'${hap_array[i]}'",hap_total/chrom_total}' ${POPULATION} >> ${OUTPUT%.txt}.tmp
done

# replace =[0/1] with allele
awk 'BEGIN{OFS="\t"} NR==FNR{ref[$3]=$4; alt[$3]=$5; next} (NR!=FNR && FNR==1){print $0} (NR!=FNR && FNR>1){str=""; found=0; n=split($2,array,","); for (i=1;i<=n;i++) {a=gensub("(.*)_[A-Z]=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); if (b==0) {str=sprintf("%s%s%s[%s]",str,(found>0?",":""),a,ref[a]); found++}; if (b==1) {str=sprintf("%s%s%s[%s]",str,(found>0?",":""),a,alt[a]); found++} }; $2=str; print $0 }' $VCF ${OUTPUT%.txt}.tmp > ${DIRECTORY}/${OUTPUT}

test -f ${OUTPUT%.txt}.tmp && rm ${OUTPUT%.txt}.tmp
