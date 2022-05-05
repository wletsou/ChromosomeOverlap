#! /bin/bash
set -e

# /home/wletsou/scripts/controls_compare.sh /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.34/phase2/haplotype_estimates.ukbb_bca_controls.chr11.68850000-69231641.txt,/scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.34/phase2/haplotype_estimates.dbgap28544_controls.chr11.68850000-69231641.txt rs74674490_G=0,rs3892895_G=1,rs74897859_T=0,rs67039008_A=1,rs55703863_C=0,rs11228565_A=0,rs116951063_T=0,rs11825015_T=0,rs74551015_A=0,rs11539762_A=0,rs7124547_T=1,rs11603814_G=1,rs79487139_T=0,rs4275647_T=0,rs11263641_C=0,rs76855139_T=0,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0 # h4 reduced

# /home/wletsou/scripts/controls_compare.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/drive4william/UKBB_chr10.1/phase2.1/haplotype_estimates.ukbb_bca_controls.chr10.122825000-123200000.txt,/research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/drive4william/UKBB_chr10.1/phase2.1/haplotype_estimates.dbgap28544_controls.chr10.122825000-123200000.txt rs4752498_T=1,rs2244506_G=0,rs4436487_C=0,rs11199800_A=1,rs10886873_T=1,rs7076500_A=1,rs4752536_G=0,rs12256422_A=0,rs7071430_T=0,rs12220114_G=1,rs10788173_T=0,rs4752550_T=0,rs11199957_T=0,rs11199964_C=0,rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0 # h5 reduced

POPULATION=$1 # comma-separated list of two haplotype_estimates files
HAPLOTYPES=$2 # file or colon-separated list of comma-separated lists of rsid_allele=[0/1]
DIRECTORY=$3
HOME_DIR=$4

module load python/3.7.0
module load python/conda/3.8.1
module load conda3/5.1.0

if [ -z $HOME_DIR ];
then
  HOME_DIR=$(echo "/home/wletsou/scripts")
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY

POPULATION=($(echo $POPULATION | perl -pne 's/([^,]+)[,]*/$1 /g'))

if [ -f $HAPLOTYPES ]
then
  HAPLOTYPES=$(cut -f1 $HAPLOTYPES | tr "\n" ":" | sed 's/:$/\n/g') # translate first column of file into colon-separated list
fi
HAPLOTYPES=($(echo $HAPLOTYPES | perl -pne 's/[:]/ /g'))

for ((i=0;i<${#HAPLOTYPES[@]};i++))
do
  echo Get counts of haplotype $((i+1)):
  echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${HAPLOTYPES[i]}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; haplotype_snps[a]=b} } NR==1{ for \(i = 2\;i\<=NF\;i++\) {if \(\$i in haplotype_snps\) {hap_col[i]=\$i} } if \(length\(hap_col\)!=n\) {print \"NA\",\"NA\"\; exit 1} } NR\>1{seen[\$1]+=1\; hap_count=1\; for \(j in hap_col\) {if \(\$j!=haplotype_snps[hap_col[j]]\) {hap_count=0}\; if \(\$j !~ /[0-9]/\) {hap_count=\"NA\"} }\; if \(hap_count\>=0\) {sum+=hap_count\; total+=1} } END{print hap_count,total}\' ${POPULATION[0]}
  pop1=($(awk 'BEGIN{OFS="\t"; n=split("'${HAPLOTYPES[i]}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } NR==1{ for (i = 2;i<=NF;i++) {if ($i in haplotype_snps) {hap_col[i]=$i} } if (length(hap_col)!=n) {print "NA","NA"; exit 1} } NR>1{seen[$1]+=1; hap_count=1; for (j in hap_col) {if ($j!=haplotype_snps[hap_col[j]]) {hap_count=0}; if ($j !~ /[0-9]/) {hap_count="NA"} }; if (hap_count>=0) {sum+=hap_count; total+=1} } END{print sum,total}' ${POPULATION[0]}))
  declare -p pop1 && printf "\n"

  echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${HAPLOTYPES[i]}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; haplotype_snps[a]=b} } NR==1{ for \(i = 2\;i\<=NF\;i++\) {if \(\$i in haplotype_snps\) {hap_col[i]=\$i} } if \(length\(hap_col\)!=n\) {print \"NA\",\"NA\"\; exit 1} } NR\>1{seen[\$1]+=1\; hap_count=1\; for \(j in hap_col\) {if \(\$j!=haplotype_snps[hap_col[j]]\) {hap_count=0}\; if \(\$j !~ /[0-9]/\) {hap_count=\"NA\"} }\; if \(hap_count\>=0\) {sum+=hap_count\; total+=1} } END{print hap_count,total}\' ${POPULATION[1]}
  pop2=($(awk 'BEGIN{OFS="\t"; n=split("'${HAPLOTYPES[i]}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } NR==1{ for (i = 2;i<=NF;i++) {if ($i in haplotype_snps) {hap_col[i]=$i} } if (length(hap_col)!=n) {print "NA","NA"; exit 1} } NR>1{seen[$1]+=1; hap_count=1; for (j in hap_col) {if ($j!=haplotype_snps[hap_col[j]]) {hap_count=0}; if ($j !~ /[0-9]/) {hap_count="NA"} }; if (hap_count>=0) {sum+=hap_count; total+=1} } END{print sum,total}' ${POPULATION[1]}))
  declare -p pop2 && printf "\n"

  echo Test difference in frequencies $(awk 'BEGIN{OFS="\t"; print '${pop1[0]}'/'${pop1[1]}','${pop2[0]}'/'${pop2[1]}'}'):
  echo python ${HOME_DIR}/controls_frequency.py \-x1 ${pop1[0]} \-n1 ${pop1[1]} \-x2 ${pop2[0]} \-n2 ${pop2[1]} \-pooled
  python ${HOME_DIR}/controls_frequency.py -x1 ${pop1[0]} -n1 ${pop1[1]} -x2 ${pop2[0]} -n2 ${pop2[1]} -pooled && printf "\n"

  echo python ${HOME_DIR}/controls_frequency2.py \-x1 ${pop1[0]} \-n1 ${pop1[1]} \-x2 ${pop2[0]} \-n2 ${pop2[1]}
  python ${HOME_DIR}/controls_frequency2.py -x1 ${pop1[0]} -n1 ${pop1[1]} -x2 ${pop2[0]} -n2 ${pop2[1]} && printf "\n"
done
