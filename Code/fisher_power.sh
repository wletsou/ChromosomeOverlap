#! /bin/bash
set -e

# /home/wletsou/scripts/fisher_power.sh /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.34/phase2/haplotype_estimates.ukbb_bca_cases.chr11.68850000-69231641.txt,/scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.34/phase2/haplotype_estimates.ukbb_bca_controls.chr11.68850000-69231641.txt /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.34/phase2/haplotype_estimates.dbgap28544_cases.chr11.68850000-69231641.txt,/scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.34/phase2/haplotype_estimates.dbgap28544_controls.chr11.68850000-69231641.txt rs74674490_G=0,rs3892895_G=1,rs74897859_T=0,rs67039008_A=1,rs55703863_C=0,rs11228565_A=0,rs116951063_T=0,rs11825015_T=0,rs74551015_A=0,rs11539762_A=0,rs7124547_T=1,rs11603814_G=1,rs79487139_T=0,rs4275647_T=0,rs11263641_C=0,rs76855139_T=0,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0 # h4 reduced

REF_POPULATION=$1 # optional comma-separated list of two haplotype_estimates files for discovery
TEST_POPULATION=$2 # comma-separated list of two haplotype_estimates files for replications
HAPLOTYPES=$3 # file or colon-separated list of comma-separated lists of rsid_allele=[0/1]
PHENOTYPES=$4 # optional phenotypes file of subjects to include
OR=$5 # OR to detect (optional if REF_POPULATION supplied, but overrides REF_POPULATION)
ALPHA=$6 # false-positive error rate
DOMINANT=$7 # string constant, TRUE if nonempty; whether to use dominant model
DIRECTORY=$6
HOME_DIR=$8

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

REF_POPULATION=($(echo $REF_POPULATION | perl -pne 's/([^,]+)[,]*/$1 /g'))
TEST_POPULATION=($(echo $TEST_POPULATION | perl -pne 's/([^,]+)[,]*/$1 /g'))

if [ -f $HAPLOTYPES ]
then
  HAPLOTYPES=$(cut -f1 $HAPLOTYPES | tr "\n" ":" | sed 's/:$/\n/g') # translate first column of file into colon-separated list
fi
HAPLOTYPES=($(echo $HAPLOTYPES | perl -pne 's/[:]/ /g'))

if [ -z $ALPHA ]
then
  ALPHA=0.05
fi

if [ ! -z $DOMINANT ]
then
  DOMINANT="-i"
fi

for ((i=0;i<${#HAPLOTYPES[@]};i++))
do
  if [ -z $OR ] && ((${#REF_POPULATION[@]}==2))
  then
    unset pop1 pop2

    echo Calculate OR:
    echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${HAPLOTYPES[i]}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; haplotype_snps[a]=b} } NR==1{ for \(i = 2\;i\<=NF\;i++\) {if \(\$i in haplotype_snps\) {hap_col[i]=\$i} } if \(length\(hap_col\)!=n\) {print \"NA\",\"NA\"\; exit 1} } NR\>1{seen[\$1]+=1\; hap_count=1\; for \(j in hap_col\) {if \(\$j!=haplotype_snps[hap_col[j]]\) {hap_count=0}\; if \(\$j !~ /[0-9]/\) {hap_count=\"NA\"} }\; if \(hap_count\>=0\) {sum+=hap_count\; total+=1} } END{print hap_count,total}\' ${REF_POPULATION[0]}
    pop1=($(awk 'BEGIN{OFS="\t"; n=split("'${HAPLOTYPES[i]}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } NR==1{ for (i = 2;i<=NF;i++) {if ($i in haplotype_snps) {hap_col[i]=$i} } if (length(hap_col)!=n) {print "NA","NA"; exit 1} } NR>1{seen[$1]+=1; hap_count=1; for (j in hap_col) {if ($j!=haplotype_snps[hap_col[j]]) {hap_count=0}; if ($j !~ /[0-9]/) {hap_count="NA"} }; if (hap_count>=0) {sum+=hap_count; total+=1} } END{print sum,total}' <( (test ! -z $PHENOTYPES && test -f $PHENOTYPES) && awk 'NR==FNR{seen[$1]; next} (FNR==1){print $0} (FNR>1 && $1 in seen){print $0}' $PHENOTYPES ${REF_POPULATION[0]} || cat ${REF_POPULATION[0]}))) # haplotype counts among cases included in phenotypes model (or all cases)
    declare -p pop1 && printf "\n" # cases

    echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${HAPLOTYPES[i]}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; haplotype_snps[a]=b} } NR==1{ for \(i = 2\;i\<=NF\;i++\) {if \(\$i in haplotype_snps\) {hap_col[i]=\$i} } if \(length\(hap_col\)!=n\) {print \"NA\",\"NA\"\; exit 1} } NR\>1{seen[\$1]+=1\; hap_count=1\; for \(j in hap_col\) {if \(\$j!=haplotype_snps[hap_col[j]]\) {hap_count=0}\; if \(\$j !~ /[0-9]/\) {hap_count=\"NA\"} }\; if \(hap_count\>=0\) {sum+=hap_count\; total+=1} } END{print hap_count,total}\' ${REF_POPULATION[1]}
    pop2=($(awk 'BEGIN{OFS="\t"; n=split("'${HAPLOTYPES[i]}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } NR==1{ for (i = 2;i<=NF;i++) {if ($i in haplotype_snps) {hap_col[i]=$i} } if (length(hap_col)!=n) {print "NA","NA"; exit 1} } NR>1{seen[$1]+=1; hap_count=1; for (j in hap_col) {if ($j!=haplotype_snps[hap_col[j]]) {hap_count=0}; if ($j !~ /[0-9]/) {hap_count="NA"} }; if (hap_count>=0) {sum+=hap_count; total+=1} } END{print sum,total}' <( (test ! -z $PHENOTYPES && test -f $PHENOTYPES) && awk 'NR==FNR{seen[$1]; next} (FNR==1){print $0} (FNR>1 && $1 in seen){print $0}' $PHENOTYPES ${REF_POPULATION[1]} || cat ${REF_POPULATION[1]}))) # haplotype counts among controls included in phenotypes model (or all controls)
    declare -p pop2 && printf "\n" # controls

    OR=$(awk 'BEGIN{printf "%0.2e\n",('${pop1[0]}'/'${pop1[1]}')/(1-'${pop1[0]}'/'${pop1[1]}')/(('${pop2[0]}'/'${pop2[1]}')/(1-'${pop2[0]}'/'${pop2[1]}'))}')
    # OR=$(awk 'BEGIN{print ('${pop1[0]}'/'${pop1[1]}')/(1-'${pop1[0]}'/'${pop1[1]}')/(('${pop2[0]}'/'${pop2[1]}')/(1-'${pop2[0]}'/'${pop2[1]}'))}')
    echo OR = $OR && printf "\n"
  elif [ -z $OR ] && ((${#REF_POPULATION[@]}!=2))
  then
    (>&2 echo "OR or reference population not supplied."; exit 1)
  fi

  unset pop1 pop2

  echo Get counts of haplotype $((i+1)):
  echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${HAPLOTYPES[i]}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; haplotype_snps[a]=b} } NR==1{ for \(i = 2\;i\<=NF\;i++\) {if \(\$i in haplotype_snps\) {hap_col[i]=\$i} } if \(length\(hap_col\)!=n\) {print \"NA\",\"NA\"\; exit 1} } NR\>1{seen[\$1]+=1\; hap_count=1\; for \(j in hap_col\) {if \(\$j!=haplotype_snps[hap_col[j]]\) {hap_count=0}\; if \(\$j !~ /[0-9]/\) {hap_count=\"NA\"} }\; if \(hap_count\>=0\) {sum+=hap_count\; total+=1} } END{print hap_count,total}\' ${TEST_POPULATION[0]}
  pop1=($(awk 'BEGIN{OFS="\t"; n=split("'${HAPLOTYPES[i]}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } NR==1{ for (i = 2;i<=NF;i++) {if ($i in haplotype_snps) {hap_col[i]=$i} } if (length(hap_col)!=n) {print "NA","NA"; exit 1} } NR>1{seen[$1]+=1; hap_count=1; for (j in hap_col) {if ($j!=haplotype_snps[hap_col[j]]) {hap_count=0}; if ($j !~ /[0-9]/) {hap_count="NA"} }; if (hap_count>=0) {sum+=hap_count; total+=1} } END{print sum,total}' <( (test ! -z $PHENOTYPES && test -f $PHENOTYPES) && awk 'NR==FNR{seen[$1]; next} (FNR==1){print $0} (FNR>1 && $1 in seen){print $0}' $PHENOTYPES ${TEST_POPULATION[0]} || cat ${TEST_POPULATION[0]}))) # haplotype counts among cases included in phenotypes model (or all cases)
  declare -p pop1 && printf "\n" # cases

  echo awk \'BEGIN{OFS=\"\\t\"\; n=split\(\"${HAPLOTYPES[i]}\",array,\",\"\)\; for \(i=1\;i\<=n\;i++\) {a=gensub\(\"\(.*\)=\(.*\)\",\"\\\\1\",\"g\",array[i]\)\; b=gensub\(\"\(.*\)=\(.*\)\",\"\\\\2\",\"g\",array[i]\)\; haplotype_snps[a]=b} } NR==1{ for \(i = 2\;i\<=NF\;i++\) {if \(\$i in haplotype_snps\) {hap_col[i]=\$i} } if \(length\(hap_col\)!=n\) {print \"NA\",\"NA\"\; exit 1} } NR\>1{seen[\$1]+=1\; hap_count=1\; for \(j in hap_col\) {if \(\$j!=haplotype_snps[hap_col[j]]\) {hap_count=0}\; if \(\$j !~ /[0-9]/\) {hap_count=\"NA\"} }\; if \(hap_count\>=0\) {sum+=hap_count\; total+=1} } END{print hap_count,total}\' ${TEST_POPULATION[1]}
  pop2=($(awk 'BEGIN{OFS="\t"; n=split("'${HAPLOTYPES[i]}'",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } NR==1{ for (i = 2;i<=NF;i++) {if ($i in haplotype_snps) {hap_col[i]=$i} } if (length(hap_col)!=n) {print "NA","NA"; exit 1} } NR>1{seen[$1]+=1; hap_count=1; for (j in hap_col) {if ($j!=haplotype_snps[hap_col[j]]) {hap_count=0}; if ($j !~ /[0-9]/) {hap_count="NA"} }; if (hap_count>=0) {sum+=hap_count; total+=1} } END{print sum,total}' <( (test ! -z $PHENOTYPES && test -f $PHENOTYPES) && awk 'NR==FNR{seen[$1]; next} (FNR==1){print $0} (FNR>1 && $1 in seen){print $0}' $PHENOTYPES ${TEST_POPULATION[1]} || cat ${TEST_POPULATION[1]}))) # haplotype counts among controls included in phenotypes model (or all controls)
  declare -p pop2 && printf "\n" # controls

  echo Power to detect OR = $OR:
  echo Rscript ${HOME_DIR}/fisher_power2.R \-x ${pop1[0]} \-n ${pop1[1]} \-y ${pop2[0]} \-m ${pop2[1]} \-o $OR \-a $ALPHA $DOMINANT
  Rscript ${HOME_DIR}/fisher_power2.R -x ${pop1[0]} -n ${pop1[1]} -y ${pop2[0]} -m ${pop2[1]} -o $OR -a $ALPHA $DOMINANT
done
