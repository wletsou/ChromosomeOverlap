#! /bin/bash
set -e

# /home/wletsou/scripts/dbgap28544_ukbb_phenotype_combine2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/drive4william/flashpcar/ukbb.drive.pcs.txt ukb.bca_dbgap28544.pheno.txt

REF_PHENOTYPES=$1 # phenotypes file for getting age and affected status
JOIN_PHENOTYPES=$2 # phenotypes file to be matched with REF_PHENOTYPES
PCS=$3 # file with principal components for all subjects in REF and JOIN
OUTPUT=$4 # optional output file name
DIRECTORY=$5
HOME_DIR=$6

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

test -f ${REF_PHENOTYPES} || (>&2 echo "Reference-phenotypes file not found."; exit 1)

test -f ${JOIN_PHENOTYPES} || (>&2 echo "Join-phenotypes file not found."; exit 1)

if [ -z $OUTPUT ]
then
  OUTPUT=pheno_combined.txt
fi

echo Make dbgap28544-like phenotype file: && printf "\n"

awk 'NR==FNR && NR==1{for (i=1;i<=NF;i++) {if ($i=="sid") {col_sid=i} }; if (!(col_sid+0>0)) {print "No sid column found"; exit 1} } NR==FNR{x=$1; $1=$col_sid; $col_sid=x; seen[$1]=$0; next} ($1 in seen){printf "%s\t%s\t%s\n",seen[$1],$2,$3}' $PCS <(awk '{print $0}' <(awk 'BEGIN{OFS="\t"} NR==1{for (i=1;i<=NF;i++) {if ($i=="ageonset") {col_age=i; age[$1]=$i}; if ($i=="BCa") {col_diag=i; diag[$1]=$i} }; if (!(col_age+0>0) || !(col_diag+0>0)) {print "Age or diagnosis column not found"; exit 1}; printf "%s\t%s\t%s\n",$1,"ageonset","BCa"} NR>1{age[$1]=$col_age; diag[$1]=$col_diag; printf "%s\t%s\t%s\n",$1,age[$1],diag[$1]}' ${REF_PHENOTYPES}) <(awk 'BEGIN{OFS="\t"} NR==1{for (i=1;i<=NF;i++) {if ($i=="age.end") {col_age=i; age[$1]=$i}; if ($i=="diag") {col_diag=i; diag[$1]=$i} }; if (!(col_age+0>0) || !(col_diag+0>0)) {print "Age or diagnosis column not found"; exit 1} } NR>1{age[$1]=$col_age; diag[$1]=$col_diag; printf "%s\t%s\t%s\n",$1,age[$1],diag[$1]}' ${JOIN_PHENOTYPES}) ) > $OUTPUT

echo Combined phenotypes file is ${OUTPUT}
