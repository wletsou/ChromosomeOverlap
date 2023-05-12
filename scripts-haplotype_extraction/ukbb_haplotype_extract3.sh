#! /bin/bash
set -e

# for extracting haplotypes from topmed-phased vcf files for UKBB and DRIVE at a common set of SNPs

# bsub -P SJLIFE -J ukbb_haplotype_extract3.chr11.69231642-69431642 -oo ukbb_haplotype_extract3.chr11.69231642-69431642.out -eo ukbb_haplotype_extract3.chr11.69231642-69431642.err -R "rusage[mem=10000]" -q large_mem "sh /home/wletsou/scripts/ukbb_haplotype_extract3.sh ukbb_bca_cases.indiv,ukbb_bca_controls.indiv,dbgap28544_cases.indiv,dbgap28544_controls.indiv ukbb_snp_list.chr11.69231642-69431642.txt chr11.qced.anno.info 11 69231642,69431642 ukbb.topmed.chr11.69231642-69431642.hg38.vcf.gz,dbgap28544_cases+dbgap28544_controls.chr11.69231642-69431642.vcf.gz"

POPULATION=$1 # commas-separated list of population indiv files
SNPS=$2 # optional comma-separated list of snps or file with one rsid snp per line
INFO=$3 # $INFO, used for converting hg38 to hg19
CHR=$4 # single chromosome number
BP_RANGE=$5 # comma-separated list from_bp,to_bp of region on chromosome CHR
VCF_FILES=$6 # comma-separated list of full paths to vcf files
DIRECTORY=$7 # PWD by default
HOME_DIR=$8 # location of program files

module load bcftools/1.10.2
module load tabix/0.2.6

if [ -z $HOME_DIR ];
then
  unset HOME_DIR
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
cd $DIRECTORY

BP_RANGE=($(echo $BP_RANGE | perl -pne 's/([0-9]+)[,]*/$1 /g'))

POPULATION=($(echo $POPULATION | perl -pne 's/[,]/ /g'))

test -f $INFO || (>&2 echo "SNP info file not found."; exit 1)

test -z $VCF_FILES && (>&2 echo "VCF files not supplied"; exit 1)
test -f $VCF_FILES || (>&2 echo "VCF files $VCF_FILES not found"; exit 1)

VCF_FILES=($(echo $VCF_FILES | sed 's/[,]/ /g'))

echo Phased haplotypes file$( (( ${#VCF_FILES[@]}==1 )) && echo " is" || echo "s are" ) ${VCF_FILES[*]}
printf "\n"

test -f vcf_list.txt && rm vcf_list.txt
touch vcf_list.txt
for ((i=0;i<${#VCF_FILES[@]};i++))
do
  echo echo ${VCF_FILES[i]} \>\> vcf_list.txt
  echo ${VCF_FILES[i]} >> vcf_list.txt
  printf "\n"
  if [ ! -f ${VCF_FILES[i]}.csi ] && [ ! -f ${VCF_FILES[i]}.tbi ] # create index if it does not exist
  then
    echo Generate .tbi index file
    echo tabix \-p vcf ${VCF_FILES[i]}
    tabix -p vcf ${VCF_FILES[i]}
    printf "\n"
  fi
done

if (( ${#VCF_FILES[@]}>1 )) && [ -f vcf_list.txt ]
then
  echo bcftools merge \-O z \-m snps \-l vcf_list.txt \-o $(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').combined.vcf.gz
  bcftools merge -O z -m snps -l vcf_list.txt -o $(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').combined.vcf.gz
  echo tabix \-p vcf $(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').combined.vcf.gz
  tabix -p vcf $(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').combined.vcf.gz
  test -f $(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').combined.vcf.gz && VCF_FILES=$(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').combined.vcf.gz
fi
test -f vcf_list.txt && rm vcf_list.txt

# specify region to extract: selected snps or entire range
if [ ! -z $SNPS ]
then
  if [ -f $SNPS ]
  then
    SNPS=$(cat $SNPS | awk '($1 ~ /rs[0-9]*/){printf "%s%s",(found>0?",":""),$1; found+=1} END{printf "\n"}') # snps supplied as a file, translated to comma-separated list, only keep rs snps
  fi
  region_str=$(awk 'BEGIN{n=split("'${SNPS}'",array,","); for (i=1;i<=n;i++) {snps[array[i]]} } ($8 in snps){printf "%schr%s:%s",(found>0?",":""),$2,$3; found+=1} END{printf "\n"}' $INFO) # get hg38 chromosome & coordinates of each snp
fi
if [ -z $region_str ]
then
  # find hg38 limits corresponding to hg19 limits, chromosome field does begin with "chr"
  range=($(awk 'NR>1{b=gensub("chr'$CHR':([0-9]*):.*","\\1","g",$7); if (b>='${BP_RANGE[0]}' && b<='${BP_RANGE[1]}') {print $3} }' $INFO | sort -gk1 | awk 'NR==1{print $0} {row=$0} END{print row}')) # prints first and last line of sorted list of within-limits hg38 positions (3rd field), based on h19 value in 7th field
  region_str=chr${CHR}:${range[0]}\-${range[1]}
fi
echo $region_str
printf "\n"

samples_file=$(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').samples
echo Get unique individuals in supplied .indiv file in order:
echo awk \'{seen[NR]=\$1\;} END{n = length\(seen\)\; for \(i=1\;i\<=n\;i++\) {print seen[i]} }\' ${POPULATION[*]} \> $samples_file # unique samples in .indiv file for population i
awk '{seen[NR]=$1;} END{n = length(seen); for (i=1;i<=n;i++) {print seen[i]} }' ${POPULATION[*]} > $samples_file # unique samples in .indiv file for population i, printed in same order as .indiv file
printf "\n"

transpose_file=haplotype_estimates_transpose.$(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
echo Get doubled header \(each subject id twice\):
echo bcftools view  \-r $region_str \--samples-file $samples_file \--force-samples \-e MAF==0 \--output-type v \--header-only $VCF_FILES \| awk \'{printf \"sjlife\\t\"\; for \(i=10\;i\<=NF\;i++\) {printf \"%s\\t%s%s\",\$i,\$i,\(i\<NF?\"\\t\":\"\\n\"\)} }\' \| tail \-1 \> $transpose_file
bcftools view -r $region_str --samples-file $samples_file --force-samples -e MAF==0 -m2 -M2 -v snps --output-type v --header-only $VCF_FILES | awk '{printf "sjlife\t"; for (i=10;i<=NF;i++) {printf "%s\t%s%s",$i,$i,(i<NF?"\t":"\n")} }' | tail -1 > $transpose_file
printf "\n"

new_vcf_file=$(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.haplotypes.hg38.vcf
echo Extract individuals from .vcf:
echo bcftools view \-r $region_str \--samples-file $samples_file \--force-samples \--output-type v \-e MAF==0 \--no-header \--output-file $new_vcf_file $VCF_FILES
bcftools view -r $region_str --samples-file $samples_file --force-samples -m2 -M2 -v snps --output-type v -e MAF==0 --no-header --output-file $new_vcf_file $VCF_FILES
printf "\n"

echo Compress and index vcf file:
echo bcftools view \-r $region_str \--samples-file $samples_file \--force-samples \--output-type z \-e MAF==0 \--output-file $(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.hg38.vcf.gz $VCF_FILES
bcftools view -r $region_str --samples-file $samples_file --force-samples -m2 -M2 -v snps --output-type z -e MAF==0 --output-file $(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.hg38.vcf.gz $VCF_FILES
printf "\n"

echo tabix \-p vcf $(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.hg38.vcf.gz
tabix -p vcf $(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').hg19_chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.hg38.vcf.gz
printf "\n"

echo Translate to hg19: # column 6 of $INFO is chr:hg38_pos:ref:alt, column 8 is rsid, column 16 is hg19 position; only take snps whose positions appear at most once (diallelic in population) in the (UKBB) INFO file
echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{b=gensub\(\"chr.*:\([0-9]\*\):.*\",\"\\\\1\",\"g\",\$6\)\; rsid[b]=\$8\; seen[b]+=1\; next} \(\$2 in rsid \&\& seen[\$2]==1\){\$3=rsid[\$2]\; print \$0}\' $INFO $new_vcf_file \> ${new_vcf_file}.tmp
awk 'BEGIN{OFS="\t"} NR==FNR{b=gensub("chr.*:([0-9]*):.*","\\1","g",$6); rsid[b]=$8; seen[b]+=1; next} ($2 in rsid && seen[$2]==1){$3=rsid[$2]; print $0}' $INFO $new_vcf_file > ${new_vcf_file}.tmp
printf "\n"

# remove SNPs whose positions appear more than once (at least triallelic in combined population); removes SNPs not diallelic in dbgap INFO file
echo awk \'NR==FNR{seen[\$2]+=1\; next} \(seen[\$2]==1\){print \$0}\' ${new_vcf_file}.tmp ${new_vcf_file}.tmp \> $new_vcf_file
awk 'NR==FNR{seen[$2]+=1; next} (seen[$2]==1){print $0}' ${new_vcf_file}.tmp ${new_vcf_file}.tmp > $new_vcf_file && printf "\n"

# test -f ${new_vcf_file}.tmp && echo mv ${new_vcf_file}.tmp $new_vcf_file && printf "\n"
# test -f ${new_vcf_file}.tmp && mv ${new_vcf_file}.tmp $new_vcf_file

echo Extract haplotypes as rsid x chromosome: # use most likely allele call
echo awk \'BEGIN{OFS=\"\\t\"} {printf \"%s_%s\\t\",\$3,\$5\; for \(i=10\;i\<=NF\;i++\) {p=gensub\(/\([0-9.]+\)[\|/]\([0-9.]+\).*/,\"\\\\1\",\"g\",\$i\)\; m=gensub\(/\([0-9.]+\)[\|/]\([0-9.]+\).*/,\"\\\\2\",\"g\",\$i\)\; printf \"%s\\t%s%s\",p,m,\(i\<NF?\"\\t\":\"\\n\"\)} }\' $new_vcf_file \>\> $transpose_file
awk 'BEGIN{OFS="\t"} {printf "%s_%s\t",$3,$5; for (i=10;i<=NF;i++) {p=gensub(/([0-9.]+)[|/]([0-9.]+).*/,"\\1","g",$i); m=gensub(/([0-9.]+)[|/]([0-9.]+).*/,"\\2","g",$i); printf "%s\t%s%s",p,m,(i<NF?"\t":"\n")} }' $new_vcf_file >> $transpose_file
printf "\n"

echo Check for alleles not present in some population:
echo awk \'NR==1{print \$0} NR\>1{for \(i=2\;i\<=NF\;i++\) {if \(\$i !\~ /[0-9]/\) {next} }\; print \$0 }\' $transpose_file \> ${transpose_file%.txt}.tmp
awk 'NR==1{print $0} NR>1{for (i=2;i<=NF;i++) {if ($i !~ /[0-9]/) {next} }; print $0 }' $transpose_file > ${transpose_file%.txt}.tmp
printf "\n"

test -f ${transpose_file%.txt}.tmp && echo mv ${transpose_file%.txt}.tmp $transpose_file && printf "\n"
test -f ${transpose_file%.txt}.tmp && mv ${transpose_file%.txt}.tmp $transpose_file

echo Remove unshared alleles from vcf file:
echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{b=gensub\(\"\(.*\)_.*\",\"\\\\1\",\"g\",\$1\)\; snp[b]\; next} \(\$3 in snp\){print \$0}\' $transpose_file $new_vcf_file \> ${new_vcf_file}.tmp
awk 'BEGIN{OFS="\t"} NR==FNR{b=gensub("(.*)_.*","\\1","g",$1); snp[b]; next} ($3 in snp){print $0}' $transpose_file $new_vcf_file > ${new_vcf_file}.tmp

test -f ${new_vcf_file}.tmp && echo mv ${new_vcf_file}.tmp $new_vcf_file && printf "\n"
test -f ${new_vcf_file}.tmp && mv ${new_vcf_file}.tmp $new_vcf_file

haplotypes_file=haplotype_estimates.$(echo ${POPULATION[*]%.indiv*} | sed 's/ /+/g').chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
echo Transpose haplotypes into chromosome x rsid:
echo awk \'BEGIN{OFS=\"\\t\"}\; {for\(j=1\;j\<=NF\;j++\) {a[NR,j]=\$j\; n_rows=NR\; n_cols=\(n_cols\<NF?NF:n_cols\)} } END{for \(j=1\;j\<=n_cols\;j++\) {for \(i=1\;i\<=n_rows\;i++\) {printf \"%s%s\",a[i,j],\(i==n_rows?\"\\n\":\"\\t\"\)} } }\' $transpose_file \> $haplotypes_file
awk 'BEGIN{OFS="\t"}; {for(j=1;j<=NF;j++) {a[NR,j]=$j; n_rows=NR; n_cols=(n_cols<NF?NF:n_cols)} } END{for (j=1;j<=n_cols;j++) {for (i=1;i<=n_rows;i++) {printf "%s%s",a[i,j],(i==n_rows?"\n":"\t")} } }' $transpose_file > $haplotypes_file
printf "\n"

for ((i=0;i<${#POPULATION[@]};i++))
do
  echo Extract ${POPULATION[i]%.indiv*} subjects as chromosome x rsid:
  echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{id[\$1]\; next} \(\$1 in id \|\| FNR==1\){print \$0}\' ${POPULATION[i]} $haplotypes_file \> haplotype_estimates.${POPULATION[i]%.indiv*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
  awk 'BEGIN{OFS="\t"} NR==FNR{id[$1]; next} ($1 in id || FNR==1){print $0}' ${POPULATION[i]} $haplotypes_file > haplotype_estimates.${POPULATION[i]%.indiv*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
  printf "\n"

  echo Transpose ${POPULATION[i]%.indiv*} haplotypes into rsid x chromosome:
  echo awk \'BEGIN{OFS=\"\\t\"}\; {for\(j=1\;j\<=NF\;j++\) {a[NR,j]=\$j\; n_rows=NR\; n_cols=\(n_cols\<NF?NF:n_cols\)} } END{for \(j=1\;j\<=n_cols\;j++\) {for \(i=1\;i\<=n_rows\;i++\) {printf \"%s%s\",a[i,j],\(i==n_rows?\"\\n\":\"\\t\"\)} } }\' haplotype_estimates.${POPULATION[i]%.indiv*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt \> haplotype_estimates_transpose.${POPULATION[i]%.indiv*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
  awk 'BEGIN{OFS="\t"}; {for(j=1;j<=NF;j++) {a[NR,j]=$j; n_rows=NR; n_cols=(n_cols<NF?NF:n_cols)} } END{for (j=1;j<=n_cols;j++) {for (i=1;i<=n_rows;i++) {printf "%s%s",a[i,j],(i==n_rows?"\n":"\t")} } }' haplotype_estimates.${POPULATION[i]%.indiv*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt > haplotype_estimates_transpose.${POPULATION[i]%.indiv*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[1]}.txt
  printf "\n"
done
