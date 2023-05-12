#! /bin/bash
set -e

# for extracting genotyped SNPs in UKBB data

# bsub -P SJLIFE -J ukbb_hybrid_haplotype3.chr11.68850000-69500000 -oo ukbb_hybrid_haplotype3.chr11.68850000-69231641.out -eo ukbb_hybrid_haplotype3.chr11.68850000-69500000.err -R "rusage[mem=20000]" "sh ukbb_hybrid_haplotype3.sh ukbb_bca_cases.indiv,ukbb_bca_controls.indiv 11 68850000,69500000 ukb.bca.hap.chr11.new.vcf.gz"

POPULATION=$1 # commas-separated .indiv files, e.g. ukbb_bca_20200116_cases.indiv or ukbb_bca_all_controls.indiv
CHR=$2 # single chromosome number
BP_RANGE=$3 # comma-separated list from_bp,to_bp of regions on chromosome CHR
VCF_FILE=$4 # vcf file
INFO=$5 # optional snp annotation file with CHR in second column, hg19 POS in third column, and SNP rsid in sixth column
DIRECTORY=$6 # PWD by default
HOME_DIR=$7 # location of program files

module load bcftools/1.10.2
module load tabix/0.2.6

printf "\n"
if [ -z $HOME_DIR ];
then
  unset HOME_DIR
fi
echo Home directory is $HOME_DIR

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
echo Directory path is $DIRECTORY
printf "\n"

# Sort chromosome regions

BP_RANGE=($(echo $BP_RANGE | perl -pne 's/([0-9]+)[,]*/$1 /g'))
IFS=$'\n'
BP_RANGE=($(echo "${BP_RANGE[*]}" | sort -gk1))
n_regions=$((${#BP_RANGE[@]}/2)) # number of regions on chromosome CHR
unset IFS

# Make array of n_regions copies of CHR
IFS=$' '
CHR=$(eval "printf ''$CHR' %.0s' {1..$n_regions}")
CHR=($(echo "${CHR[*]}"))
unset IFS

(( $((2*${#CHR[@]}))==${#BP_RANGE[@]} )) || (>&2 echo "Chromosome list and bp_range list do not agree.  Make sure each bp_range has a chromosome associated with it."; exit 1)

POPULATION=($(echo $POPULATION | perl -pne 's/[,]/ /g'))

test -z $VCF_FILE && (>&2 echo "VCF file not supplied"; exit 1)
test -f $VCF_FILE || (>&2 echo "VCF file $VCF_FILE not found"; exit 1)

echo Phased haplotypes files is $VCF_FILE
printf "\n"

if [ ! -f ${VCF_FILE}.csi ] && [ ! -f ${VCF_FILE}.tbi ] # create index if it does not exist
then
  echo Modify header:
  echo bcftools view -h ${VCF_FILE} \> ${VCF_FILE/.vcf.gz/.header.txt}
  echo awk \'\$0 \~ /\#\#FILTER/{\$0=sprintf\(\"%s\\n\#\#FILTER=\<ID=NA,Description=\\\"NA\\\"\>\",\$0\)}1\' ${VCF_FILE/.vcf.gz/.header.txt} \> ${VCF_FILE/.vcf.gz/.header_new.txt}
  echo bcftools reheader -h ${VCF_FILE/.vcf.gz/.header_new.txt} -o ${VCF_FILE/.vcf.gz/.new.vcf.gz} ${VCF_FILE}
  echo test -f ${VCF_FILE/.vcf.gz/.header.txt} \&\& rm ${VCF_FILE/.vcf.gz/.header.txt}
  echo test -f ${VCF_FILE/.vcf.gz/.header_new.txt} \&\& rm ${VCF_FILE/.vcf.gz/.header_new.txt}
  bcftools view -h ${VCF_FILE} > ${VCF_FILE/.vcf.gz/.header.txt}
  awk '$0 ~ /##FILTER/{$0=sprintf("%s\n##FILTER=<ID=NA,Description=\"NA\">",$0)}1' ${VCF_FILE/.vcf.gz/.header.txt} > ${VCF_FILE/.vcf.gz/.header_new.txt}
  bcftools reheader -h ${VCF_FILE/.vcf.gz/.header_new.txt} -o ${VCF_FILE/.vcf.gz/.new.vcf.gz} ${VCF_FILE}
  test -f ${VCF_FILE/.vcf.gz/.header.txt} && rm ${VCF_FILE/.vcf.gz/.header.txt}
  test -f ${VCF_FILE/.vcf.gz/.header_new.txt} && rm ${VCF_FILE/.vcf.gz/.header_new.txt}
  printf "\n"

  VCF_FILE=${VCF_FILE/.vcf.gz/.new.vcf.gz}
  echo Modified vcf file is $VCF_FILE
  printf "\n"

  echo Generate .tbi index file:
  echo tabix \-p vcf $VCF_FILE
  tabix -p vcf $VCF_FILE

  printf "\n"
fi

region_str=${CHR[i]}:${BP_RANGE[0]}\-${BP_RANGE[$((2*n_regions-1))]}
for ((i=1;i<${#CHR[@]};i++))
do
  region_str=${region_str},${CHR[i]}:${BP_RANGE[$((2*i))]}\-${BP_RANGE[$((2*i+1))]} # string for region specification in bcftools view
done
echo $region_str
printf "\n"

samples_file=$(echo ${POPULATION[*]%.*} | sed 's/ /+/g').samples
echo Get unique individuals in supplied .indiv file in order:
echo awk \'{seen[NR]=\$1\;} END{n = asorti\(seen,a\)\; for \(i=1\;i\<=n\;i++\) {print seen[a[i]]} }\' ${POPULATION[*]} \> $samples_file # unique samples in .indiv file for population i
awk '{seen[NR]=$1;} END{n = asorti(seen,a); for (i=1;i<=n;i++) {print seen[a[i]]} }' ${POPULATION[*]} > $samples_file # unique samples in .indiv file for population i, printed in same order as .indiv file
printf "\n"

transpose_file=haplotype_estimates_transpose.$(echo ${POPULATION[*]%.*} | sed 's/ /+/g').chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[$((2*n_regions-1))]}.txt
echo Get doubled header \(each subject id twice\):
echo bcftools view \-r ${region_str} \--samples-file $samples_file \--force-samples -m2 -M2 -v snps \--output-type v \--header-only $VCF_FILE \| awk \'{printf \"sjlife\\t\"\; for \(i=10\;i\<=NF\;i++\) {printf \"%s\\t%s%s\",\$i,\$i,\(i\<NF?\"\\t\":\"\\n\"\)} }\' \| tail \-1 \> $transpose_file
bcftools view -r ${region_str} --samples-file $samples_file --force-samples -m2 -M2 -v snps --output-type v --header-only $VCF_FILE | awk '{printf "sjlife\t"; for (i=10;i<=NF;i++) {printf "%s\t%s%s",$i,$i,(i<NF?"\t":"\n")} }' | tail -1 > $transpose_file
printf "\n"

new_vcf_file=$(echo ${POPULATION[*]%.*} | sed 's/ /+/g').chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[$((2*n_regions-1))]}.haplotypes.vcf
echo Extract individuals from .vcf:
echo bcftools view \-r ${region_str} \--samples-file $samples_file \--force-samples -m2 -M2 -v snps \--output-type v \--no-header \--output-file $new_vcf_file $VCF_FILE
bcftools view -r ${region_str} --samples-file $samples_file --force-samples -m2 -M2 -v snps --output-type v --no-header  --output-file $new_vcf_file $VCF_FILE
printf "\n"

if [ ! -z $INFO ] && [ -f $INFO ]
then
  echo Translate rsids:
  echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR \&\& \$2==${CHR[0]}{rsid[\$3]=\$6\; next} \(\$2 in rsid\){\$3=rsid[\$2]\; print \$0}\' $INFO $new_vcf_file \> ${new_vcf_file/.txt/.tmp}
  awk 'BEGIN{OFS="\t"} NR==FNR && $2=='${CHR[0]}'{rsid[$3]=$6; next} ($2 in rsid){$3=rsid[$2]; print $0}' $INFO $new_vcf_file > ${new_vcf_file}.tmp && printf "\n"
  test -f ${new_vcf_file}.tmp && echo mv ${new_vcf_file}.tmp $new_vcf_file
  test -f ${new_vcf_file}.tmp && mv ${new_vcf_file}.tmp $new_vcf_file && printf "\n"
fi

echo Compress and index vcf file:
echo bcftools view \-r $region_str \--samples-file $samples_file \--force-samples -m2 -M2 -v snps \--output-type z \--output-file $(echo ${POPULATION[*]%.*} | sed 's/ /+/g').chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[$((2*n_regions-1))]}.vcf.gz $VCF_FILE
bcftools view -r $region_str --samples-file $samples_file --force-samples -m2 -M2 -v snps --output-type z --output-file $(echo ${POPULATION[*]%.*} | sed 's/ /+/g').chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[$((2*n_regions-1))]}.vcf.gz $VCF_FILE
printf "\n"

echo tabix \-p vcf $(echo ${POPULATION[*]%.*} | sed 's/ /+/g').chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[$((2*n_regions-1))]}.vcf.gz
tabix -p vcf $(echo ${POPULATION[*]%.*} | sed 's/ /+/g').chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[$((2*n_regions-1))]}.vcf.gz
printf "\n"

echo awk \'BEGIN{FS=\"\\t\"\; OFS=\"\\t\"} {rsid=gensub\(/\([^,]\*\),\(.\*\)/,\"\\\\1\",\"g\",\$3\)\; \$3=rsid}1\' $new_vcf_file \> ${new_vcf_file}.tmp # convert $3 from rsid,chr to $3=rsid
awk 'BEGIN{FS="\t"; OFS="\t"} {rsid=gensub(/([^,]*),(.*)/,"\\1","g",$3); $3=rsid}1' $new_vcf_file > ${new_vcf_file}.tmp # convert $3 from rsid,chr to $3=rsid
printf "\n"

test -f ${new_vcf_file}.tmp && echo mv ${new_vcf_file}.tmp $new_vcf_file && printf "\n"
test -f ${new_vcf_file}.tmp && mv ${new_vcf_file}.tmp $new_vcf_file

echo Remove SNPs that are not biallelic:
echo awk \'NR==FNR{seen[\$3]+=1\; next} \(\$3 in seen\){if \(seen[\$3]==1\) {print \$0} }\' $new_vcf_file $new_vcf_file \> ${new_vcf_file}.tmp
awk 'NR==FNR{seen[$3]+=1; next} ($3 in seen){if (seen[$3]==1) {print $0} }' $new_vcf_file $new_vcf_file > ${new_vcf_file}.tmp
printf "\n"

test -f ${new_vcf_file}.tmp && echo mv ${new_vcf_file}.tmp $new_vcf_file && printf "\n"
test -f ${new_vcf_file}.tmp && mv ${new_vcf_file}.tmp $new_vcf_file

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

haplotypes_file=haplotype_estimates.$(echo ${POPULATION[*]%.*} | sed 's/ /+/g').chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[$((2*n_regions-1))]}.txt
echo Transpose haplotypes into chromosome x rsid:
echo awk \'BEGIN{OFS=\"\\t\"}\; {for\(j=1\;j\<=NF\;j++\) {a[NR,j]=\$j\; n_rows=NR\; n_cols=\(n_cols\<NF?NF:n_cols\)} } END{for \(j=1\;j\<=n_cols\;j++\) {for \(i=1\;i\<=n_rows\;i++\) {printf \"%s%s\",a[i,j],\(i==n_rows?\"\\n\":\"\\t\"\)} } }\' $transpose_file \> $haplotypes_file
awk 'BEGIN{OFS="\t"}; {for(j=1;j<=NF;j++) {a[NR,j]=$j; n_rows=NR; n_cols=(n_cols<NF?NF:n_cols)} } END{for (j=1;j<=n_cols;j++) {for (i=1;i<=n_rows;i++) {printf "%s%s",a[i,j],(i==n_rows?"\n":"\t")} } }' $transpose_file > $haplotypes_file
printf "\n"

if ((${#POPULATION[@]}>1))
then
  for ((i=0;i<${#POPULATION[@]};i++))
  do
    echo Extract ${POPULATION[i]} subjects as chromosome x rsid:
    echo awk \'BEGIN{OFS=\"\\t\"} NR==FNR{id[\$1]\; next} \(\$1 in id \|\| FNR==1\){print \$0}\' ${POPULATION[i]} $haplotypes_file \> haplotype_estimates.${POPULATION[i]%.*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[$((2*n_regions-1))]}.txt
    awk 'BEGIN{OFS="\t"} NR==FNR{id[$1]; next} ($1 in id || FNR==1){print $0}' ${POPULATION[i]} $haplotypes_file > haplotype_estimates.${POPULATION[i]%.*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[$((2*n_regions-1))]}.txt
    printf "\n"

    echo Transpose ${POPULATION[i]} haplotypes into rsid x chromosome:
    echo awk \'BEGIN{OFS=\"\\t\"}\; {for\(j=1\;j\<=NF\;j++\) {a[NR,j]=\$j\; n_rows=NR\; n_cols=\(n_cols\<NF?NF:n_cols\)} } END{for \(j=1\;j\<=n_cols\;j++\) {for \(i=1\;i\<=n_rows\;i++\) {printf \"%s%s\",a[i,j],\(i==n_rows?\"\\n\":\"\\t\"\)} } }\' haplotype_estimates.${POPULATION[i]%.*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[$((2*n_regions-1))]}.txt \> haplotype_estimates_transpose.${POPULATION[i]%.*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[$((2*n_regions-1))]}.txt
    awk 'BEGIN{OFS="\t"}; {for(j=1;j<=NF;j++) {a[NR,j]=$j; n_rows=NR; n_cols=(n_cols<NF?NF:n_cols)} } END{for (j=1;j<=n_cols;j++) {for (i=1;i<=n_rows;i++) {printf "%s%s",a[i,j],(i==n_rows?"\n":"\t")} } }' haplotype_estimates.${POPULATION[i]%.*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[$((2*n_regions-1))]}.txt > haplotype_estimates_transpose.${POPULATION[i]%.*}.chr${CHR[0]}.${BP_RANGE[0]}-${BP_RANGE[$((2*n_regions-1))]}.txt
    printf "\n"
  done
fi
