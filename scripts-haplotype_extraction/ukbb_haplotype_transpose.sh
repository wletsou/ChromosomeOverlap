#! /bin/bash
set -e

# transpose a haplotype_estimates file from chromosomes x SNP to SNP x chromosome, or vice versa

# bsub -P SJLIFE -J ukbb_haplotype_transpose -oo ukbb_haplotype_transpose.out -eo ukbb_haplotype_transpose.err -R "rusage[mem=10000]" -q "large_mem" "/home/wletsou/scripts/ukbb_haplotype_transpose.sh haplotype_estimates.ukbb_bca_cases+controls.chr11.69231642-69431642.txt 1"

FILE=$1 # haplotypes file to be transposed into rsid x chromosome
HEADER=$2 # whether (1) or not (0) the transposed file should have a header
DIRECTORY=$3 # working directory

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
echo Directory path is $DIRECTORY
printf "\n"

filename=${FILE##*/} # trim directory up to rightmost forward slash "/"
filename_temp=${filename/estimates./estimates_transpose.}
if [ "$filename_temp" == "$filename" ] # file is already transposed
then
  filename=${filename/estimates_transpose./estimates.}
else
  filename=${filename/estimates./estimates_transpose.}
fi
echo Output file name is $filename.
printf "\n"

test -z $HEADER && HEADER=1 # include header by default
! [[ $HEADER =~ ^[01]+$ ]] && (>&2 echo "Error: HEADER not 0 or 1"; exit 1)
echo Header will$( (($HEADER==0)) && echo " not" || echo "") be included.
printf "\n"

start=$(date +%s) # starting time in seconds
echo Transpose haplotypes into rsid x chromosome:
echo awk \'BEGIN{OFS=\"\\t\"}\; {for\(j=1\;j\<=NF\;j++\) {a[NR,j]=\$j\; n_rows=NR\; n_cols=\(n_cols\<NF?NF:n_cols\)} } END{for \(j=1\;j\<=n_cols\;j++\) {for \(i=1\;i\<=n_rows\;i++\) {printf \"%s%s\",a[i,j],\(i==n_rows?\"\\n\":\"\\t\"\)} } }\' ${FILE} \> ${DIRECTORY}/${filename}
awk 'BEGIN{OFS="\t"}; {for(j=1+'$((1-HEADER))';j<=NF;j++) {a[NR,j]=$j; n_rows=NR; n_cols=(n_cols<NF?NF:n_cols)} } END{for (j=1+'$((1-HEADER))';j<=n_cols;j++) {for (i=1;i<=n_rows;i++) {printf "%s%s",a[i,j],(i==n_rows?"\n":"\t")} } }' ${FILE} > ${DIRECTORY}/${filename}
printf "\n"
duration=$(( $(date +%s)-start )) # time finished less time started, in seconds
echo Transposition completed in $duration second$( (($duration==1)) && echo "" || echo "s" ).
printf "\n"
