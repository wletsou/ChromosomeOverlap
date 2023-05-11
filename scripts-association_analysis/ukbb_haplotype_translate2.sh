#! /bin/bash
set -e

POPULATION=$1 # haplotype_estimates file for population with subject id in first column and rsids in columns 2 and beyond
HAPLOTYPES=$2 # colon-separated list of haplotypes, or a file with haplotypes in the first column.  Can be of the form column_allele,... or rsid_allele=[0/1]
RANGE=$3 # Range (ll,ul) of haplotypes to do, from 0 to n_haplotypes - 1.  OR a period-separated list "STEP_SIZE.STEP_NO" defining the range "(i-1) * step,i * step - 1"
DIRECTORY=$4 # current working directory
HOME_DIR=$5 # location of scripts to be run

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
char_s=$(($(awk 'NR==1{print NF-1}' ${POPULATION[0]} | wc -c)-1)) # number of characters in the total number snps

if [ -f $HAPLOTYPES ]
then
  HAPLOTYPES=($(cut -f1 $HAPLOTYPES)) # first column of file, assuming no header
else
  HAPLOTYPES=($(echo $HAPLOTYPES | perl -pne 's/[:]/ /g'))
fi
n_haplotypes=${#HAPLOTYPES[@]}
char_h=$(( $(echo $n_haplotypes | wc -c)-1)) # number of characters in the total number haplotypes

if [ -z $RANGE ]
then
  RANGE="${n_haplotypes}.1"
fi

RANGE=($(echo $RANGE | perl -pne 's/([0-9]+)[.]+/$1 /g'))
if ((${#RANGE[@]}>1)) # length=1 if RANGE is a comma-separated list, =2 if in the form step_size.step_no
then
  ll=$(($((RANGE[1]-1))*$((RANGE[0])))) # (i - 1) * step_size
  ul=$(($((RANGE[1]*RANGE[0]))-1)) # i * step_size
  RANGE=$(echo ${ll},${ul})
fi

# two-element array of the lower and upper ranges of sigma0-tuples to sample
RANGE=($(echo $RANGE | perl -pne 's/([0-9]+)[,]+/$1 /g'))
echo Supplied range is \(${RANGE[0]},${RANGE[1]}\).
printf "\n"

((${RANGE[0]}>$((n_haplotypes-1)))) && (>&2 echo "Range exceeds total number of haplotypes."; exit 1)
((${RANGE[1]}>$((n_haplotypes-1)))) && RANGE[1]=$((n_haplotypes-1))
STEP_SIZE=$((${RANGE[1]}-${RANGE[0]}+1))
echo Corrected range is \(${RANGE[0]},${RANGE[1]}\).
printf "\n"

test -f ${2%.*}.$(eval "printf '%0${char_h}d-%0${char_h}d' ${RANGE[0]} ${RANGE[1]}").translated.txt && rm ${2%.*}.$(eval "printf '%0${char_h}d-%0${char_h}d' ${RANGE[0]} ${RANGE[1]}").translated.txt && echo rm ${2%.*}.$(eval "printf '%0${char_h}d-%0${char_h}d' ${RANGE[0]} ${RANGE[1]}").translated.txt && printf "\n"
test -f $2 && (touch ${2%.*}.$(eval "printf '%0${char_h}d-%0${char_h}d' ${RANGE[0]} ${RANGE[1]}").translated.txt && echo touch ${2%.*}.$(eval "printf '%0${char_h}d-%0${char_h}d' ${RANGE[0]} ${RANGE[1]}").translated.txt && printf "\n") || echo ""
for ((i=${RANGE[0]};i<=${RANGE[1]};i++))
do
  echo Haplotype $((i+1)) of ${STEP_SIZE}:
  if [[ ${HAPLOTYPES[i]} =~ ^[0-9].*_[0-9]* ]] # haplotype as a list of columns (and values) corresponding to fields of supplied haplotype_estimates file
  then
    TRANSLATED[i]=$(awk 'BEGIN{n=split("'${HAPLOTYPES[i]}'",array,","); for (i=1;i<=n;i++) {b=gensub("0*([1-9][0-9]*).*","\\1","g",array[i]); a=gensub("[0-9]*_([0-9]*)","\\1","g",array[i]); snps[i]=b+1; alleles[i]=a;} } NR==1{for (i=1;i<=n;i++) {printf "%s%s=%s%s",(i>1?",":""),$snps[i],alleles[i],(i==n?"\n":"")} }' ${POPULATION[0]}) # translate to rsid_allele=[0/1] list
  elif [[ ${HAPLOTYPES[i]} =~ ^rs.*_.*=.* ]] # haplotype as a list of rsids (with values) corresponding to fields of supplied haplotype_estimates file
  then
    TRANSLATED[i]=$(awk 'BEGIN{n=split("'${HAPLOTYPES[i]}'",array,","); for (i=1;i<=n;i++) {b=gensub("(rs[0-9]*_[A-Z]*)=[0-9]*","\\1","g",array[i]); a=gensub("rs[0-9]*_[A-Z]*=([0-9]*)","\\1","g",array[i]); snps[b]=a} } NR==1{ for (i=2;i<=NF;i++) {if ($i in snps) {found+=1; printf "%s%0.'${char_s}'d_%s%s",(found>1?",":""),i-1,snps[$i],(found<n?"":"\n")} } }' ${POPULATION[0]}) # translate to column_allele[0/1] list
  fi
  if [ -f $2 ] # if haplotypes input is a file
  then
    # file with haplotype pattern translated
    awk 'BEGIN{OFS="\t"} ($1=="'${HAPLOTYPES[i]}'"){$1="'${TRANSLATED[i]}'"; print $0}' $2 >> ${2%.*}.$(eval "printf '%0${char_h}d-%0${char_h}d' ${RANGE[0]} ${RANGE[1]}").translated.txt
  else
    echo ${TRANSLATED[i]} # print result if not saving to file
  fi
done
