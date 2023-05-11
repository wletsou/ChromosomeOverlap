#! /bin/bash
set -e

# Forms all SIGMA-tuples of the rows of FILE and looks for common alleles.

# sh /home/wletsou/scripts/ChromosomeOverlap_iteration.v2.sh "001_0,002_0,003_0,005_0,006_0,008_1,009_0,010_0,012_0,014_0,016_0,017_1,018_0,019_0,021_0,022_0,023_0,024_0,026_0,029_0,030_0,031_0,034_0,037_0,041_0,043_0,046_0,047_0,050_0,051_0,054_0,055_0,057_0,060_0,061_0,063_0,064_0,068_0,069_0,070_0,072_0,073_0,074_0,075_0,079_0,081_0,082_0,084_0,085_1,089_0,090_0,095_0,096_0,097_0,098_0,099_0,103_0,104_0,105_0,107_0,108_0,110_0,111_0,112_0,113_0,114_0,116_0,117_0,118_0,120_0,122_0,124_0,125_0,127_0,129_0,130_0,131_0,132_0,134_0,135_0,136_0,140_0,143_0,145_0,147_0,149_0,150_0,152_0,153_0,154_0,155_0,156_0,159_0,160_0,161_0,163_0,164_0,167_0,169_0,170_0,172_0,173_0,178_0,180_0,182_0,183_0,186_0,188_0,189_0,190_0,192_0,195_0,196_0,198_0,199_0,202_0,205_0,206_0,210_0,213_0,215_0,216_0,217_0,218_0,220_0,221_0,222_0,223_0,224_0,225_0,227_0,229_0,231_0,232_0,234_0,236_0,237_0,238_0,240_0,241_0,242_0,243_0,245_0,246_0,247_0,249_0,250_0,252_0,253_0,255_0,258_0,259_0,260_0,261_0,263_0,264_0,265_0,269_0,271_0,272_0,275_0,276_0,277_0,280_0,281_0,282_0,284_0,285_0,288_1,289_0,292_0,293_0,294_0,295_0,297_1,298_0,299_0,301_0,302_1,303_0,304_0,305_1,306_0,309_0,311_0,313_0,314_0,317_1,319_1,320_0,321_0,325_0,326_0,327_0,328_0,329_0,330_0,332_0,333_0,334_0,335_0,338_0,342_0,343_0,344_0,345_0,347_0,348_0,350_0,351_0,352_0,353_0,354_0,355_0,356_0,357_0,358_0,360_0,361_0,362_0,363_0,364_0,366_0,367_0,369_0,372_0,373_0,374_0,377_0,380_0,381_0,384_0,385_0,387_0,388_0,389_0,390_0,391_0,393_0,394_0,395_0,399_0,400_0,401_0,402_0,403_0,404_0,405_0,406_0,408_0,409_0,410_0,411_0,414_0,415_0,416_0,417_0,419_0,420_0,421_0,422_0,423_0,427_0,429_0,430_0,431_0,432_0,433_0,437_0,439_0,440_0,441_0,442_0,444_0,445_0,448_0,449_0,450_0,451_0,453_0,454_0,456_0,457_0,458_0,459_0,460_0,461_0,462_0,466_0,468_0,470_0,471_0,473_0,475_0,476_0,479_0,481_0,482_0,483_0,484_0,485_0,486_0,487_0,488_0,489_0,490_0,491_0,492_0,493_0,494_0,495_0,496_0,498_0,500_0,502_0,503_0,506_0,508_0,509_0,510_0,511_0,512_0,513_0,516_0,518_0,520_0,524_0,526_0,527_0,528_0,529_0,531_0,535_0,539_0,543_0,546_0,547_0,548_0,551_0,552_0,553_0,554_0,556_0,559_0,560_0,562_0,563_0,564_0,567_0,568_0,570_0,572_0,573_0,574_0,575_0,579_0,580_0,581_0,582_0,584_0,586_0,587_0,588_0,589_0,590_0,591_0,594_0,596_0,597_0,599_0,600_0,603_0,604_0,605_1,606_0,607_0,609_0,610_0,612_0,614_1,616_0,618_1,621_0,622_0,624_0,625_0,627_1,628_0,629_0,630_0,631_0,632_0,633_0,634_0,635_0,636_0,637_0,638_0,639_0,640_0,641_0,642_0,643_0,644_0,645_0,646_0,647_0,648_0,649_0,651_0,653_0,656_0,657_0,663_0,664_0,666_0,670_0,671_0,672_0,673_0,674_0,675_0,676_0,679_0,682_0,684_0,685_0,686_1,687_0,688_0,689_1,690_0,692_1,693_1,697_0,701_0,702_0,703_0,704_0,705_0,706_0,707_0,709_0,711_0,712_1,713_1,714_1,716_1,717_0,719_1,720_1,721_0,722_0,723_1,724_0,726_0,727_0,728_0,730_0,731_0,732_1,733_0,734_1,735_0,736_1,737_0,738_1,739_1,740_1,741_1,742_0,745_0,746_0,748_0,749_0,750_0,751_0,752_0,753_0,754_0,755_1,756_0,757_0,761_0,762_1,763_0,765_0,766_0,767_1,768_0,770_0,772_0,773_0,775_0,776_0,778_0,779_0,780_0,782_0,783_0,784_0,786_0,787_0,789_0,790_0,791_0,792_0,793_0,794_0,795_0,796_0,799_0,800_0,801_0,803_0,804_0,806_0,807_0,808_0,809_0,811_0,813_0,814_0,816_0,818_0,819_0,820_0,821_0,822_0,823_0,824_0,825_0,830_0,831_0,832_0,833_0,834_0,835_0,837_0,839_0,840_0,841_0,842_0,843_0,844_0,845_0,846_0,847_0,848_0,849_0,850_0,851_0,853_0,854_0,855_0,856_0,858_0,861_0,862_0,864_0,865_0,869_0,872_0,873_0,875_0,877_0,879_0,882_0,883_0,884_0,885_0,886_0,887_0,888_0,892_0,894_0,895_0,896_0,898_0,899_0,900_0,904_0,905_0,906_0,912_0,914_0,917_0,919_0,922_0,925_0,929_0,930_0,931_0,933_0,936_0" Pattern_combined.Iteration000.chr8.128682464-129616034_2,j.txt "" 1000000 1 "h2" "/scratch_space/wletsou/sjlife/GWAS/CCS_chr8.3" "/home/wletsou/scripts"

HAPLOTYPE=$1 # haplotype to be included in every overlap
FILE=$2 # Pattern_combined file with first field counts and second and further fields patterns
TRANSPOSE_FILE=$3 # optional transposed version of FILE; searched for and/or created if empty
RANGE=$4 # Range (comma-separated list, lower,upper) of intersections (of one group of SIGMA rows with all others) to do, from 0 to (n_rows choose sigma) - 1.  OR a period-separated list "step_size.step_number" defining the range "(i-1) * step,i * step - 1"
STEP_SIZE=$5 # number of overlaps to do in one step
SIGMA=$6 # number of rows to overlap
NAME=$7 # optional tag for output files Overlap_tuples.NAME..., by iteration
DIRECTORY=$8 # output directory to store Output.Tuples.Iteration00x.Range files
HOME_DIR=$9 # location of program files

if [ -z $HOME_DIR ]
then
  unset HOME_DIR
fi

if [ -z $DIRECTORY ];
then
  DIRECTORY=$PWD
fi
echo Overlap files stored in $DIRECTORY
printf "\n"

module load R/3.6.1

if [ -z $SIGMA ];
then
  SIGMA=1 # number of rows to overlap with a given row
fi

[ -z $FILE ] && (>&2 echo "Pattern_combined file not supplied."; exit 1)
[ ! -f $FILE ] && (>&2 echo "Pattern_combined file not found."; exit 1)

if [ -z $NAME ]
then
  unset NAME
fi

n_rows=$(awk 'END{print NR}' $FILE) # number of patterns
num0=$n_rows
if (($num0>=$SIGMA))
then
  let num=1;
  let den=1
  for ((i=1; i<=$SIGMA; i++));
  do
    let num=$num*$((num0-i+1)) # error if negative
    let den=$den*$i
  done
  let n_tuples=($num/$den)
else
  n_tuples=0
fi
(($n_rows==0)) && n_tuples=1 # Correction for one way of arranging 0 things
echo $n_tuples total combination$((($n_tuples==1)) || echo "s").

if [ ! -z $RANGE ]; # user-provided comma-separated list lower,upper limits on combinations to perform in this round
then
  RANGE=($(echo $RANGE | perl -pne 's/([0-9]+)[.]+/$1 /g'))
  if ((${#RANGE[@]}>1)) # length=1 if TUPLE_RANGE is a comma-separated list
  then
    # For the case that RANGE is of the form STEP_SIZE.STEP
    ll=$(($((RANGE[1]-1))*${RANGE[0]})) # (i - 1) * step_size
    ul=$((RANGE[1]*RANGE[0]-1)) # i * step_size
  elif ((${#RANGE[@]}==2))
  then
    # For the case that RANGE is of the form LOWER,UPPER
    RANGE=($(echo $RANGE | perl -pne 's/([0-9]+)[,]+/$1 /g'))
    ll=${RANGE[0]} # starting tuple index
    ul=${RANGE[1]} # unding tuple index
  else
    (>&2 echo "Invalid range."; exit 1)
  fi
  (( $ll > $((n_tuples-1)) )) && ll=$((n_tuples-1))
  (( $ul > $((n_tuples-1)) )) && ul=$((n_tuples-1))
  (( $ll > $ul )) && (>&2 echo "Invalid range."; exit 1)
  RANGE=(${ll} ${ul})
else
  RANGE=(0 $((n_tuples-1))) # compute all combinations if none supplied
fi
echo Supplied range is \(${RANGE[0]},${RANGE[1]}\)
printf "\n"

((${RANGE[0]}>$((n_tuples-1)))) && (>&2 echo "Range exceeds total number of tuples."; exit 1)
((${RANGE[1]}>$((n_tuples-1)))) && RANGE[1]=$((n_tuples-1))
char=$(( $( echo $( (( $n_tuples>${RANGE[1]} )) && echo $n_tuples || echo ${RANGE[1]} ) | wc -m)-1)) # number of characters in the number of tuples
n_fields=$(awk 'BEGIN{x=0} (NF>x){x=NF} END{print x}' $FILE)

if [ -z $STEP_SIZE ]
then
  STEP_SIZE=$((${RANGE[1]}-${RANGE[0]}+1))
fi

if [ -z $TRANSPOSE_FILE ]
then
  TRANSPOSE_FILE=${FILE%.*}.transpose.${FILE##*.}
  exec 3<>$TRANSPOSE_FILE
  # (
  if flock -x -n 3
  then
    echo Obtained lock.
    echo wc \-l $TRANSPOSE_FILE
    if [ ! -s $TRANSPOSE_FILE ]
    then
      echo Transpose input file:
      echo awk \'BEGIN{OFS=\"\\t\"}\; {for\(j=1\;j\<=NF\;j++\) {a[NR,j]=\$j\; n_rows=NR\; n_cols=\(n_cols\<NF?NF:n_cols\)} } END{for \(j=1\;j\<=n_cols\;j++\) {for \(i=1\;i\<=n_rows\;i++\) {printf \"%s%s\",a[i,j],\(i==n_rows?\"\\n\":\"\\t\"\)} } }\' $FILE \> $TRANSPOSE_FILE # transpose to counts in first row, snp patterns by bar groups in subsequent rows
      awk 'BEGIN{OFS="\t"}; {for(j=1;j<=NF;j++) {a[NR,j]=$j; n_rows=NR; n_cols=(n_cols<NF?NF:n_cols)} } END{for (j=1;j<=n_cols;j++) {for (i=1;i<=n_rows;i++) {printf "%s%s",a[i,j],(i==n_rows?"\n":"\t")} } }' $FILE > $TRANSPOSE_FILE && printf "\n"
    fi
  else
    while ! flock -x -n 3
    do
      echo Did not obtain lock.
      sleep 0.1 && printf "\n"
    done
  fi
  #) 3<>$TRANSPOSE_FILE
  # 3>$TRANSPOSE_FILE # transpose file, but only once https://www.karltarvas.com/2021/05/17/bash-using-flock-to-ensure-parallel-scripts-perform-an-action-only-once.html
fi

echo wc \-l $TRANSPOSE_FILE
wc -l $TRANSPOSE_FILE && printf "\n"
echo awk \'{print NF}\' $TRANSPOSE_FILE
awk '{print NF}' $TRANSPOSE_FILE && printf "\n"

file_ID=${FILE%.*} # drop file extension
file_ID=${file_ID##*/} # delete to the left of the last /
echo File identifier is $file_ID
printf "\n"

test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && echo rm Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && printf "\n"
test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && rm Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
echo touch Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
touch Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt

if (($n_rows > 0))
then
  for ((k=1;k<=$(($((${RANGE[1]}-${RANGE[0]}+1+STEP_SIZE-1))/STEP_SIZE));k++)) # https://stackoverflow.com/questions/2395284/round-a-divided-number-in-bash
  do
    range[0]=$((${RANGE[0]}+(k-1)*STEP_SIZE))
    ((${range[0]}>${RANGE[1]})) && range[0]=${RANGE[1]}
    range[1]=$((${RANGE[0]}+k*STEP_SIZE-1))
    ((${range[1]}>${RANGE[1]})) && range[1]=${RANGE[1]}
    declare -p range
    echo Rscript ${HOME_DIR/%/\/}index2combo2.R I=${range[0]} n=$n_rows sigma=$SIGMA
    Rscript ${HOME_DIR/%/\/}index2combo2.R I=${range[0]} n=$n_rows sigma=$SIGMA
    out=($(Rscript ${HOME_DIR/%/\/}index2combo2.R I=${range[0]} n=$n_rows sigma=$SIGMA)) # get the tuple corresponding to combination ${RANGE[0]} in the list of combinations from 0 to (n_rows choose SIGMA) - 1
    printf "\n"

    echo First joining tuple \(of rows\) is:
    declare -p out
    printf "\n"

    str=$(echo "awk 'BEGIN{OFS=\"\\t\"; n0=split(\"$HAPLOTYPE\",array,\",\"); for (j=1;j<=n0;j++) {snps_0[array[j]]}; delete array; i0=1} NR>1{ k=${range[0]}; ")
    for ((i=1;i<${#out[@]};i++))
    do
      str=${str}$(echo "if (i$((i-1))>${out[$((i-2))]}) {ll=i$((i-1))+1} else {ll=${out[$((i-1))]}}; for (i$i=ll;i$i<=$n_rows;i$i++) {n$i=split(\$i$i,array,\",\"); delete snps_$i; for (j=1;j<=n$i;j++) {snps_$i[array[j]]}; delete array; ")
    done
    str=${str}$(echo "if (i$((${#out[@]}-1))>${out[$((${#out[@]}-2))]}) {ll=i$((${#out[@]}-1))+1} else {ll=${out[$((${#out[@]}-1))]}}; for (i${#out[@]}=ll;i${#out[@]}<=$n_rows;i${#out[@]}++) { n${#out[@]}=split(\$i${#out[@]},snps_${#out[@]},\",\"); idx=0; printf \"%s\t\", 1; for (j=1;j<=n${#out[@]};j++) { if (")
    for ((i=1;i<${#out[@]};i++))
    do
      if ((i>1))
      then
        str=$str$(echo " && ")
      fi
      str=$str$(echo "snps_${#out[@]}[j] in snps_$i")
    done
    if [ ! -z $HAPLOTYPE ] && (( $(awk 'BEGIN{n0=split("'$HAPLOTYPE'",array,","); print n0}')>0 )) # if a non-zero-length haplotype was supplied, then add it to the comparrison
    then
      if ((i>1))
      then
        str=$str$(echo " && ")
      fi
      str=$str$(echo "snps_${#out[@]}[j] in snps_0")
    fi
    str=${str}$(echo ") {printf \"%s%s\",(idx==0?\"\":\",\"),snps_${#out[@]}[j]; idx+=1 } ")
    for ((i=${#out[@]};i>1;i--))
    do
      str=${str}$(echo "} ")
    done
    if ((${#out[@]}==1))
    then
      str=${str}$(echo "} ")
    fi
    str=${str}$(echo "if (idx==0) {printf \"0\"}; k+=1; printf \"\n\"; if (k>${range[1]}) {next} ")
    for ((i=${#out[@]};i>=1;i--))
    do
      str=${str}$(echo "} ")
    done
    str=${str}$(echo "}' $TRANSPOSE_FILE >> Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt")
    echo $str
    printf "\n"
    start=$(date +%s.%N)
    eval $str
    finish=$(date +%s.%N)
    duration=$(awk 'BEGIN{printf "%0.2f\n", '$finish'-'$start'}')
    echo Overlap took $duration second$( [ "$duration" == "1.00" ] && echo "" || echo "s").
    printf "\n"
  done
fi

if (($n_fields>2))
then
  # Make a tab-delimined list of columns 2,3,...,n_fields and split into an array. Sort the array and print
  echo Sort fields of Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt:
  str=$(echo "awk 'BEGIN{ OFS=\"\t\" } { split(")
  for ((i=2; i<$n_fields; i++));
  do
    str=${str}$(echo "\$$i\"\t\"")
  done
  str=${str}$(echo "\$$n_fields,array,\"\t\"); asort(array); { print \$1,")
  for ((i=1; i<$(($n_fields-1)); i++));
  do
    str=${str}$(echo "array[$i],")
  done
  str=${str}$(echo "array[$(($n_fields-1))] } }' Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt > Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp ")
  echo $str
  start=$(date +%s.%N)
  printf "\n"
  output=$(eval "$str" 2> /dev/null) # eliminate error message
  finish=$(date +%s.%N)
  duration=$(awk 'BEGIN{printf "%0.2f\n", '$finish'-'$start'}')
  echo Sorting fields took $duration second$( [ "$duration" == "1.00" ] && echo "" || echo "s").
  if (( $(echo $output | wc -w) != 0 )) # check if output is 0-length (i.e., an error)
  then
    echo cat Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt > Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp
    #cat Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt > Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp
  fi
  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && echo mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && printf "\n"
  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
fi

while ! test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
do
  echo Waiting for Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt...
  sleep 10
done

# ! test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && (&>2 echo "Could not find Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt"; exit 1)

# Sort SNPs within a field
sort_snps=0
if (($sort_snps==1))
then
  echo Sort SNPs within each group:
  start=$(date +%s.%N)
  str=$(echo "awk 'BEGIN{ OFS=\"\t\" } { split(\$2,array2,\",\"); n=asort(array2); idx=0; for (j=1;j<=n;j++) {if (idx==0) {\$2=array2[j]; idx++} else {\$2=\$2\",\"array2[j]} } }")
  for ((i=3; i<=$n_fields; i++));
  do
    str=${str}$(echo " { split(\$$i,array$i,\",\"); n=asort(array$i); idx=0; for (j=1;j<=n;j++) {if (idx==0) {\$$i=array$i[j]; idx++} else {\$$i=\$$i\",\"array$i[j]} } }")
  done
  str=${str}$(echo " { print }' Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt > Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp")
  echo $str
  eval $str
  finish=$(date +%s.%N)
  duration=$(awk 'BEGIN{printf "%0.2f\n", '$finish'-'$start'}')
  echo Sorting SNPs took $duration second$( [ "$duration" == "1.00" ] && echo "" || echo "s").
  printf "\n"
  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && echo mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && printf "\n"
  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
fi

while ! test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
do
  echo Waiting for Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt...
  sleep 10
done

if (($n_fields>2))
then
  start=$(date +%s.%N)
  echo \"Slide over\" fields to the right into empty fields to the left:
  echo awk \'BEGIN{ OFS=\"\\t\" } { for \(i=1\; i\<NF\; i++\) { if \(\$i==\"\"\) { \$i=\$\(i+1\)\; \$\(i+1\)=\"\" } } } { print }\' Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt \> Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp
  awk 'BEGIN{ OFS="\t" } { for (i=1; i<NF; i++) { if ($i=="") { $i=$(i+1); $(i+1)="" } } } { print }' Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt > Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp
  printf "\n"
  finish=$(date +%s.%N)
  duration=$(awk 'BEGIN{printf "%0.2f\n", '$finish'-'$start'}')
  echo Sliding over took $duration second$( [ "$duration" == "1.00" ] && echo "" || echo "s").
  printf "\n"
  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && echo mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && printf "\n"
fi

while ! test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
do
  echo Waiting for Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt...
  sleep 10
done

# Sum the first column of rows with matching patterns in all of columns 2,3,...,n_fields
echo Sum counts of matching rows:
start=$(date +%s.%N)
str=$(echo "awk 'BEGIN{OFS=\"\t\"} { array[")
for ((i=2; i<$n_fields; i++));
do
  str=${str}$(echo "\$$i\"\t\"")
done
str=${str}$(echo "\$$n_fields]+=\$1 } END{ for (row in array) { print array[row],row } }' Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt > Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp")
echo $str
eval "$str"
finish=$(date +%s.%N)
duration=$(awk 'BEGIN{printf "%0.2f\n", '$finish'-'$start'}')
echo Summing took $duration second$( [ "$duration" == "1.00" ] && echo "" || echo "s").
printf "\n"

while ! test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp
do
  echo Waiting for Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp...
  sleep 10
done

test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && echo mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && printf "\n"
test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt

while ! test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
do
  echo Waiting for Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt...
  sleep 10
done

# eliminate null rows
check_zeros=1
if (($check_zeros==1))
then
  echo Eliminate rows with all 0\'s:
  start=$(date +%s.%N)
  str=$(echo "awk 'BEGIN{OFS=\"\t\"} { if (")
  for ((i=2; i<$n_fields; i++))
  do
    str=${str}$(echo "\$$i != 0 &&")
  done
  str=${str}$(echo "\$$n_fields != 0) { print \$0 } }' Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt > Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp")
  echo $str
  eval "$str"
  finish=$(date +%s.%N)
  duration=$(awk 'BEGIN{printf "%0.2f\n", '$finish'-'$start'}')
  echo Elimination took $duration second$( [ "$duration" == "1.00" ] && echo "" || echo "s").
  printf "\n"

  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && echo mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && printf "\n"
  test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp && mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.tmp Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
fi

while ! test -f Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt
do
  echo Waiting for Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt...
  sleep 10
done

# Move to output directory with "temp" dropped from name
f_name=$(printf "Overlap_tuples${NAME/#/.}.%0.${char}d-%0.${char}d.%s.txt" ${RANGE[0]} ${RANGE[1]} ${file_ID})
# check if error during mv operation
if [[ "${PWD/%/\/}Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt" != "${DIRECTORY/%/\/}${f_name}" ]]
then
  test -f  Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && echo mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt ${DIRECTORY/%/\/}${f_name} && printf "\n"
  test -f  Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt && mv Overlap_tuples.${RANGE[0]}-${RANGE[1]}.${file_ID}.txt ${DIRECTORY/%/\/}${f_name}
fi

while ! test -f ${DIRECTORY/%/\/}${f_name}
do
  echo Waiting for ${DIRECTORY/%/\/}${f_name}...
  sleep 10
done
