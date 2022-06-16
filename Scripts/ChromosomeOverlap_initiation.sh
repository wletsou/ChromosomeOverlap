#! /bin/bash

# Overlaps the variable chromosome(s) given by the column(s) INDEX of the (transposed) HAPLOTYPES_FILE with every column j of the HAPLOTYPES file which is greater than the maximum column in INDEX

HAPLOTYPES_FILE=$1 # Matrix of snp genotypes (rows) x subjects (columns), no header, but first column is rsids
INDEX=$2 # comma-separated list i1,...i-(sigma-1),j, where i are the fixed chromosomes and j is the variable chromosome; set j to the maximum chromosome index and it will be replaced; example: 2,100 overlaps column 2 with every column up to 100; 2,3,100 overlaps columns 2 and 3 with every column
TRANSPOSE=$3 # whether to transpose input file from subjects (rows) x genotypes (columns); set to any non-null value

file_ID=${HAPLOTYPES_FILE%.*} #get identifiers from HAPLOTYPES_FILE before the file extension
echo fileID is ${file_ID}
printf "\n"
if [ ! -z $TRANSPOSE ];
then
  n_fields=$(awk 'END{print NF}' $HAPLOTYPES_FILE)
  for i in `seq 2 $n_fields`
  do
    cut -f$i $HAPLOTYPES_FILE | paste -s # paste column as a row
  done > subjects_${file_ID}_${INDEX}.txt
  n_subjects=$(awk 'END{print NF-1}' subjects_${file_ID}_${INDEX}.txt)
  n_snps=$(awk 'END{print NR}' subjects_${file_ID}_${INDEX}.txt)
else
  n_subjects=$(awk 'END{print NF}' $HAPLOTYPES_FILE)
  n_snps=$(awk 'END{print NR}' $HAPLOTYPES_FILE)
fi

# Generate all combinations of the INDEX set separated by either a comma (,) or a bar (|)
list=$(echo $INDEX | perl -pne 's/([0-9]+)[,]*/$1 /g')
index_list=($list)
len=${#index_list[@]}
str=$(perl -e "print '{$INDEX}{\",\",\"|\"}'x$len")
str=$(echo $str | perl -pne 's/{[^{]+}$//g')

eval str=($str)
for partition in "${str[@]}"
do
  echo $partition
done > stirling_${file_ID}_${INDEX}.txt

echo List of possible partitionings of $INDEX:
cat stirling_${file_ID}_${INDEX}.txt
printf "\n"

# Prune the combinations that are not increasing by and within each bar group; e.g. keep 1,3|2 but not 3,1|2 (uses lexigraphical ordering, 100 < 11)
code=$(echo stirling_${file_ID}_${INDEX}.txt | awk '{printf "awk -F\",\" \x27 BEGIN{OFS=\",\"} { for (i=1;i<NF;i++) if (gensub(\"[|][0-9]+\",\"\",\"g\",$(i+1))<=gensub(\"[0-9]+[|]\",\"\",\"g\",$i)) $(i+1)=\"\"} {print}\x27 %s | awk -F\"|\" \x27 BEGIN{OFS=\"|\"} { for (i=1;i<NF;i++) if (gensub(\",[0-9]+\",\"\",\"g\",$(i+1))<=gensub(\",[0-9]+\",\"\",\"g\",$i)) $(i+1)=\"\"} {print}\x27 > stirling_'${file_ID}'_'${INDEX}'_temp.txt",$0}')
echo $code
printf "\n"
eval "$code"

# Index values in the INDEX set by their position
INDEX_ORDER=$(seq 1 $len)
INDEX_ORDER=$(echo $INDEX_ORDER | perl -pne 's/[ ]+/,/g')
declare -p INDEX_ORDER

echo Prune combinations that do not sum to INDEX total:
expr=$(echo $INDEX_ORDER | perl -pne 's/([0-9]+)[,]*/\$$1+/g')
expr=$(echo $expr | perl -pne 's/[+]$//g')
val=$(echo $INDEX | perl -pne 's/([0-9]+)[,]*/$1+/g')
val=$(echo $val | perl -pne 's/[+]$//g')
let val=$val
str=$(echo ${expr}==${val})
code1=$(echo $str | awk '{printf "awk -F\"[,|]\" \x27%s {print}\x27 stirling_'${file_ID}'_'${INDEX}'_temp.txt > stirling_'${file_ID}'_'${INDEX}'_temp1.txt",$0}')
echo $code1
eval "$code1"
printf "\n"
cat stirling_${file_ID}_${INDEX}_temp1.txt
printf "\n"

echo Prune combinations that do not multiply to INDEX product:
expr=$(echo $INDEX_ORDER | perl -pne 's/([0-9]+)[,]*/\$$1*/g')
expr=$(echo $expr | perl -pne 's/[*]$//g')
val=$(echo $INDEX | perl -pne 's/([0-9]+)[,]*/$1*/g')
val=$(echo $val | perl -pne 's/[*]$//g')
let val=$val
str=$(echo ${expr}==${val})
code2=$(echo $str | awk '{printf "awk -F\"[,|]\" \x27%s {print}\x27 stirling_'${file_ID}'_'${INDEX}'_temp1.txt > stirling_'${file_ID}'_'${INDEX}'_subsets.txt",$0}')
echo $code2
eval "$code2"
printf "\n"
cat stirling_${file_ID}_${INDEX}_subsets.txt
printf "\n"

n_combos=$(awk 'END{print NR}' stirling_${file_ID}_${INDEX}_subsets.txt)
INDEX_SEP=$(echo $INDEX | perl -pne 's/,/ /g')
# Find the maximum index in the input; it will be replaced in combinations
max=0;
for i in ${INDEX_SEP[*]}
do
  (( i > max )) && max=$i
done
echo Maximum index is $max.
printf "\n"
INDEX_NEW=$(echo $INDEX | perl -pne "s/$max/j/g") # replace maximum index with j to be joined

# Join groups of subjects and look for unique snp patterns
for row in `seq 1 $n_combos` # get stirling combination from file
do

  code_str=$(echo "awk 'NR==$row{print \$0}' stirling_${file_ID}_${INDEX}_subsets.txt")
  echo $code_str
  printf "\n"
  pattern=$(eval "$code_str")
  pattern=$(echo $pattern | perl -pne "s/[|]/+/g")
  echo Row $row from stirling_${file_ID}_${INDEX}_subsets.txt is $pattern
  code_str=$(echo "awk 'NR==$row{print \$0}' stirling_${file_ID}_${INDEX}_subsets.txt > combo_${file_ID}_${INDEX}_${pattern}_temp.txt")
  echo $code_str
  eval "$code_str"
  printf "\n"

  COMBO=$(awk -F\"\t\" "BEGIN{
    OFS=\"|\"
  } {
    split(\$0,array,\"|\");
    for (i in array) {
      split(array[i],a,\",\");
      asort(a);
      idx=0;
      for (j in a) {
        if (idx==0) {
          \$i=a[j];
          idx++
        } else {
          \$i=\$i\",\"a[j]
        }
      }
    }
  } {
    print
  }" combo_${file_ID}_${INDEX}_${pattern}_temp.txt) # sort numerically within bar group
  echo $COMBO > combo_${file_ID}_${INDEX}_${pattern}_temp.txt
  COMBO=$(awk -F\"\t\" "BEGIN{
    OFS=\"|\"
  } {
    split(\$0,array,\"|\");
    for (i in array) {
      split(array[i],a,\",\");
      b[sprintf(\"%0.8d\",a[1])]=array[i]
    }
    asorti(b,c);
    col=1;
    for (i in c) {
      \$col=b[c[i]];
      col++
    }
  } {
    print
  }" combo_${file_ID}_${INDEX}_${pattern}_temp.txt) # sort the bar groups based on the value of the first comma-separated entry

  echo Sorted combo is $COMBO
  COMBO_NEW=$(echo $COMBO | perl -pne "s/$max/j/g") # replace the maximum index with j to be joined
  COMBO_file=$(echo $COMBO_NEW | perl -pne "s/[|]/+/g") # replace bars | with plus + for file naming

  BAR_COMBO=$(echo $COMBO | perl -pne 's/\|/ /g') # separate by bars |
  for partition in $BAR_COMBO # break combination into groups by bars |
  do
    echo Bar partition is $partition
    partition_new=$(echo $partition | perl -pne "s/$max/j/g") # has a j instead of max
    echo partion_new=$partition_new
    partition_file=$(echo $partition_new | perl -pne "s/[|]/+/g") # replace bars | with + for file naming
    echo partition_file=$partition_file
    PARTITION_SEP=$(echo $partition | perl -pne 's/,/ /g') # bar group without commas
    PARTITION_SEP=$(echo $PARTITION_SEP | perl -pne "s/$max//g") # remove commas and repace max with j, sorted as array; max=j is deleted
    echo PARTITION_SEP=$PARTITION_SEP
    array=($PARTITION_SEP)
    declare -p array
    printf "\n"

    MAXIND=0
    for ((i=0;i<${#array[@]};i++));
    do
      if ((${array[i]} > $MAXIND));
      then
        MAXIND=${array[i]}
      fi
    done
    echo Maximum index in list is $MAXIND
    PARTITION_JOIN=$(echo $PARTITION_SEP | perl -pne 's/ /,/g') # put commas back in for file naming
    echo PARTITION_JOIN=$PARTITION_JOIN
    printf "\n"

    if [ -f subjects_${file_ID}_${INDEX}.txt ];
    then
      echo cp \--remove-destination subjects_${file_ID}_${INDEX}.txt subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
      cp --remove-destination subjects_${file_ID}_${INDEX}.txt subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
    else
      echo cp \--remove-destination $HAPLOTYPES_FILE subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
      cp --remove-destination $HAPLOTYPES_FILE subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
    fi
    printf "\n"
    echo Index combination for the current bar group is $PARTITION_JOIN.  Checking to see if it combines with j.
    echo Joined set is $partition_new
    printf "\n"
    # After replacing the max index with j, perform the joining operation at each column j
    if [ ! -z $PARTITION_JOIN ] && [ $PARTITION_JOIN != $partition_new ];
    then
      echo Join index with j
      echo Join $PARTITION_SEP with j
      printf "\n"
      # get the values of the index columns and prepend them to every column j exceeding the maximum index.  Then assign NA to all others, including the index set
      awk -F"\t" "BEGIN{
        OFS=\"\t\";
      } {
        split(\"$PARTITION_SEP\",array,\" \");
        n=asort(array)
        val[\$1]=\$array[1]
        for (i = 2; i <= n; i++) {
          val[\$1]=val[\$1]\$array[i]
        } for (i = 2; i <= NF; i++) {
          if (i > $MAXIND) {
            \$i=val[\$1]\$i
          } else {
            \$i=\"NA\"
          }
        } {
          print
        }
      }" subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt > subjects_${file_ID}_${INDEX}_temp.txt

      echo write file subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
      echo cp \--remove-destination subjects_${file_ID}_${INDEX}_temp.txt subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
      cp --remove-destination subjects_${file_ID}_${INDEX}_temp.txt subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
      printf "\n"
      if [ -f subjects_${file_ID}_${INDEX}_temp.txt ];
      then
        echo rm subjects_${file_ID}_${INDEX}_temp.txt
        rm subjects_${file_ID}_${INDEX}_temp.txt
        printf "\n"
      fi
      #done
    elif [ ! -z $PARTITION_JOIN ] && [ $PARTITION_JOIN == $partition_new ]; #index is not joined to j
    then
      echo make $n_subjects copies of $PARTITION_JOIN
      printf "\n"
      # all columns are the same, i.e., the non-j index columns merged
      awk -F"\t" "BEGIN{
        OFS=\"\t\";
      } {
        split(\"$PARTITION_SEP\",array,\" \");
        n=asort(array)
        val[\$1]=\$array[1]
        for (i = 2; i <= n; i++) {
          val[\$1]=val[\$1]\$array[i]
        } for (i = 2; i <= NF; i++) {
          if (i > $MAXIND) {
            \$i=val[\$1]
          } else {
            \$i=\"NA\"
          }
        } {
          print
        }
      }" subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt > subjects_${file_ID}_${INDEX}_temp.txt

      echo Write file subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
      echo cp \--remove-destination subjects_${file_ID}_${INDEX}_temp.txt subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
      cp --remove-destination subjects_${file_ID}_${INDEX}_temp.txt subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
      printf "\n"
      if [ -f subjects_${file_ID}_${INDEX}_temp.txt ];
      then
        echo rm subjects_${file_ID}_${INDEX}_temp.txt
        rm subjects_${file_ID}_${INDEX}_temp.txt
        printf "\n"
      fi

      index_array=($PARTITION_SEP)

    elif [ -z $PARTITION_JOIN ];
    then
      echo Evaluate for j alone:
      echo cp \--remove-destination subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt subjects_${file_ID}_${INDEX}_temp.txt
      cp --remove-destination subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt subjects_${file_ID}_${INDEX}_temp.txt
      # echo awk \'BEGIN{OFS=\"\\t\"} {for \(i=2\;i\<=NF\;i++\) {if \(i \<= $MAXIND\) {\$i=\"NA\"} }\; print \$0}\' subjects_${file_ID}_${INDEX}_temp.txt \> subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
      # awk 'BEGIN{OFS="\t"} {for (i=2;i<=NF;i++) {if (i <= '$MAXIND') {$i="NA"} }; print $0}' subjects_${file_ID}_${INDEX}_temp.txt > subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt && printf "\n"
      echo write file subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
      cp --remove-destination subjects_${file_ID}_${INDEX}_temp.txt subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
      printf "\n"
      if [ -f subjects_${file_ID}_${INDEX}_temp.txt ];
      then
        echo rm subjects_${file_ID}_${INDEX}_temp.txt
        rm subjects_${file_ID}_${INDEX}_temp.txt
        printf "\n"
      fi
    fi
    if [ -f subjects_${file_ID}_${INDEX}_temp.txt ];
    then
      echo rm -f subjects_${file_ID}_${INDEX}_temp.txt
      rm -f subjects_${file_ID}_${INDEX}_temp.txt
      printf "\n"
    fi
    echo What does this file subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt look like?
    head subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
    printf "\n"

    char=$(($(cat subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt | wc -l | wc -m)-1)) # number of characters to use to print each SNP index
    # Look in each column for the lines corresponding to the unique snp sharing patterns, e.g., 210 at line 3 means snp 3 is homozygous in subject 1, heterozygous in 2, and absent in 3. Use gensub to remove leading commas.  Outputs a three-column table of chromosome index j, pattern type, and lines with the pattern.  Each j occurs once for every pattern it participates in.  The line numbers form the SNP patterns used downstream
    echo awk \-F\"\\t\" \'BEGIN{OFS=\"\\t\"} NR\>=1{ for \(i=2\;i\<= NF\;i++\) {snps[i\,\$i]=sprintf\(\"%s\,%0.${char}d\"\,snps[i\,\$i]\,NR\)} } END{for \(i in snps\) {split\(i,idx,SUBSEP\)\; printf \"%s\\t%s\\t%s\n\", idx[1],idx[2],gensub\(\"^[,]+\",\"\",\"g\",snps[i]\)} }\' subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt \> subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}_temp.txt
    awk -F"\t" 'BEGIN{OFS="\t"} NR>=1{ for (i=2;i<=NF;i++) {snps[i,$i]=sprintf("%s,%0.'${char}'d",snps[i,$i],NR)} } END{for (i in snps) {split(i,idx,SUBSEP); printf "%s\t%s\t%s\n", idx[1],idx[2],gensub("^[,]+","","g",snps[i])} }' subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt > subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}_temp.txt # creates an array indexed by column nuumber and pattern type of the row numbers with the same pattern type
    printf "\n"

    echo The SNPs \(line nos, col 3\) and their pattern \(col 2\) in each chromosome tuple \(col 1\):
    echo head subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}_temp.txt
    head subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}_temp.txt
    printf "\n"

    echo Write snps in the same state in all chromosomes of the tuple:
    echo awk \-F\"\\t\" \'{a=\$2\; b=\$2\; if \(gsub\(\"[^1]\",\"\",a\)==0 \|\| gsub\(\"[^0]\",\"\",b\)==0\) {print \$0} }\' subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}_temp.txt \| sort \-k1,1 \| awk \'BEGIN{OFS=\"\\t\"} {tuple_pattern[\$1]=tuple_pattern[\$1]\(length\(tuple_pattern[\$1]\)==0?\"\":\",\"\)\$2\; tuple_snps[\$1]=sprintf\(\"%s%s%s\"\,tuple_snps[\$1]\,\(length\(tuple_snps[\$1]\)==0?\"\":\"\,\"\),gensub\(/\([,]\|\$\)/\,sprintf\(\"_%s%s\",substr\(\$2,1,1\),\"\\\\1\"\)\,\"g\"\,\$3\)\)}\; END{for \(i in tuple_pattern\) {split\(tuple_snps[i],snps,\",\"\)\; n=asort\(snps,dest\)\; str=dest[1]\; for \(j=2\;j\<=n\;j++\) {str=sprintf\(\"%s\,%s\"\,str\,dest[j]\)}\; print i,tuple_pattern[i],str} }\' \> snps_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
    awk -F"\t" '{a=$2; b=$2; if (gsub("[^1]","",a)==0 || gsub("[^0]","",b)==0) {print $0} }' subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}_temp.txt | sort -k1,1 | awk 'BEGIN{OFS="\t"} {tuple_pattern[$1]=tuple_pattern[$1](length(tuple_pattern[$1])==0?"":",")$2; tuple_snps[$1]=sprintf("%s%s%s",tuple_snps[$1],(length(tuple_snps[$1])==0?"":","),gensub(/([,]|$)/,sprintf("_%s%s",substr($2,1,1),"\\1"),"g",$3))}; END{for (i in tuple_pattern) {split(tuple_snps[i],snps,","); n=asort(snps,dest); str=dest[1]; for (j=2;j<=n;j++) {str=sprintf("%s,%s",str,dest[j])}; print i,tuple_pattern[i],str} }' > snps_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt # check if 0 characters are replaced when non-1's (or non-0's) are deleted; record SNPs (line numbers) where event occurs for each tuple of chromosomes; gensub returns 0 when no substitutions are made, so only print lines in which no 1's or no 0's are replaced in the SNP pattern, i.e., those whch are all-1 or all-0
    printf "\n"

    echo Chromosome tuples \(col 1\) with sharing pattern\(s\) \(col 2\) among the indicated SNPs \(col 3\):
    echo head snps_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
    head snps_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
    printf "\n"
    if [ -f subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}_temp.txt ];
    then
      echo rm subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}_temp.txt
      rm subjects_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}_temp.txt
      printf "\n"
    fi
  done

  # Look for files with the same partitioning in order to find snps specific to each bar grouping
  f_count=0;
  for file in snps_${file_ID}_${INDEX_NEW}_${COMBO_file}*.txt
  do
    echo $file
    let f_count=$f_count+1
  done
  let f_count=$f_count-1
  n_bars=$(echo $COMBO_NEW | awk "{print gsub(\"[|]\",\"\",\$0)}")
  echo Number of bars in the pattern is $n_bars and file count is $f_count beyond the first
  printf "\n"
  #Assemble awk code

  if (( $f_count > 0 ));
  then
    # first block: look for snps in file 1 but not file 2; store in array "found", indexed by line name; update line names to those also in file2
    str=$(echo "awk -F\"\t\" 'BEGIN{OFS=FS} FILENAME==ARGV[1]{line[\$1]=\$1; pattern[\$1]=\$2; array1[\$1]=\$3; next}; FILENAME==ARGV[2]{ if (\$1 in line) {split(array1[\$1],snps1,\",\"); split(\$3,snps2,\",\"); idx = 0; for (i1 in snps1) {count = 0; for (i2 in snps2) {if (snps1[i1]==snps2[i2]) { count++; break }} if (count == 0) {if (idx == 0) {found[\$1]=snps1[i1]; idx++} else {found[\$1]=found[\$1]\",\"snps1[i1];}} if (idx == 0) {found[\$1]=\" \"} }} next}" )
    for i in `seq 3 $f_count`
    do
      # middle blocks: update "found" array with snps not in files 3 to n-1; update line names to those also in file3, etc
      str=${str}$(echo " FILENAME==ARGV[$i]{ if (\$1 in line) {line[\$1]=\$1; split(found,snps_prev,\",\"); split(\$3,snps$i,\",\"); idx = 0; for (i in snps_prev) {count = 0; for (i$i in snps$i) {if (snps_prev[i] == snps$i[i$i]) { count++ }} if (count == 0) {if (idx == 0) {found[\$1]=snps_prev[i]; idx++} else {found[\$1]=found[\$1]\",\"snps_prev[i];}} if (idx == 0) {found[\$1]=\" \"} }} next}")
    done
    # last block: print snps in "found" array not in file n; exception is for n=2 when there is no need for "found" array
    let f_count=$f_count+1
    if (( $f_count == 2 ));
    then
      str=$(echo "awk -F\"\t\" 'BEGIN{OFS=FS} FILENAME==ARGV[1]{line[\$1]=\$1; pattern[\$1]=\$2; array1[\$1]=\$3; next}; (\$1 in line) {n=split(array1[\$1],snps1,\",\"); m=split(\$3,snps2,\",\"); asort(snps1); asort(snps2); idx = 0; for (i1 =1; i1<=n; i1++) {count = 0; for (i2=1; i2<=m; i2++) {if (snps1[i1]==snps2[i2]) { count++; break }} if (count == 0) {if (idx == 0) {\$2=pattern[\$1]; \$3=snps1[i1]; idx++} else {\$2=pattern[\$1]; \$3=\$3\",\"snps1[i1];}} if (idx == 0) {\$2=pattern[\$1]; \$3=\"\";} } print \$0}'" )
    else
      #last step prints values in file1 not in any of the others, i.e., checks if the snps previously "found" to be unique ("snps_prev") are still unique in the last file; also prints line number and pattern type
      str=${str}$(echo " FILENAME==ARGV[$f_count]{ if (\$1 in line) {n=split(found[\$1],snps_prev,\",\"); m=split(\$3,snps$f_count,\",\"); asort(snps_prev); asort(snps$f_count); idx = 0; for (i=1; i<=n; i++) {count = 0; for (i$f_count=1; i$f_count<=m; i$f_count++) {if (snps_prev[i] == snps$f_count[i$f_count]) { count++; break }} if (count == 0) {if (idx == 0) {\$2=pattern[\$1]; \$3=snps_prev[i]; idx++} else {\$2=pattern[\$1]; \$3=\$3\",\"snps_prev[i];}} if (idx == 0) {\$2=pattern[\$1]; \$3=\"\";} } print \$0 } }'")
    fi

    #loop over all ways of selecting a bar group to be file1
    echo Find unique SNPs in each file:
    for partition in $BAR_COMBO
    do
      partition_new=$(echo $partition | perl -pne "s/$max/j/g")
      partition_file=$(echo $partition_new | perl -pne "s/[|]/+/g")
      str1=${str}$(echo " snps_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt")
      for file in snps_${file_ID}_${INDEX_NEW}_${COMBO_file}*.txt
      do
        if [ ! $file == "snps_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt" ];
        then
          str1=${str1}$(echo " $file")
        fi
      done
      str1=${str1}$(echo " > snps_unique_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt")
      code3=$(echo $str1 | awk '{printf "%s",$0}') && printf "\n"

      echo Write unique snps in snps_unique_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
      echo $code3
      eval "$code3"
      printf "\n"

      # remove lines with fewer than the req'd number of chromosomes grouped, e.g., remove "1" lines when the rest are "11"
      n_comma=$(echo $partition_file | tr -cd "," | wc -c)
      let n_comma=$n_comma+1
      code4=$(echo "awk -F\"\t\" 'BEGIN{OFS=FS} { pattern=gensub(\"([^,]*).*\",\"\\\\1\",\"g\",\$2); if (length(pattern) == $n_comma) print \$0 }' snps_unique_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt > snps_unique_temp_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt") # get length of second field (pattern type) before the first comma, e.g., pattern 00,11 (indicating some shared SNPs are 0 and some SNPs are 1 in a chromosome double) returns length 2, but 0,1 returns length 1
      echo Remove lines with fewer than $n_comma chromosome$( (($n_comma!=1)) && echo "s" || echo "" ) in pattern:
      echo $code4
      eval "$code4"
      printf "\n"

      if [ -f snps_unique_temp_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt ]
      then
        echo cat snps_unique_temp_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt \> snps_unique_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
        echo rm snps_unique_temp_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
        cat snps_unique_temp_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt > snps_unique_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
        rm snps_unique_temp_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt

        echo head snps_unique_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
        head snps_unique_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt
        printf "\n"
      fi
    done

    # for exclusive patterns, make sure the snp appears in each file (with a different allele)
    echo Find alleles that are different in each file:
    str=$(echo "awk 'FILENAME==ARGV[1]{n=split(\$3,array,\",\"); for (i=1;i<=n;i++) {b=gensub(\"([0-9]*)_[0-9]*\",\"\\\\1\",\"g\",array[i]); snps1[NR][b]} next}")
    for i in `seq 2 $((f_count-1))`
    do
      str=${str}$(echo " FILENAME==ARGV[$i]{if (length(snps$((i-1))[FNR])>0) {n=split(\$3,array,\",\"); for (i=1;i<=n;i++) {b=gensub(\"([0-9]*)_[0-9]*\",\"\\\\1\",\"g\",array[i]); if (b in snps$((i-1))[FNR]) {snps$i[FNR][b]} } } next }")
    done
    str=${str}$(echo " FILENAME==ARGV[$f_count]{if (length(snps$((f_count-1))[FNR])>0) {str=\"\"; found=0; n=split(\$3,array,\",\"); for (i=1;i<=n;i++) {b=gensub(\"([0-9]*)_[0-9]*\",\"\\\\1\",\"g\",array[i]); if (b in snps$((f_count-1))[FNR]) {str=sprintf(\"%s%s%s\",str,(found>0?\",\":\"\"),array[i]); found++} } if (found>0) {\$3=str; print \$0} } next }'") && printf "\n"

    # loop over all segments of the partition of the COMBO and print alleles of snps that appear in all segments
    for partition in $BAR_COMBO
    do
      partition_new=$(echo $partition | perl -pne "s/$max/j/g")
      partition_file=$(echo $partition_new | perl -pne "s/[|]/+/g")
      str1=$str
      for file in snps_unique_${file_ID}_${INDEX_NEW}_${COMBO_file}*.txt
      do
        if [ ! $file == "snps_unique_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt" ];
        then
          str1=${str1}$(echo " $file")
        fi
      done
      str1=${str1}$(echo " snps_unique_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.txt > snps_unique_${file_ID}_${INDEX_NEW}_${COMBO_file}_${partition_file}.tmp")
      code5=$(echo $str1)
      echo $code5
      eval "$code5"
      printf "\n"
    done

    for file in snps_unique_${file_ID}_${INDEX_NEW}_${COMBO_file}*.tmp
    do
      test -f $file && echo mv $file ${file%.*}.txt
      test -f $file && mv $file ${file%.*}.txt && printf "\n"
    done
  else
    echo Write unique snps in ${file//snps/snps_unique}
    cp $file ${file//snps/snps_unique}
    printf "\n"
  fi
  for file in snps_${file_ID}_${INDEX_NEW}_${COMBO_file}*.txt
  do
    if [ -f $file ];
    then
      echo rm $file
      rm $file
    fi
  done
  printf "\n"
done
for file in stirling_${file_ID}_${INDEX}*
do
  if [ -f $file ];
  then
    echo rm $file
    rm $file
  fi
done

# Delete subjects files when done

for file in subjects_${file_ID}_${INDEX_NEW}*.txt
do
  if [ -f $file ]; then
    echo rm $file
    rm $file
  fi
done
for file in subjects_${file_ID}_${INDEX}*.txt
do
  if [ -f $file ];
  then
    echo rm $file
    rm $file
  fi
done

for file in combo_${file_ID}_${INDEX}*.txt
do
  if [ -f $file ];
  then
    echo rm $file
    rm $file
  fi
done
