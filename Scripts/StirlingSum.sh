#! /bin/bash
#Compute sum over k of the Stirling numbers
# S(n,k) = (1 / k!) * sum_{i=0}^{k} (-1)^i * (k choose i) * (k - i)^n
# of the second kind http://mathworld.wolfram.com/StirlingNumberoftheSecondKind.html
#Gives the total number of partitions of n objects such that no subset is empty
n=$1

S=0
i=0
Ski=0;
for ((k=1; k<$n; k++))
do
   #k^n https://stackoverflow.com/questions/13111967/shell-scripting-raise-to-the-power
  k_factorial=$k;
  for ((j=$((k-1)); j>0; j--))
  do
    k_factorial=$((k_factorial*j))
  done
  #echo k! = $k_factorial
  let i=0
  Ski=$((k**n));
  Sk=$Ski
  #echo term $i is $Ski
  for ((i=1; i<$k; i++))
  do
    num0=$k;
    num=1;
    den=1;
    for ((j=1; j<=$i; j++))
    do
      let num=$num*$((num0-j+1))
      let den=$den*$j
    done
    k_choose_i=$((num/den))
    let x=$((k_choose_i*$(($((k-i))**n))*$(($((0-1))**i)))) #(k-1)^n https://stackoverflow.com/questions/13111967/shell-scripting-raise-to-the-power
    let Ski=$((Ski+k_choose_i*$(($((k-i))**n))*$(($((0-1))**i))))
    #echo term $i is $x, subtotal $Ski
  done

  Sk=$((Ski/k_factorial))
  #echo S\($n,$k\) = $Sk
  S=$((S+Sk))
done
if (($k==$n)) && (($((i+1))>1))
then
  Sk=1 #for 0^n in the last term
  #echo term 1 is 1
  #echo S\($n,$n\) = $Sk
  S=$((S+Sk))
fi
echo $S
