#! /bin/bash
set -e

# finds the multiindex i_1,i_2,...,i_SIGMA corresponding to a linear index L_INDEX in the triangular array of N-choose-SIGMA combinations without repeats (i_1 < i_2 < ... < i_SIGMA)

# sh /home/wletsou/scripts/index2combo2.sh

L_INDEX=$1 # linear index ranging from 0 to (N-choose-SIGMA - 1)
N=$2 # total number of objects
SIGMA=$3 # size of the combination
REPEATS=$4 # nonempty if i_1 <= i_2 <= ... <= i_SIGMA

[ -z $L_INDEX ] && (>&2 echo "Index not supplied"; exit 1)
[ -z $N ] && (>&2 echo "Total not supplied"; exit 1)
[ -z $SIGMA ] && (>&2 echo "Combo size not supplied"; exit 1)
(( $SIGMA > $N )) && (>&2 echo "Invalid combo size"; exit 1)

L_ARRAY+=($L_INDEX) # initiate array of linear indices in a triangular array of SIGMA-i+1 dimensions

max=$(awk 'BEGIN{print binom('$N','$SIGMA')} function binom(n,x) {num=1; den=1; for (i=1;i<=x;i++) {num=num*(n-i+1); den=den*i}; return(num/den)}')
(( $L_INDEX>=$max )) && (>&2 echo "Index too large"; exit 1)

if [ ! -z $REPEATS ] # like evaluating i_1 < i_2 < ... < i_SIGMA when there are SIGMA-1 more objects (just enough so i_2,...,i_SIGMA can all be equal to i_1)
then
  N=$((N+SIGMA-1))
fi

multiindex+=(0) # "i_0"

for ((i=1;i<=$SIGMA;i++))
do
  k=0
  epsilon=0
  while (($((epsilon+0))<=0))
  do
    epsilon=$(awk 'BEGIN{print binom('$N'-'$((multiindex[$((i-1))]+0))','$SIGMA'-'$i'+1)-binom('$N'-'$((multiindex[$((i-1))]+0))'-'$k','$SIGMA'-'$i'+1)-'${L_ARRAY[$((i-1))]}'} function binom(n,x) {num=1; den=1; for (i=1;i<=x;i++) {num=num*(n-i+1); den=den*i}; return(num/den)}') # number of indices including L[i-1] and beyond in the triangular array of SIGMA-i dimensions
    k=$((k+1))
  done
  L_NEW=$(awk 'BEGIN{print binom('$N'-'$((multiindex[$((i-1))]+0))'-'$((k-1))','$SIGMA'-'$i')-'$epsilon'} function binom(n,x) {num=1; den=1; for (i=1;i<=x;i++) {num=num*(n-i+1); den=den*i}; return(num/den)}') # Number of indices in the triangular array of SIGMA-i dimensions less the error beyond (and including) L[i-1] is the new index L[i]
  L_ARRAY+=($L_NEW)
  multiindex+=($((multiindex[$((i-1))]+k-1))) # total number of layers that had to be removed
done

if [ ! -z $REPEATS ] # i_k starts at k; reduce i_k by k-1 so that the first multiindex is 1,1,...,1
then
  for ((i=1;i<=$((SIGMA+1));i++)) # includes i_0
  do
    multiindex[$((i-1))]=$((multiindex[$((i-1))]-(i-1)+1))
  done
fi

echo ${multiindex[@]:1:${#multiindex[@]}} # multiindex i_1,i_2,...i_SIGMA corresponding to linear index L_INDEX
