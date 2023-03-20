# Chromosome Overlap #

This document outlines how to perform Chromosome Overlap in an HPC environment (LSF).  The goal is to obtain a list of closed haplotype patterns from (filtered meta-)chromosomes from a population.

### Initiation ###

Your input should be a <kbd>.txt</kbd> file of the following form:

<table>
  <tr>
    <th>sid</th>
    <th>rs0000001_A</th>
    <th>rs0000002_G</th>
    <th>rs0000003_C</th>
    <th>...</th>
  </tr>
  <tr>
    <th>0001</th>
    <td>0</td>
    <td>0</td>
    <td>0</td>
    <td>...</td>
  </tr>
  <tr>
    <th>0001</th>
    <td>0</td>
    <td>1</td>
    <td>0</td>
    <td>...</td>
  </tr>
  <tr>
    <th>0002</th>
    <td>1</td>
    <td>0</td>
    <td>0</td>
    <td>...</td>
  </tr>
  <tr>
    <th>0002</th>
    <td>1</td>
    <td>0</td>
    <td>0</td>
    <td>...</td>
  </tr>
 </table>
 
Each row is chromosome from subject sid; note that each has two chromosomes with the same id.  Alternate alleles of SNP </kbd>rsid_Alt</kbd> are indicated with a 1 and reference alleles with a 0.

Find all unique chromosomes in <kbd>file.txt</kbd> using

```
awk 'NR==1{char=length(NF-1)} NR>1{printf "1\t"; for (i=2;i<=NF;i++) {printf "%s%0"char"d_%s",(i>2?",":""),i-1,$i}; printf "\n"}' file.txt | awk 'BEGIN{OFS="\t"} {seen[$2]+=$1} END{for (i in seen) {print seen[i],i} }' > Pattern_combined.Iteration000.${NAME}_${PATTERN}.txt
```

where <kbd>${NAME}</kbd> is a location label (what chromosome and position do your haplotypes come from?) and <kbd>${PATTERN}</kbd> is a suffix indicating the type of overlap (e.g., <kbd>2,j</kbd> indicates pairwise overlap of and <kbd>2,3,j</kbd> triple-wise).  These names will be used for combining the results of different jobs.  The output <kbd>Pattern_combined.Iteration000.${NAME}_${PATTERN}.txt</kbd> will be a two-column list of counts of meta-chromosomes:

<table>
<tr>
  <td>2</td>
  <td>01_1,02_0,04_0,...</td>
</tr>
<tr>
  <td>1</td>
  <td>01_1,02_1,03_0,...</td>
</tr>
</table>

Each meta-chromosome is a list of SNPs (the coumn number in your haplotype file) and their corresponding alleles (0 or 1).

Perform the initial overlap with at least 50 overlaps per job by running

```
bsub -P ChromosomeOverlap -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.err -R "rusage[mem=512]" "sh ChromosomeOverlap_iteration_sub_parallel.v3.sh \"${NAME}\" 2 \"${Pattern}\" 50 0"
```

The overlaps are partioned to nodes using <kbd>index2combo2.R</kbd>, which assigns a unique number to each possible overlap of <kbd>$SIGMA</kbd> chromosomes.  Here, <kbd>50</kbd> indicates that the total number of overlaps per job should be 50, which value doubles until the total number of jobs required to complete all pairwise overlaps is below a threshold <kbd>MAX_JOBS</kbd>; <kbd>2</kbd> indicates that the overlaps should consist of two chromosomes each; and <kbd>0</kbd> indicates that no overlaps have been completed yet (corresponding to Iteration000).  The output of this step is a file <kbd>Pattern_combined.Iteration001.${NAME}_${PATTERN}.txt</kbd> containing counts of meta-chromosomes and their corresponding SNP-alleles.  The program automatically halts at iteration 0 to allow the user to apply p-value-based filtering.

### Filtering ###
