# Chromosome Overlap #

This document outlines how to perform Chromosome Overlap in an HPC environment (LSF).  The goal is to obtain a list of closed haplotype patterns from (filtered meta-)chromosomes from a population consisting of cases and controls.

### Initiation ###

Your input should be a file <kbd>haplotype_estimates.cases.${NAME}.txt</kbd> of the following form:

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
 
Each row is chromosome from subject sid in the cases population; another file should exist for the controls.  Alternate alleles of SNP <kbd>rsid_Alt</kbd> are indicated with a 1 and reference alleles with a 0.

Find all unique chromosomes in <kbd>haplotype_estimates.cases.${NAME}.txt</kbd> using

```
awk 'NR==1{char=length(NF-1)} NR>1{printf "1\t"; for (i=2;i<=NF;i++) {printf "%s%0"char"d_%s",(i>2?",":""),i-1,$i}; printf "\n"}' haplotype_estimates.cases.${NAME}.txt | awk 'BEGIN{OFS="\t"} {seen[$2]+=$1} END{for (i in seen) {print seen[i],i} }' > Pattern_combined.Iteration000.${NAME}_${PATTERN}.txt
```

where <kbd>${NAME}</kbd> is a location label (what chromosome and position do your haplotypes come from?) and <kbd>${PATTERN}</kbd> is a suffix indicating the type of overlap (e.g., <kbd>2,j</kbd> indicates pairwise overlap of and <kbd>2,3,j</kbd> triple-wise).  These names will be used for combining the results of different jobs.  The output <kbd>Pattern_combined.Iteration000.${NAME}_${PATTERN}.txt</kbd> will be a two-column table of counts of meta-chromosomes:

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

Each meta-chromosome is a list of SNPs (column numbers in your haplotype file) and their corresponding alleles (0 or 1).

Perform the initial overlap with at least 50 overlaps per job by running

```
bsub -P ChromosomeOverlap -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.err -R "rusage[mem=512]" "sh ChromosomeOverlap_iteration_sub_parallel.v3.sh \"${NAME}\" 2 \"${Pattern}\" 50 0"
```

The overlaps are partioned to nodes using <kbd>index2combo2.R</kbd>, which assigns a unique number to each possible overlap of <kbd>$SIGMA</kbd> chromosomes.  Here, <kbd>50</kbd> indicates that the total number of overlaps per job should be 50, which value doubles until the total number of jobs required to complete all pairwise overlaps is below a threshold <kbd>MAX_JOBS</kbd>; <kbd>2</kbd> indicates that the overlaps should consist of two chromosomes each; and <kbd>0</kbd> indicates that no overlaps have been completed yet (corresponding to Iteration000).  The output of this step is a file <kbd>Pattern_combined.Iteration001.${NAME}\_${PATTERN}.txt</kbd>, containing counts of meta-chromosomes and their corresponding SNP-alleles, and another file <kbd>Closed_patterns_stats.${NAME}_${PATTERN}.txt</kbd>, which keeps track of all the closed patterns found so far and the iteration at which they were discovered.  The program automatically halts at iteration 0 to allow the user to apply p-value-based filtering.

### Filtering ###

Get the segregation of each pattern discovered at iteration 1 in <kbd>Closed_patterns_stats.${NAME}_${PATTERN}.txt</kbd> between cases and controls using 

```
bsub -P ChromosomeOverlap -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.err -R "rusage[mem=256]" "sh ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.cases.{NAME}.txt,haplotype_estimates.controls.${NAME}.txt Closed_patterns_stats.${NAME}_${PATTERN}.txt 50 \"1\" \"Iteration001.${NAME}\""
```

where again <kbd>50</kbd> is the initial number of overlaps assigned to each job and <kbd>1</kbd> instructs the programto only take patterns in <kbd>Closed_patterns_stats.${NAME}_${PATTERN}.txt</kbd> discovered at iteration 1.  The counts are then fed to <kbd>ChromosomeOverlap_fisher_exact.R</kbd> for computation of the p-values.  The result is a file <kbd>fisher_exact.Iteration001.${NAME}.patterns_0000001-9999999.txt</kbd> with OR in the fourth field and p-value in the fifth.  Filtering can then be applied to update the <kbd>Pattern_combined</kbd> file using

```
awk 'NR==FNR{seen[$1]; next} ($2 in seen){print $0}' <(awk '(NR>1 && $5+0<1e-9 && $4+0>1){print $1}' fisher_exact.Iteration001.${NAME}.patterns_0000001-9999999.txt) Pattern_combined.Iteration001.chr11.69231642-69431642_2,j.txt > Pattern_combined.Iteration001.${NAME}_${PATTERN}.tmp
mv Pattern_combined.Iteration001.${NAME}_${PATTERN}.txt Pattern_combined_old.Iteration001.${NAME}_${PATTERN}.txt
mv Pattern_combined.Iteration001.${NAME}_${PATTERN}.tmp Pattern_combined.Iteration001.${NAME}_${PATTERN}.txt
```

which selects risk-increasing patterns with p-values less than <kbd>1e-9</kbd>.  The unfiltered list is saved as <kbd>Pattern_combined_old.Iteration001.${NAME}_${PATTERN}.txt</kbd> in case the user should want to try a different threshold, in which case he should first run

```
Pattern_combined_old.${NAME}_${PATTERN}.txt Pattern_combined.${NAME}_${PATTERN}.txt
```

### Iteration ###

Once the filtered patterns are selected, run

```
bsub -P ChromosomeOverlap -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.${NAME} -oo ChromosomeOverlap_iteration_sub_parallel.v3.${NAME}.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.${NAME}.err -R "rusage[mem=512]" "sh ChromosomeOverlap_iteration_sub_parallel.v3.sh \"${NAME}\" 2 \"${PATTERN}\" 50 1"
```

to initiate the overlaps starting from iteration 1.  All new patterns discovered at each iteration will be added to <kbd>Closed_patterns_stats.${NAME}_${PATTERN}.txt</kbd>.  The program ceases when no new patterns are discovered, typically after five iterations have been completed.  The patterns can be assessed for association with disease using Fisher's exact test

```
bsub -P ChromosomeOverlap -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001=005.err -R "rusage[mem=256]" "sh ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.cases.{NAME}.txt,haplotype_estimates.controls.${NAME}.txt Closed_patterns_stats.${NAME}_${PATTERN}.txt 50 \"\" \"Iteration001-005.${NAME}\""
```

This code finds the counts of all patterns in <kbd>Closed_patterns_stats.${NAME}_${PATTERN}.txt</kbd>, starting from iteration 0.  But if only new patterns are desired, one can create a subset file

```
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1){print $0}' Pattern_combined.Iteration001.${NAME}_${PATTERN}.txt Closed_patterns_stats.${NAME}_${PATTERN}.txt > Closed_patterns_stats.Iteration001-005.${NAME}_${PATTERN}.txt

```

on which <kbd>ChromosomeOverlap_iteration_sub_parallel.v3.sh</kbd> can be run.
