# Chromosome Overlap scripts #

* Use <kbd>ChromosomeOverlap_iteration_sub_parallele.v3.sh</kbd> on a list of unique contiguous haplotypes (Iteration000 chromosomes) to form all pairwise overlaps
* Use <kbd>ChromosomeOverlap_haplotype_count_sub.v5.sh</kbd> to get the counts of each meta-chromosome from Iteration001 in cases and controls and determine the Fisher's exact test p-values
* Use <kbd>ChromosomeOverlap_iteration_sub_parallele.v3.sh</kbd> again starting from Iteration001 to complete the overlaps on a filtered set of meta-chromosomes
* Use <kbd>index2combo2.R</kbd> to find the unique multiindex i_1,i_2,...,i_sigma corresponding to the I<sup>th</sup> sigma-combination of n objects in a sigma-dimensional triangular array