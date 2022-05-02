# ChromosomeOverlap

Chromosome Overlap is a haplotype pattern-mining algorithm for use on phased genotype data.  Chromosome Overlap runs in an HPC environment on Platform LSF.  The program has the following features.

## Triangular-array algorithm

For a given number n of items, there are n-choose-k k-tuples of items.  What k-tuple is the i<sup>th</sup>?  The Rscript index2combo2.R finds the answer.

## Initial round of overlaps

Form all σ-tuples of chromosomes from affected individuals and look for shared allele patterns.  Filter by p-value in comparison with a control group.

## Iterated overlap

Overlap pairs of meta-chromosomes to find *closed* patterns.  Look for patterns that appear for the last time in any given iteration.

## Model-fitting

Test the closed patterns for association with the outcome.
