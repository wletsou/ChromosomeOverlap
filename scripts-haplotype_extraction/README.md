# Haplotype extraction #

* Use <kbd>ukbb_hybrid_haplotyp3.sh</kbd> to extract genotyped SNPs from the phased UKBB data
* Use <kbd>ukbb_topmed_merge.sh</kbd> to create a TOPMed-phased UKBB vcf file
* Use <kbd>dbgap_haplotype_merge.sh</kbd> to create a TOPMed-phased DRIVE vcf file
* <kbd>ukbb_submission.v2.sh</kbd> is a wrapper script that extracts haplotypes from the UKBB and DRIVE vcf files at a common set of SNPs
* Use <kbd>ukbb_haplotype_extract3.sh</kbd> to create haplotype_estimates files for cases and controls
* Use <kbd>ukbb_haplotype_permute.v2</kbd> to randomly permute cases and controls such that the number of hetero- and homozygotes is preserved in each group

For theses scripts, the UKBB phenotype files have subject id in the first column and case/control status in the third column, and the DRIVE phenotype files have subject id in the first column and case/control status in the second column.