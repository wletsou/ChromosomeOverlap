mkdir UKBB_chr22.6
cd UKBB_chr22.6

mkdir /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_genotyped_snps
cd ukbb_genotyped_snps

awk 'BEGIN{OFS="\t"} ($3==1){print $1,$1}' /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_cases.indiv
awk 'BEGIN{OFS="\t"} ($3==0){print $1,$1}' /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_controls.indiv

bsub -P SJLIFE -J ukbb_hybrid_haplotype3.chr22.28000000-30000000 -oo ukbb_hybrid_haplotype3.chr22.28000000-30000000.out -eo ukbb_hybrid_haplotype3.chr22.28000000-30000000.err -R "rusage[mem=20000]" "sh /home/wletsou/scripts/ukbb_hybrid_haplotype3.sh ukbb_bca_cases.indiv,ukbb_bca_controls.indiv 22 28000000,30000000 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/biobank/gds/ukb.bca.hap/ukb.bca.hap.chr22.new.vcf.gz /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/biobank/gds/ukb.bca.hap/ukb.bca.hap.qc.anno.txt"

cd /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6

awk 'BEGIN{OFS="\t"} ($3==1){print $1,$1}' /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_cases.indiv
awk 'BEGIN{OFS="\t"} ($3==0){print $1,$1}' /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_controls.indiv

# get genotyped SNPs in region
awk 'BEGIN{OFS="\t"} ($2>=28145327 && $2<=28270372){print $3}' /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_genotyped_snps/ukbb_bca_cases+ukbb_bca_controls.chr22.28000000-30000000.haplotypes.vcf > ukbb_snp_list.chr22.28145327-28270372.txt

# check that SNPs have good "imputation" quality
awk 'NR==FNR{snp[$1]; next} ($8 in snp && $15==1){print $8}' ukbb_snp_list.chr22.28145327-28270372.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr22.qced.anno.info > ukbb_snp_list.chr22.28145327-28270372.tmp
test -f ukbb_snp_list.chr22.28145327-28270372.tmp && mv ukbb_snp_list.chr22.28145327-28270372.tmp ukbb_snp_list.chr22.28145327-28270372.txt

# vcf file for UKBB
bsub -P SJLIFE -J ukbb_topmed_merge.chr22.28145327-28270372 -oo ukbb_topmed_merge.chr22.28145327-28270372.out -eo ukbb_topmed_merge.chr22.28145327-28270372.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_topmed_merge.sh ukbb_snp_list.chr22.28145327-28270372.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr22.qced.anno.info 22 28145327,28270372 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/raw/chr22"

# vcf file for dbGaP
bsub -P SJLIFE -J dbgap_haplotype_merge.chr22.28145327-28270372 -oo dbgap_haplotype_merge.chr22.28145327-28270372.out -eo dbgap_haplotype_merge.chr22.28145327-28270372.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/dbgap_haplotype_merge.sh dbgap28544_cases,dbgap28544_controls /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/chr22.qced.info ukbb_snp_list.chr22.28145327-28270372.txt 22 28145327,28270372 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/vcfs/dbgap28544.topmedr2.cleaned.hg38.chr22.vcf.gz"

# extract haplotypes for cases and controls
bsub -P SJLIFE -J ukbb_bca_overlap.chr22.28145327-28270372 -oo ukbb_bca_overlap.chr22.28145327-28270372.out -eo ukbb_bca_overlap.chr22.28145327-28270372.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_submission.v2.sh ukbb_bca_cases,ukbb_bca_controls dbgap28544_cases,dbgap28544_controls \"\" 22 28145327,28270372 \"\" \"\" ukbb_snp_list.chr22.28145327-28270372.txt ukbb.topmed.hg19_chr22.28145327-28270372.hg38.vcf.gz,dbgap28544_cases+dbgap28544_controls.hg19_chr22.28145327-28270372.hg38.vcf.gz /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr22.qced.anno.info \"\" \"\" \"\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6 /home/wletsou/scripts"

awk 'NR==1{char=length(NF-1)} NR>1{printf "1\t"; for (i=2;i<=NF;i++) {printf "%s%0"char"d_%s",(i>2?",":""),i-1,$i}; printf "\n"}' haplotype_estimates.ukbb_bca_cases.chr22.28145327-28270372.txt | awk 'BEGIN{OFS="\t"} {seen[$2]+=$1} END{for (i in seen) {print seen[i],i} }' > Pattern_combined.Iteration000.chr22.28145327-28270372_2,j.txt # starting with all chromosomes, not just pairs

# First round of overlaps:
bsub -P SJLIFE -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr22.28145327-28270372 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr22.28145327-28270372.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr22.28145327-28270372.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ChromosomeOverlap_iteration_sub_parallel.v3.sh chr22.28145327-28270372 2 2,j 50 0 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6 /home/wletsou/scripts"

# p-values for first round of overlaps
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr22.28145327-28270372.txt,haplotype_estimates.ukbb_bca_controls.chr22.28145327-28270372.txt Closed_patterns_stats.chr22.28145327-28270372_2,j.txt 50 \"1\" \"Iteration001.chr22.28145327-28270372\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6 /home/wletsou/scripts"

# remove older results if staring again
rm -r *Iteration00[2-9]*
awk '($1<2){print $0}' Closed_patterns_stats.chr22.28145327-28270372_2,j.txt > Closed_patterns_stats.chr22.28145327-28270372_2,j.tmp
mv Closed_patterns_stats.chr22.28145327-28270372_2,j.tmp Closed_patterns_stats.chr22.28145327-28270372_2,j.txt

test -f Pattern_combined_old.Iteration001.chr22.28145327-28270372_2,j.txt && mv Pattern_combined_old.Iteration001.chr22.28145327-28270372_2,j.txt Pattern_combined.Iteration001.chr22.28145327-28270372_2,j.txt
awk 'NR==FNR{seen[$1]; next} ($2 in seen){print $0}' <(awk '(NR>1 && $5+0<2e-12 && $4+0>1){print $1}' fisher_exact.Iteration001.chr22.28145327-28270372.patterns_0000001-7759059.txt) Pattern_combined.Iteration001.chr22.28145327-28270372_2,j.txt > Pattern_combined.Iteration001.chr22.28145327-28270372_2,j.tmp
# awk 'NR==FNR{seen[$1]; next} ($2 in seen){print $0}' <(awk '(NR>1 && $5+0<5e-13 && $4+0>1){print $0}' fisher_exact.Iteration001.chr22.28145327-28270372.patterns_0000001-7759059.txt | awk 'BEGIN{array[0]} {delete array; n=split($1,array,","); if (n<snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]+0 || snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]+0==0) {snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]=n; seen[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]=$1} } END{for (i in seen) {print seen[i]} }') Pattern_combined.Iteration001.chr22.28145327-28270372_2,j.txt
mv Pattern_combined.Iteration001.chr22.28145327-28270372_2,j.txt Pattern_combined_old.Iteration001.chr22.28145327-28270372_2,j.txt
mv Pattern_combined.Iteration001.chr22.28145327-28270372_2,j.tmp Pattern_combined.Iteration001.chr22.28145327-28270372_2,j.txt

# Begin iterations:
bsub -P SJLIFE -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr22.28145327-28270372 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr22.28145327-28270372.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr22.28145327-28270372.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ChromosomeOverlap_iteration_sub_parallel.v3.sh chr22.28145327-28270372 2 2,j 50 1 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6 /home/wletsou/scripts"

# (filtered) patterns appearing for the first time at Iterations 001-005
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1){print $0}' Pattern_combined.Iteration001.chr22.28145327-28270372_2,j.txt Closed_patterns_stats.chr22.28145327-28270372_2,j.txt > Closed_patterns_stats.Iteration001-005.chr22.28145327-28270372_2,j.txt

# p-values for Iterations 001-005
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr22.28145327-28270372.txt,haplotype_estimates.ukbb_bca_controls.chr22.28145327-28270372.txt Closed_patterns_stats.Iteration001-005.chr22.28145327-28270372_2,j.txt 50 \"\" \"Iteration001-005.chr22.28145327-28270372\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6 /home/wletsou/scripts"

# top closed patterns
awk '($5+0<1e-9 && $4+0>1){print $0}' fisher_exact.Iteration001-005.chr22.28145327-28270372.patterns_000001-398251.txt > fisher_exact.Iteration001-005.chr22.28145327-28270372.patterns_000001-398251.Results.txt

# tranlate haplotype to rsid
bsub -P SJLIFE -J ukbb_haplotype_translate2_sub -oo ukbb_haplotype_translate2_sub.out -eo ukbb_haplotype_translate2_sub.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr22.28145327-28270372.txt fisher_exact.Iteration001-005.chr22.28145327-28270372.patterns_000001-398251.Results.txt 50"

# precompute haplotype counts
bsub -P SJLIFE -J ukbb_haplotype_model9_sub.ukbb_discovery -oo ukbb_haplotype_model9_sub.v2.ukbb_discovery.out -eo ukbb_haplotype_model9_sub.v2.ukbb_discovery.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt fisher_exact.Iteration001-005.chr22.28145327-28270372.patterns_000001-398251.Results.translated.txt \"\" \"\" \"\" 22 28145327,28270372 50 \"ukbb_discovery\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6 /home/wletsou/scripts"

bsub -P SJLIFE -J ukbb_haplotype_model9_sub.drive_replication -oo ukbb_haplotype_model9_sub.v2.drive_replication.out -eo ukbb_haplotype_model9_sub.v2.drive_replication.err -R "rusage[mem=256]" -w "done(178595692)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr22.28145327-28270372.txt fisher_exact.Iteration001-005.chr22.28145327-28270372.patterns_000001-398251.Results.translated.txt \"\" \"\" \"\" 22 28145327,28270372 50 \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6 /home/wletsou/scripts"

# forward selection
bsub -P SJLIFE -J ukbb_haplotype_model9_iterate.v2.discovery -oo ukbb_haplotype_model9_iterate.v2.discovery.out -eo ukbb_haplotype_model9_iterate.v2.discovery.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.v2.sh ukbb_discovery.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt \"\" \"\" 1e-5 ukbb_discovery.Significant_patterns.txt \"ukbb_discovery\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6 /home/wletsou/scripts"

bsub -P SJLIFE -J dbgap28544_haplotype_model9_iterate.v2.replication -oo dbgap28544_haplotype_model9_iterate.v2.replication.out -eo dbgap28544_haplotype_model9_iterate.v2.replication.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/dbgap28544_haplotype_model9_iterate.v2.sh drive_replication.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt \"\" \"\" 1e-5 drive_replication.Significant_patterns.txt \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6 /home/wletsou/scripts"

# replication
awk -v alpha=0.05 'BEGIN{OFS = "\t"} NR==FNR && $2 > 1 && $5 < 1e-5{seen[$1]=$0; n+=1; next} ($1 in seen && $2 > 1 && $5 < alpha){m+=1} END{printf("%s of %s replicated at p < %0.2e\n",m + 0,n + 0,alpha)}' ukbb_discovery.Conditional_haplotype_effects.h1.txt drive_replication.Conditional_haplotype_effects.h1.txt
# 29080 of 31857 replicated at p < 5.00e-02

# including contiguous haplotypes (Iteration000)

# (filtered) patterns appearing for the first time at Iterations 000,001-005
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1 || $1 == 0){print $0}' Pattern_combined.Iteration001.chr22.28145327-28270372_2,j.txt Closed_patterns_stats.chr22.28145327-28270372_2,j.txt > Closed_patterns_stats.Iteration000+Iteration001-005.chr22.28145327-28270372_2,j.txt

# p-values for Iterations 000,001-005
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr22.28145327-28270372.txt,haplotype_estimates.ukbb_bca_controls.chr22.28145327-28270372.txt Closed_patterns_stats.Iteration000+Iteration001-005.chr22.28145327-28270372_2,j.txt 50 \"\" \"Iteration000+Iteration001-005.chr22.28145327-28270372\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6 /home/wletsou/scripts"

# top closed patterns
awk '($5+0<1e-9 && $4+0>1){print $0}' fisher_exact.Iteration000+Iteration001-005.chr22.28145327-28270372.patterns_000001-402472.txt > fisher_exact.Iteration000+Iteration001-005.chr22.28145327-28270372.patterns_000001-402472.Results.txt

# tranlate haplotype to rsid
bsub -P SJLIFE -J ukbb_haplotype_translate2_sub -oo ukbb_haplotype_translate2_sub.Iteration000+Iteration001-005.out -eo ukbb_haplotype_translate2_sub.Iteration000+Iteration001-005.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr22.28145327-28270372.txt fisher_exact.Iteration000+Iteration001-005.chr22.28145327-28270372.patterns_000001-402472.Results.txt 50"

# precompute haplotype counts
bsub -P SJLIFE -J ukbb_haplotype_model9_sub.ukbb_discovery.Iteration000+Iteration001-005 -oo ukbb_haplotype_model9_sub.v2.ukbb_discovery.Iteration000+Iteration001-005.out -eo ukbb_haplotype_model9_sub.v2.ukbb_discovery.Iteration000+Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt fisher_exact.Iteration000+Iteration001-005.chr22.28145327-28270372.patterns_000001-402472.Results.translated.txt \"\" \"\" \"\" 22 28145327,28270372 50 \"ukbb_discovery.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6 /home/wletsou/scripts"

bsub -P SJLIFE -J ukbb_haplotype_model9_sub.drive_replication.Iteration000+Iteration001-005 -oo ukbb_haplotype_model9_sub.v2.drive_replication.Iteration000+Iteration001-005.out -eo ukbb_haplotype_model9_sub.v2.drive_replication.Iteration000+Iteration001-005.err -R "rusage[mem=256]" -w "done(190034103)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr22.28145327-28270372.txt fisher_exact.Iteration000+Iteration001-005.chr22.28145327-28270372.patterns_000001-402472.Results.translated.txt \"\" \"\" \"\" 22 28145327,28270372 50 \"drive_replication.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6 /home/wletsou/scripts"

# forward selection
bsub -P SJLIFE -J ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-005 -oo ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-005.out -eo ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-005.err -w "done(190037260)" -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.v2.sh ukbb_discovery.Iteration000+Iteration001-005.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt \"\" \"\" 1e-5 ukbb_discovery.Iteration000+Iteration001-005.Significant_patterns.txt \"ukbb_discovery.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6 /home/wletsou/scripts"

bsub -P SJLIFE -J dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-005 -oo dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-005.out -eo dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-005.err -w "done(190428534)" -R "rusage[mem=256]" "sh /home/wletsou/scripts/dbgap28544_haplotype_model9_iterate.v2.sh drive_replication.Iteration000+Iteration001-005.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt \"\" \"\" 1e-5 drive_replication.Iteration000+Iteration001-005.Significant_patterns.txt \"drive_replication.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6 /home/wletsou/scripts"

# replication
awk -v alpha=0.05 'BEGIN{OFS = "\t"} NR==FNR && $2 > 1 && $5 < 1e-5{seen[$1]=$0; n+=1; next} ($1 in seen && $2 > 1 && $5 < alpha){m+=1} END{printf("%s of %s replicated at p < %0.2e\n",m + 0,n + 0,alpha)}' ukbb_discovery.Iteration000+Iteration001-005.Conditional_haplotype_effects.h1.txt drive_replication.Iteration000+Iteration001-005.Conditional_haplotype_effects.h1.txt
# 29081 of 31858 replicated at p < 5.00e-02

# reduction
sh /home/wletsou/scripts/haplotype_rpart.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0" "" 1 # h1
#  reduced = rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0

awk 'BEGIN{OFS="\t"} ("rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0" ~ $3){print $1,$2,$3,$4,$5}' ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr22.28145327-28270372.haplotypes.hg38.vcf # h1 rows of vcf file

# chr22   27749339        rs41277849      G       A
# chr22   27753172        rs5762341       A       C
# chr22   27754121        rs4822935       A       G
# chr22   27755685        rs4820776       A       G
# chr22   27759416        rs2267113       C       T
# chr22   27760800        rs2014274       A       G
# chr22   27763546        rs74562526      T       G
# chr22   27767358        rs7292322       G       A
# chr22   27773533        rs1297593       C       A
# chr22   27774178        rs1807510       A       C
# chr22   27774229        rs2005397       C       T
# chr22   27775522        rs8142823       C       T
# chr22   27782526        rs45512704      C       T
# chr22   27783835        rs1297595       T       C
# chr22   27785448        rs382819        A       C
# chr22   27786522        rs424345        G       A
# chr22   27791408        rs45462093      A       G
# chr22   27791785        rs45574432      T       C
# chr22   27802203        rs45540534      T       G
# chr22   27808809        rs11704546      T       C
# chr22   27810924        rs11705555      A       C
# chr22   27874384        rs5997320       T       G

# discovery
sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0" "" "" "" "" h1.discovery
# LRT p value = 8.587212e-15 on 1 degree(s) of freedom with HR = 2.583979e+00 (2.099458e+00 to 3.180321e+00), frequency = 2.005148e-03 (726) = 4.993896e-03/1.848590e-03 (90/636), p = 4.056955e-15

# DRIVE replication
sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0" "" "" "" "" h1.replication
# LRT p value = 4.608411e-03 on 1 degree(s) of freedom with OR = 1.386348e+00 (1.103208e+00 to 1.742158e+00), frequency = 2.854768e-03 (316) = 3.259713e-03/2.373230e-03 (196/120), p = 6.553832e-03

# lead haplotype
awk 'BEGIN{haplotype="11111101111111111101111101111111111011111111101111"; n=split(haplotype,array,"")} NR==FNR{hap[$3]=array[FNR]; ref[$3]=$4; alt[$3]=$5; next} ($3 in hap){ printf "%s%s_%s=%s",(seen+0>0?",":""),$3,($5==alt[$3]?alt[$3]:ref[$3]),($5==alt[$3]?hap[$3]:1-hap[$3]); seen+=1 } END{printf "\n"}' <(awk '($2 >=28145327 && $2<=28270372 ){print $0}' /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_genotyped_snps/ukbb_bca_cases+ukbb_bca_controls.chr22.28000000-30000000.haplotypes.vcf) ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr22.28145327-28270372.haplotypes.hg38.vcf
# rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs118009109_T=0,rs4822935_G=1,rs2283844_A=1,rs4820776_G=1,rs2267106_A=0,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs8136505_T=0,rs74562526_G=0,rs16985531_G=0,rs7292322_A=0,rs17477193_G=0,rs2283846_G=0,rs5997301_G=1,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs4822939_G=0,rs5752634_G=0,rs2073780_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45443199_T=0,rs12166473_G=1,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs5752639_G=0,rs5752640_A=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0

awk 'BEGIN{OFS="\t"} ("rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs118009109_T=0,rs4822935_G=1,rs2283844_A=1,rs4820776_G=1,rs2267106_A=0,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs8136505_T=0,rs74562526_G=0,rs16985531_G=0,rs7292322_A=0,rs17477193_G=0,rs2283846_G=0,rs5997301_G=1,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs4822939_G=0,rs5752634_G=0,rs2073780_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45443199_T=0,rs12166473_G=1,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs5752639_G=0,rs5752640_A=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0" ~ $3){print $1,$2,$3,$4,$5}' ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr22.28145327-28270372.haplotypes.hg38.vcf # h38 rows of vcf file

# chr22   27749339        rs41277849      G       A
# chr22   27750289        rs41277851      G       A
# chr22   27750747        rs45583937      G       A
# chr22   27753172        rs5762341       A       C
# chr22   27753786        rs118009109     C       T
# chr22   27754121        rs4822935       A       G
# chr22   27754827        rs2283844       G       A
# chr22   27755685        rs4820776       A       G
# chr22   27755837        rs2267106       G       A
# chr22   27757871        rs75851536      C       A
# chr22   27759416        rs2267113       C       T
# chr22   27760800        rs2014274       A       G
# chr22   27762426        rs8136505       C       T
# chr22   27763546        rs74562526      T       G
# chr22   27765715        rs16985531      A       G
# chr22   27767358        rs7292322       G       A
# chr22   27769113        rs17477193      T       G
# chr22   27769346        rs2283846       A       G
# chr22   27770414        rs5997301       C       G
# chr22   27772264        rs2073784       C       T
# chr22   27773533        rs1297593       C       A
# chr22   27774178        rs1807510       A       C
# chr22   27774229        rs2005397       C       T
# chr22   27775522        rs8142823       C       T
# chr22   27776589        rs4822939       T       G
# chr22   27777582        rs5752634       A       G
# chr22   27779174        rs2073780       C       T
# chr22   27779244        rs45552035      T       C
# chr22   27782526        rs45512704      C       T
# chr22   27783835        rs1297595       T       C
# chr22   27785448        rs382819        A       C
# chr22   27786522        rs424345        G       A
# chr22   27786548        rs45510597      G       A
# chr22   27787665        rs45443199      G       T
# chr22   27789464        rs12166473      T       G
# chr22   27791408        rs45462093      A       G
# chr22   27791785        rs45574432      T       C
# chr22   27796994        rs45589739      C       G
# chr22   27798759        rs45554336      C       T
# chr22   27799727        rs201186821     C       A
# chr22   27800488        rs200030766     C       T
# chr22   27802203        rs45540534      T       G
# chr22   27804188        rs5752639       A       G
# chr22   27807762        rs5752640       G       A
# chr22   27808809        rs11704546      T       C
# chr22   27810924        rs11705555      A       C
# chr22   27813879        rs117447690     G       A
# chr22   27818698        rs16985573      C       T
# chr22   27874384        rs5997320       T       G


# discovery
sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs118009109_T=0,rs4822935_G=1,rs2283844_A=1,rs4820776_G=1,rs2267106_A=0,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs8136505_T=0,rs74562526_G=0,rs16985531_G=0,rs7292322_A=0,rs17477193_G=0,rs2283846_G=0,rs5997301_G=1,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs4822939_G=0,rs5752634_G=0,rs2073780_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45443199_T=0,rs12166473_G=1,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs5752639_G=0,rs5752640_A=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0" "" "" "" "" h38.discovery
# LRT p value = 4.815529e-12 on 1 degree(s) of freedom with HR = 2.885349e+00 (2.233563e+00 to 3.727334e+00), frequency = 1.187622e-03 (430) = 3.273776e-03/1.078344e-03 (59/371), p = 2.938599e-12

# DRIVE replication
sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs118009109_T=0,rs4822935_G=1,rs2283844_A=1,rs4820776_G=1,rs2267106_A=0,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs8136505_T=0,rs74562526_G=0,rs16985531_G=0,rs7292322_A=0,rs17477193_G=0,rs2283846_G=0,rs5997301_G=1,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs4822939_G=0,rs5752634_G=0,rs2073780_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45443199_T=0,rs12166473_G=1,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs5752639_G=0,rs5752640_A=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0" "" "" "" "" h38.replication
# LRT p value = 1.307269e-05 on 1 degree(s) of freedom with OR = 1.915797e+00 (1.415681e+00 to 2.592589e+00), frequency = 1.788747e-03 (198) = 2.278473e-03/1.206392e-03 (137/61), p = 2.253350e-05

# LD with lead haplotype h38
sh /home/wletsou/scripts/ukbb_haplotype_haplotype_LD2.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0" "rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs118009109_T=0,rs4822935_G=1,rs2283844_A=1,rs4820776_G=1,rs2267106_A=0,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs8136505_T=0,rs74562526_G=0,rs16985531_G=0,rs7292322_A=0,rs17477193_G=0,rs2283846_G=0,rs5997301_G=1,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs4822939_G=0,rs5752634_G=0,rs2073780_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45443199_T=0,rs12166473_G=1,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs5752639_G=0,rs5752640_A=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6
# Haplotype 0 vs. haplotype 0:
# 0.00200515      0.00118762      0.00118762
# rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0     rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs118009109_T=0,rs4822935_G=1,rs2283844_A=1,rs4820776_G=1,rs2267106_A=0,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs8136505_T=0,rs74562526_G=0,rs16985531_G=0,rs7292322_A=0,rs17477193_G=0,rs2283846_G=0,rs5997301_G=1,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs4822939_G=0,rs5752634_G=0,rs2073780_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45443199_T=0,rs12166473_G=1,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs5752639_G=0,rs5752640_A=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0   0.591802   1

# LD with GWAS hit
mkdir ukbb_imputed_snps
cp dbgap28544_cases.indiv ukbb_imputed_snps/dbgap28544_cases.indiv
cp dbgap28544_controls.indiv ukbb_imputed_snps/dbgap28544_controls.indiv
cp ukbb_genotyped_snps/ukbb_bca_cases.indiv ukbb_imputed_snps/ukbb_bca_cases.indiv
cp ukbb_genotyped_snps/ukbb_bca_controls.indiv ukbb_imputed_snps/ukbb_bca_controls.indiv

cd ukbb_imputed_snps

# vcf file for UKBB
bsub -P SJLIFE -J ukbb_topmed_merge.chr22.28145327-28270372 -oo ukbb_topmed_merge.chr22.28145327-28270372.out -eo ukbb_topmed_merge.chr22.28145327-28270372.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_topmed_merge.sh \"rs41277849,rs41277851,rs45583937,rs5762341,rs118009109,rs4822935,rs2283844,rs4820776,rs2267106,rs75851536,rs2267113,rs2014274,rs8136505,rs74562526,rs16985531,rs7292322,rs17477193,rs2283846,rs5997301,rs2073784,rs1297593,rs1807510,rs2005397,rs8142823,rs4822939,rs5752634,rs2073780,rs45552035,rs45512704,rs1297595,rs16985561,rs16985561,rs382819,rs424345,rs45510597,rs45443199,rs12166473,rs45462093,rs45574432,rs45589739,rs45554336,rs201186821,rs200030766,rs45540534,rs5752639,rs5752640,rs11704546,rs11705555,rs117447690,rs16985573,rs5997320,rs62237573\" /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr22.qced.anno.info 22 28145327,28270372 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/raw/chr22 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps /home/wletsou/scripts"

# vcf file for dbGaP
bsub -P SJLIFE -J dbgap_haplotype_merge.chr22.28145327-28270372 -oo dbgap_haplotype_merge.chr22.28145327-28270372.out -eo dbgap_haplotype_merge.chr22.28145327-28270372.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/dbgap_haplotype_merge.sh dbgap28544_cases,dbgap28544_controls /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/chr22.qced.info \"rs41277849,rs41277851,rs45583937,rs5762341,rs118009109,rs4822935,rs2283844,rs4820776,rs2267106,rs75851536,rs2267113,rs2014274,rs8136505,rs74562526,rs16985531,rs7292322,rs17477193,rs2283846,rs5997301,rs2073784,rs1297593,rs1807510,rs2005397,rs8142823,rs4822939,rs5752634,rs2073780,rs45552035,rs45512704,rs1297595,rs16985561,rs16985561,rs382819,rs424345,rs45510597,rs45443199,rs12166473,rs45462093,rs45574432,rs45589739,rs45554336,rs201186821,rs200030766,rs45540534,rs5752639,rs5752640,rs11704546,rs11705555,rs117447690,rs16985573,rs5997320,rs62237573\" 22 28145327,28270372 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/vcfs/dbgap28544.topmedr2.cleaned.hg38.chr22.vcf.gz"

# extract imputed SNPs in range
bsub -P SJLIFE -J ukbb_haplotype_extract3.chr22.28145327-28270372 -oo ukbb_haplotype_extract3.chr22.28145327-28270372.out -eo ukbb_haplotype_extract3.chr22.28145327-28270372.err -R "rusage[mem=20000]" -q large_mem "sh /home/wletsou/scripts/ukbb_haplotype_extract3.sh ukbb_bca_cases.indiv,ukbb_bca_controls.indiv,dbgap28544_cases.indiv,dbgap28544_controls.indiv \"rs41277849,rs41277851,rs45583937,rs5762341,rs118009109,rs4822935,rs2283844,rs4820776,rs2267106,rs75851536,rs2267113,rs2014274,rs8136505,rs74562526,rs16985531,rs7292322,rs17477193,rs2283846,rs5997301,rs2073784,rs1297593,rs1807510,rs2005397,rs8142823,rs4822939,rs5752634,rs2073780,rs45552035,rs45512704,rs1297595,rs16985561,rs16985561,rs382819,rs424345,rs45510597,rs45443199,rs12166473,rs45462093,rs45574432,rs45589739,rs45554336,rs201186821,rs200030766,rs45540534,rs5752639,rs5752640,rs11704546,rs11705555,rs117447690,rs16985573,rs5997320,rs62237573\" /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr22.qced.anno.info 22 28145327,28270372 ukbb.topmed.hg19_chr22.28145327-28270372.hg38.vcf.gz,dbgap28544_cases+dbgap28544_controls.hg19_chr22.28145327-28270372.hg38.vcf.gz /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps /home/wletsou/scripts"

awk 'FILENAME==ARGV[1]{seen[$1]; next} FILENAME==ARGV[2]{seen[$1]; next} ($1 in seen || FNR==1){print $0}' ukbb_bca_cases.indiv ukbb_bca_controls.indiv haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.chr22.28145327-28270372.txt > haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt

awk 'FILENAME==ARGV[1]{seen[$1]; next} FILENAME==ARGV[2]{seen[$1]; next} ($1 in seen || FNR==1){print $0}' dbgap28544_cases.indiv dbgap28544_controls.indiv haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.chr22.28145327-28270372.txt > haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr22.28145327-28270372.txt

# LD, h1 (reduced) vs GWAS hit
sh /home/wletsou/scripts/ukbb_haplotype_haplotype_LD2.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0" "rs62237573_T=1" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps
# Haplotype 0 vs. haplotype 0:
# 0.00200515      0.00980755      0.00117105
# rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0     rs62237573_T=1  0.068216  0.579902

# effect of GWAS hit
sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs62237573_T=1" "" "" "" "" rs62237573.ukbb.discovery /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps /home/wletsou/scripts
# LRT p value = 1.674185e-05 on 1 degree(s) of freedom with HR = 1.347548e+00 (1.183747e+00 to 1.534015e+00), frequency = 9.807550e-03 (3551) = 1.292864e-02/9.644059e-03 (233/3318), p = 2.744970e-05

sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr22.28145327-28270372.txt "rs62237573_T=1" "" "" "" "" rs62237573.drive.replication /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps /home/wletsou/scripts
# LRT p value = 9.836504e-06 on 1 degree(s) of freedom with OR = 1.330243e+00 (1.170819e+00 to 1.511373e+00), frequency = 9.187656e-03 (1017) = 1.039449e-02/7.752551e-03 (625/392), p = 4.386159e-06

# effect of GWAS hit conditioned on h1 
sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs62237573_T=1" "1" "rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0" "1" "" h1.rs62237573.ukbb.discovery /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps /home/wletsou/scripts
# LRT p value = 3.063318e-02 on 1 degree(s) of freedom with HR = 3.078682e+00 (2.398129e+00 to 3.952365e+00), frequency = 1.171051e-03 (424) = 3.440240e-03/1.052185e-03 (62/362)

# effect of h1 conditioned on GWAS hit 
sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0" "1" "rs62237573_T=1" "1" "" rs62237573.h1.ukbb.discovery /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps /home/wletsou/scripts
# LRT p value = 4.505531e-10 on 1 degree(s) of freedom with HR = 3.076402e+00 (2.396354e+00 to 3.949437e+00), frequency = 1.171051e-03 (424) = 3.440240e-03/1.052185e-03 (62/362))

# effect of GWAS hit together with h1 
sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0:rs62237573_T=1" "" "" "" "" h1+rs62237573.ukbb.discovery /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps /home/wletsou/scripts
# LRT p value = 3.941604e-02 on 1 degree(s) of freedom with HR = 1.162131e+00 (1.010267e+00 to 1.336824e+00), frequency = 9.807550e-03 (3551) = 1.292864e-02/9.644059e-03 (233/3318)

sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs62237573_T=1:rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0" "" "" "" "" rs62237573+h1.ukbb.discovery /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps /home/wletsou/scripts
# LRT p value = 1.237728e-11 on 1 degree(s) of freedom with HR = 2.363786e+00 (1.888377e+00 to 2.958883e+00), frequency = 2.005148e-03 (726) = 4.993896e-03/1.848590e-03 (90/636)

# GWAS hit joined to h1
sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0,rs62237573_T=1" "" "" "" "" rs62237573,h1.ukbb.discovery /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps /home/wletsou/scripts
# LRT p value = 9.886850e-14 on 1 degree(s) of freedom with HR = 3.074147e+00 (2.394603e+00 to 3.946532e+00), frequency = 1.171051e-03 (424) = 3.440240e-03/1.052185e-03 (62/362), p = 4.806714e-14

sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0,rs62237573_T=1" "" "" "" "" rs62237573,h1.drive.replication /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps /home/wletsou/scripts
# LRT p value = 3.873884e-06 on 1 degree(s) of freedom with OR = 2.018460e+00 (1.480257e+00 to 2.752347e+00), frequency = 1.743577e-03 (193) = 2.261841e-03/1.127284e-03 (136/57), p = 6.132735e-06

# LD of h1 + GWAS hit with h38
sh /home/wletsou/scripts/ukbb_haplotype_haplotype_LD2.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0,rs62237573_T=1" "rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs118009109_T=0,rs4822935_G=1,rs2283844_A=1,rs4820776_G=1,rs2267106_A=0,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs8136505_T=0,rs74562526_G=0,rs16985531_G=0,rs7292322_A=0,rs17477193_G=0,rs2283846_G=0,rs5997301_G=1,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs4822939_G=0,rs5752634_G=0,rs2073780_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45443199_T=0,rs12166473_G=1,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs5752639_G=0,rs5752640_A=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps
# Haplotype 0 vs. haplotype 0:
# 0.00117105      0.00118762      0.00106886
# rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0,rs62237573_T=1      rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs118009109_T=0,rs4822935_G=1,rs2283844_A=1,rs4820776_G=1,rs2267106_A=0,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs8136505_T=0,rs74562526_G=0,rs16985531_G=0,rs7292322_A=0,rs17477193_G=0,rs2283846_G=0,rs5997301_G=1,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs4822939_G=0,rs5752634_G=0,rs2073780_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45443199_T=0,rs12166473_G=1,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs5752639_G=0,rs5752640_A=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0      0.821262        0.912632

# LD of GWAS hit with h38
sh /home/wletsou/scripts/ukbb_haplotype_haplotype_LD2.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs62237573_T=1" "rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs118009109_T=0,rs4822935_G=1,rs2283844_A=1,rs4820776_G=1,rs2267106_A=0,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs8136505_T=0,rs74562526_G=0,rs16985531_G=0,rs7292322_A=0,rs17477193_G=0,rs2283846_G=0,rs5997301_G=1,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs4822939_G=0,rs5752634_G=0,rs2073780_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45443199_T=0,rs12166473_G=1,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs5752639_G=0,rs5752640_A=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps
# Haplotype 0 vs. haplotype 0:
# 0.00980755      0.00118762      0.00106886
# rs62237573_T=1  rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs118009109_T=0,rs4822935_G=1,rs2283844_A=1,rs4820776_G=1,rs2267106_A=0,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs8136505_T=0,rs74562526_G=0,rs16985531_G=0,rs7292322_A=0,rs17477193_G=0,rs2283846_G=0,rs5997301_G=1,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs4822939_G=0,rs5752634_G=0,rs2073780_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45443199_T=0,rs12166473_G=1,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs5752639_G=0,rs5752640_A=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0      0.0970246       0.89901


# effect of h1,rs62237573_T=1 and h38
sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0,rs62237573_T=1:rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs118009109_T=0,rs4822935_G=1,rs2283844_A=1,rs4820776_G=1,rs2267106_A=0,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs8136505_T=0,rs74562526_G=0,rs16985531_G=0,rs7292322_A=0,rs17477193_G=0,rs2283846_G=0,rs5997301_G=1,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs4822939_G=0,rs5752634_G=0,rs2073780_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45443199_T=0,rs12166473_G=1,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs5752639_G=0,rs5752640_A=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0" "" "" "" "" h1,rs62237573+h38.ukbb.discovery /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps /home/wletsou/scripts

# effect of h1, and h38
sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0:rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs118009109_T=0,rs4822935_G=1,rs2283844_A=1,rs4820776_G=1,rs2267106_A=0,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs8136505_T=0,rs74562526_G=0,rs16985531_G=0,rs7292322_A=0,rs17477193_G=0,rs2283846_G=0,rs5997301_G=1,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs4822939_G=0,rs5752634_G=0,rs2073780_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45443199_T=0,rs12166473_G=1,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs5752639_G=0,rs5752640_A=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0" "" "" "" "" h1+h38.ukbb.discovery /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps /home/wletsou/scripts

# effect of GWAS hit and h38
sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr22.28145327-28270372.txt "rs62237573_T=1:rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs118009109_T=0,rs4822935_G=1,rs2283844_A=1,rs4820776_G=1,rs2267106_A=0,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs8136505_T=0,rs74562526_G=0,rs16985531_G=0,rs7292322_A=0,rs17477193_G=0,rs2283846_G=0,rs5997301_G=1,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs4822939_G=0,rs5752634_G=0,rs2073780_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45443199_T=0,rs12166473_G=1,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs5752639_G=0,rs5752640_A=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0" "" "" "" "" rs62237573+h38.ukbb.discovery /scratch_space/wletsou/sjlife/GWAS/UKBB_chr22.6/ukbb_imputed_snps /home/wletsou/scripts

sh /home/wletsou/scripts/vcf2bed.v3.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr22.28145327-28270372.haplotypes.hg38.vcf "rs62237573_T=1" "" "rs62237573.bed" "0" # rs62237573, with alleles, no header

awk 'BEGIN{OFS="\t"} {$9="0,0,0"}1' rs62237573.bed > rs62237573.tmp && mv rs62237573.tmp rs62237573.bed # blue to black

sh /home/wletsou/scripts/vcf2bed.v3.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr22.28145327-28270372.haplotypes.hg38.vcf "rs41277849_A=0,rs5762341_C=0,rs4822935_G=1,rs4820776_G=1,rs2267113_T=0,rs2014274_G=0,rs74562526_G=0,rs7292322_A=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45462093_G=0,rs45574432_C=0,rs45540534_G=0,rs11704546_C=1,rs11705555_C=0,rs5997320_G=0" "" "h1.bed" "0" # h1, with alleles, no header

sh /home/wletsou/scripts/vcf2bed.v3.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr22.28145327-28270372.haplotypes.hg38.vcf "rs41277849_A=0,rs41277851_A=0,rs45583937_A=0,rs5762341_C=0,rs118009109_T=0,rs4822935_G=1,rs2283844_A=1,rs4820776_G=1,rs2267106_A=0,rs75851536_A=0,rs2267113_T=0,rs2014274_G=0,rs8136505_T=0,rs74562526_G=0,rs16985531_G=0,rs7292322_A=0,rs17477193_G=0,rs2283846_G=0,rs5997301_G=1,rs2073784_T=0,rs1297593_A=1,rs1807510_C=0,rs2005397_T=0,rs8142823_T=0,rs4822939_G=0,rs5752634_G=0,rs2073780_T=0,rs45552035_C=0,rs45512704_T=0,rs1297595_C=1,rs382819_C=1,rs424345_A=0,rs45510597_A=0,rs45443199_T=0,rs12166473_G=1,rs45462093_G=0,rs45574432_C=0,rs45589739_G=0,rs45554336_T=0,rs201186821_A=0,rs200030766_T=0,rs45540534_G=0,rs5752639_G=0,rs5752640_A=0,rs11704546_C=1,rs11705555_C=0,rs117447690_A=0,rs16985573_T=0,rs5997320_G=0" "" "h38.bed" "0" # h38, with alleles, no header

module load ucsc/112922

fetchChromSizes hg38 > hg38.sizes # chrom sizes for bigBed file

# convert bed to bigBed
bedToBigBed rs62237573.bed hg38.sizes rs62237573.bb
bedToBigBed h1.bed hg38.sizes h1.bb
bedToBigBed h38.bed hg38.sizes h38.bb

# move to ClusterHome
cp rs62237573.bed /home/wletsou/Chromosome_Overlap_results/UKBB_chr22.6/ukbb_imputed_snps/rs62237573.bed
cp h1.bed /home/wletsou/Chromosome_Overlap_results/UKBB_chr22.6/ukbb_imputed_snps/h1.bed
cp h38.bed /home/wletsou/Chromosome_Overlap_results/UKBB_chr22.6/ukbb_imputed_snps/h38.bed

cp rs62237573.bb /home/wletsou/Chromosome_Overlap_results/UKBB_chr22.6/ukbb_imputed_snps/rs62237573.bb
cp h1.bb /home/wletsou/Chromosome_Overlap_results/UKBB_chr22.6/ukbb_imputed_snps/h1.bb
cp h38.bb /home/wletsou/Chromosome_Overlap_results/UKBB_chr22.6/ukbb_imputed_snps/h38.bb