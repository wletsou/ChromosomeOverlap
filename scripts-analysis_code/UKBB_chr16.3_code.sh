mkdir /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3
cd /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3

# extract UKBB-genotyped SNPs in entire TAD

mkdir /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/ukbb_genotyped_snps
cd ukbb_genotyped_snps

awk 'BEGIN{OFS="\t"} ($3==1){print $1,$1}' /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_cases.indiv
awk 'BEGIN{OFS="\t"} ($3==0){print $1,$1}' /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_controls.indiv

bsub -P SJLIFE -J ukbb_hybrid_haplotype3.chr16.52000000-54000000 -oo ukbb_hybrid_haplotype3.chr16.52000000-54000000.out -eo ukbb_hybrid_haplotype3.chr16.52000000-54000000.err -R "rusage[mem=20000]" "sh /home/wletsou/scripts/ukbb_hybrid_haplotype3.sh ukbb_bca_cases.indiv,ukbb_bca_controls.indiv 16 52000000,54000000 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/biobank/gds/ukb.bca.hap/ukb.bca.hap.chr16.new.vcf.gz"

cd /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3

# vcf and haplotypes for Phase 1

# get genotyped SNPs in region
awk 'BEGIN{OFS="\t"} ($2>=52519188 && $2<=52679188){print $3}' /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/ukbb_genotyped_snps/ukbb_bca_cases+ukbb_bca_controls.chr16.52000000-54000000.haplotypes.vcf > ukbb_snp_list.chr16.52519188-52679188.txt

# check that SNPs have good "imputation" quality
awk 'NR==FNR{snp[$1]; next} ($8 in snp && $15==1){print $8}' ukbb_snp_list.chr16.52519188-52679188.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr16.qced.anno.info > ukbb_snp_list.chr16.52519188-52679188.tmp
test -f ukbb_snp_list.chr16.52519188-52679188.tmp && mv ukbb_snp_list.chr16.52519188-52679188.tmp ukbb_snp_list.chr16.52519188-52679188.txt

awk 'BEGIN{OFS="\t"} ($3==1){print $1,$1}' /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_cases.indiv
awk 'BEGIN{OFS="\t"} ($3==0){print $1,$1}' /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_controls.indiv

# vcf file for UKBB
bsub -P SJLIFE -J ukbb_topmed_merge.chr16.52519188-52679188 -oo ukbb_topmed_merge.chr16.52519188-52679188.out -eo ukbb_topmed_merge.chr16.52519188-52679188.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_topmed_merge.sh ukbb_snp_list.chr16.52519188-52679188.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr16.qced.anno.info 16 52519188,52679188 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/raw/chr16"

# vcf file for dbGaP
bsub -P SJLIFE -J dbgap_haplotype_merge.chr16.52519188-52679188 -oo dbgap_haplotype_merge.chr16.52519188-52679188.out -eo dbgap_haplotype_merge.chr16.52519188-52679188.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/dbgap_haplotype_merge.sh dbgap28544_cases,dbgap28544_controls /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/chr16.qced.info ukbb_snp_list.chr16.52519188-52679188.txt 16 52519188,52679188 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/vcfs/dbgap28544.topmedr2.cleaned.hg38.chr16.vcf.gz"

# extract haplotypes for cases and controls
bsub -P SJLIFE -J ukbb_bca_overlap.chr16.52519188-52679188 -oo ukbb_bca_overlap.chr16.52519188-52679188.out -eo ukbb_bca_overlap.chr16.52519188-52679188.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_submission.v2.sh ukbb_bca_cases,ukbb_bca_controls dbgap28544_cases,dbgap28544_controls \"\" 16 52519188,52679188 \"\" \"\" ukbb_snp_list.chr16.52519188-52679188.txt ukbb.topmed.hg19_chr16.52519188-52679188.hg38.vcf.gz,dbgap28544_cases+dbgap28544_controls.hg19_chr16.52519188-52679188.hg38.vcf.gz /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr16.qced.anno.info \"\" \"\" \"\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3 /home/wletsou/scripts"

awk 'NR==1{char=length(NF-1)} NR>1{printf "1\t"; for (i=2;i<=NF;i++) {printf "%s%0"char"d_%s",(i>2?",":""),i-1,$i}; printf "\n"}' haplotype_estimates.ukbb_bca_cases.chr16.52519188-52679188.txt | awk 'BEGIN{OFS="\t"} {seen[$2]+=$1} END{for (i in seen) {print seen[i],i} }' > Pattern_combined.Iteration000.chr16.52519188-52679188_2,j.txt # starting with all chromosomes, not just pairs

# First round of overlaps:
bsub -P SJLIFE -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr16.52519188-52679188 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr16.52519188-52679188.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr16.52519188-52679188.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ChromosomeOverlap_iteration_sub_parallel.v3.sh chr16.52519188-52679188 2 2,j 50 0 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3 /home/wletsou/scripts"

# p-values for first round of overlaps
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr16.52519188-52679188.txt,haplotype_estimates.ukbb_bca_controls.chr16.52519188-52679188.txt Closed_patterns_stats.chr16.52519188-52679188_2,j.txt 50 \"1\" \"Iteration001.chr16.52519188-52679188\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3 /home/wletsou/scripts"

# remove older results if staring again
rm -r *Iteration00[2-9]*
awk '($1<2){print $0}' Closed_patterns_stats.chr16.52519188-52679188_2,j.txt > Closed_patterns_stats.chr16.52519188-52679188_2,j.tmp
mv Closed_patterns_stats.chr16.52519188-52679188_2,j.tmp Closed_patterns_stats.chr16.52519188-52679188_2,j.txt

test -f Pattern_combined_old.Iteration001.chr16.52519188-52679188_2,j.txt && mv Pattern_combined_old.Iteration001.chr16.52519188-52679188_2,j.txt Pattern_combined.Iteration001.chr16.52519188-52679188_2,j.txt
awk 'NR==FNR{seen[$1]; next} ($2 in seen){print $0}' <(awk '(NR>1 && $5+0<1e-16 && $4+0>1){print $1}' fisher_exact.Iteration001.chr16.52519188-52679188.patterns_0000001-3458020.txt) Pattern_combined.Iteration001.chr16.52519188-52679188_2,j.txt > Pattern_combined.Iteration001.chr16.52519188-52679188_2,j.tmp
mv Pattern_combined.Iteration001.chr16.52519188-52679188_2,j.txt Pattern_combined_old.Iteration001.chr16.52519188-52679188_2,j.txt
mv Pattern_combined.Iteration001.chr16.52519188-52679188_2,j.tmp Pattern_combined.Iteration001.chr16.52519188-52679188_2,j.txt

# Begin iterations:
bsub -P SJLIFE -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr16.52519188-52679188 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr16.52519188-52679188.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr16.52519188-52679188.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ChromosomeOverlap_iteration_sub_parallel.v3.sh chr16.52519188-52679188 2 2,j 50 1 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3 /home/wletsou/scripts"

# (filtered) patterns appearing for the first time at Iterations 001-005
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1){print $0}' Pattern_combined.Iteration001.chr16.52519188-52679188_2,j.txt Closed_patterns_stats.chr16.52519188-52679188_2,j.txt > Closed_patterns_stats.Iteration001-005.chr16.52519188-52679188_2,j.txt

# p-values for Iterations 001-005
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr16.52519188-52679188.txt,haplotype_estimates.ukbb_bca_controls.chr16.52519188-52679188.txt Closed_patterns_stats.Iteration001-005.chr16.52519188-52679188_2,j.txt 50 \"\" \"Iteration001-005.chr16.52519188-52679188\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3 /home/wletsou/scripts"

# top closed patterns
awk '($5+0<1e-30 && $4+0>1){print $0}'  fisher_exact.Iteration001-005.chr16.52519188-52679188.patterns_000001-269933.txt >  fisher_exact.Iteration001-005.chr16.52519188-52679188.patterns_000001-269933.Results.txt

# tranlate haplotype to rsid
bsub -P SJLIFE -J ukbb_haplotype_translate2_sub -oo ukbb_haplotype_translate2_sub.out -eo ukbb_haplotype_translate2_sub.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr16.52519188-52679188.txt fisher_exact.Iteration001-005.chr16.52519188-52679188.patterns_000001-269933.Results.txt 50"

# precompute haplotype counts
bsub -P SJLIFE -J ukbb_haplotype_model9_sub.ukbb_discovery -oo ukbb_haplotype_model9_sub.v2.ukbb_discovery.out -eo ukbb_haplotype_model9_sub.v2.ukbb_discovery.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52519188-52679188.txt fisher_exact.Iteration001-005.chr16.52519188-52679188.patterns_000001-269933.Results.translated.txt \"\" \"\" \"\" 16 52519188,52679188 50 \"ukbb_discovery\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3 /home/wletsou/scripts"

bsub -P SJLIFE -J ukbb_haplotype_model9_sub.drive_replication -oo ukbb_haplotype_model9_sub.v2.drive_replication.out -eo ukbb_haplotype_model9_sub.v2.drive_replication.err -R "rusage[mem=256]" -w "done(179460547)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52519188-52679188.txt fisher_exact.Iteration001-005.chr16.52519188-52679188.patterns_000001-269933.Results.translated.txt \"\" \"\" \"\" 16 52519188,52679188 50 \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3 /home/wletsou/scripts"

# forward selection
bsub -P SJLIFE -J ukbb_haplotype_model9_iterate.v2.discovery -oo ukbb_haplotype_model9_iterate.v2.discovery.out -eo ukbb_haplotype_model9_iterate.v2.discovery.err -R "rusage[mem=256]" -w "done(179460557)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.v2.sh ukbb_discovery.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt \"\" \"\" 1e-5 ukbb_discovery.Significant_patterns.txt \"ukbb_discovery\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3 /home/wletsou/scripts"

bsub -P SJLIFE -J dbgap28544_haplotype_model9_iterate.v2.replication -oo dbgap28544_haplotype_model9_iterate.v2.replication.out -eo dbgap28544_haplotype_model9_iterate.v2.replication.err -w "done(179460560)" -R "rusage[mem=256]" "sh /home/wletsou/scripts/dbgap28544_haplotype_model9_iterate.v2.sh drive_replication.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt \"\" \"\" 1e-5 drive_replication.Significant_patterns.txt \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3 /home/wletsou/scripts"

# Including contiguous haplotypes (Iteration000)

# (filtered) patterns appearing for the first time at Iterations 000,001-005
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1 || $1==0){print $0}' Pattern_combined.Iteration001.chr16.52519188-52679188_2,j.txt Closed_patterns_stats.chr16.52519188-52679188_2,j.txt > Closed_patterns_stats.Iteration000+Iteration001-005.chr16.52519188-52679188_2,j.txt

# p-values for Iterations 000,001-005
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr16.52519188-52679188.txt,haplotype_estimates.ukbb_bca_controls.chr16.52519188-52679188.txt Closed_patterns_stats.Iteration000+Iteration001-005.chr16.52519188-52679188_2,j.txt 50 \"\" \"Iteration000+Iteration001-005.chr16.52519188-52679188\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3 /home/wletsou/scripts"

# top closed patterns
awk '($5+0<1e-30 && $4+0>1){print $0}'  fisher_exact.Iteration000+Iteration001-005.chr16.52519188-52679188.patterns_000001-272791.txt >  fisher_exact.Iteration000+Iteration001-005.chr16.52519188-52679188.patterns_000001-272791.Results.txt

# minimum p-value for contiguous patterns
awk 'BEGIN{min=1} NR==FNR{seen[$2]; next} ($1 in seen){ if ($5<min+0) {min=$5} } END{print min}' Pattern_combined.Iteration000.chr16.52519188-52679188_2,j.txt fisher_exact.Iteration000+Iteration001-005.chr16.52519188-52679188.patterns_000001-272791.txt
# 0.000321207713384453

# tranlate haplotype to rsid
bsub -P SJLIFE -J ukbb_haplotype_translate2_sub.Iteration000+Iteration001-005 -oo ukbb_haplotype_translate2_sub.Iteration000+Iteration001-005.out -eo ukbb_haplotype_translate2_sub.Iteration000+Iteration001-005.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr16.52519188-52679188.txt fisher_exact.Iteration000+Iteration001-005.chr16.52519188-52679188.patterns_000001-272791.Results.txt 50"

# precompute haplotype counts
bsub -P SJLIFE -J ukbb_haplotype_model9_sub.ukbb_discovery.Iteration000+Iteration001-005 -oo ukbb_haplotype_model9_sub.v2.ukbb_discovery.Iteration000+Iteration001-005.out -eo ukbb_haplotype_model9_sub.v2.ukbb_discovery.Iteration000+Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52519188-52679188.txt fisher_exact.Iteration000+Iteration001-005.chr16.52519188-52679188.patterns_000001-272791.Results.translated.txt \"\" \"\" \"\" 16 52519188,52679188 50 \"ukbb_discovery.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3 /home/wletsou/scripts"

bsub -P SJLIFE -J ukbb_haplotype_model9_sub.drive_replication.Iteration000+Iteration001-005 -oo ukbb_haplotype_model9_sub.v2.drive_replication.Iteration000+Iteration001-005.out -eo ukbb_haplotype_model9_sub.v2.drive_replication.Iteration000+Iteration001-005.err -R "rusage[mem=256]" -w "done(179460547)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52519188-52679188.txt fisher_exact.Iteration000+Iteration001-005.chr16.52519188-52679188.patterns_000001-272791.Results.translated.txt \"\" \"\" \"\" 16 52519188,52679188 50 \"drive_replication.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3 /home/wletsou/scripts"

# forward selection
bsub -P SJLIFE -J ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-005 -oo ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-005.out -eo ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-005.err -R "rusage[mem=256]" -w "done(179460557)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.v2.sh ukbb_discovery.Iteration000+Iteration001-005.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt \"\" \"\" 1e-5 ukbb_discovery.Iteration000+Iteration001-005.Significant_patterns.txt \"ukbb_discovery.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3 /home/wletsou/scripts"

bsub -P SJLIFE -J dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-005 -oo dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-005.out -eo dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-005.err -w "done(179460560)" -R "rusage[mem=256]" "sh /home/wletsou/scripts/dbgap28544_haplotype_model9_iterate.v2.sh drive_replication.Iteration000+Iteration001-005.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt \"\" \"\" 1e-5 drive_replication.Iteration000+Iteration001-005.Significant_patterns.txt \"drive_replication.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3 /home/wletsou/scripts"

awk -v alpha=0.05 'BEGIN{OFS = "\t"} NR==FNR && $2 > 1 && $5 < 1e-5{seen[$1]=$0; n+=1; next} ($1 in seen && $2 > 1 && $5 < alpha){m+=1} END{printf("%s of %s replicated at p < %0.2e\n",m + 0,n + 0,alpha)}' ukbb_discovery.Conditional_haplotype_effects.h1.txt drive_replication.Conditional_haplotype_effects.h1.txt
# 29078 of 29078 replicated at p < 5.00e-02

# reduction
sh /home/wletsou/scripts/haplotype_rpart.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52519188-52679188.txt "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs45625634_C=0,rs45496893_T=0,rs34750829_A=0,rs9931232_A=1,rs45477396_C=0,rs75571494_A=0,rs4594251_T=0,rs4784227_T=1,rs11864809_C=0,rs118162666_A=0,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "" 1 # h1
# reduced = rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0

sh /home/wletsou/scripts/vcf2bed.v3.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr16.52519188-52679188.haplotypes.hg38.vcf "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "" "h1.bed" "0" # h1, with alleles, no header

awk 'BEGIN{OFS="\t"} ("rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" ~ $3){print $1,$2,$3,$4,$5}' ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr16.52519188-52679188.haplotypes.hg38.vcf
# h1 rows of vcf file

# chr16   52486414        rs78338509      A       G
# chr16   52487716        rs12447430      T       C
# chr16   52494687        rs45577538      G       T
# chr16   52502436        rs8046985       G       A
# chr16   52506283        rs45542333      G       T
# chr16   52507901        rs78268044      C       T
# chr16   52514240        rs45454402      A       G
# chr16   52525814        rs45512493      A       G
# chr16   52535625        rs34750829      G       A
# chr16   52538920        rs9931232       G       A
# chr16   52563932        rs4594251       C       T
# chr16   52565276        rs4784227       C       T
# chr16   52603296        rs45584434      G       A
# chr16   52603720        rs45560737      C       A
# chr16   52604455        rs45575339      C       T
# chr16   52606294        rs45587544      G       A
# chr16   52626191        rs78841172      T       C
# chr16   52626967        rs76000465      T       C
# chr16   52633030        rs111925335     A       C
# chr16   52634678        rs75127968      G       A

# GWAS hit
sh /home/wletsou/scripts/vcf2bed.v3.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr16.52519188-52679188.haplotypes.hg38.vcf "rs4784227_T=1" "" "rs4784227.bed" "0" # rs4784227, with alleles, no header

awk 'BEGIN{OFS="\t"} {$9="0,0,0"}1' rs4784227.bed > rs4784227.tmp && mv rs4784227.tmp rs4784227.bed # blue to black

module load ucsc/112922

fetchChromSizes hg38 > hg38.sizes # chrom sizes for bigBed file

# convert bed to bigBed
bedToBigBed rs4784227.bed hg38.sizes rs4784227.bb
bedToBigBed h1.bed hg38.sizes h1.bb

# move to ClusterHome
cp rs4784227.bed /home/wletsou/Chromosome_Overlap_results/UKBB_chr16.3/rs4784227.bed
cp h1.bed /home/wletsou/Chromosome_Overlap_results/UKBB_chr16.3/h1.bed

cp rs4784227.bb /home/wletsou/Chromosome_Overlap_results/UKBB_chr16.3/rs4784227.bb
cp h1.bb /home/wletsou/Chromosome_Overlap_results/UKBB_chr16.3/h1.bb

# discovery
sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52519188-52679188.txt "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "" "" "" "" h1.discovery
# LRT p value = 1.234118e-41 on 1 degree(s) of freedom with HR = 1.265392e+00 (1.223759e+00 to 1.308441e+00), frequency = 2.139874e-01 (77478) = 2.549662e-01/2.118409e-01 (4595/72883), p = 1.898249e-41

# DRIVE replication
sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52519188-52679188.txt "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "" "" "" "" h1.replication
# LRT p value = 8.520011e-43 on 1 degree(s) of freedom with OR = 1.219175e+00 (1.184992e+00 to 1.254345e+00), frequency = 2.276407e-01 (25198) = 2.430814e-01/2.092793e-01 (14616/10582), p = 7.866135e-41

# discovery
sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52519188-52679188.txt "rs4784227_T=1" "" "" "" "" rs4784227.discovery
# LRT p value = 1.170830e-40 on 1 degree(s) of freedom with HR = 1.252320e+00 (1.212312e+00 to 1.293649e+00), frequency = 2.389303e-01 (86509) = 2.808234e-01/2.367358e-01 (5061/81448), p = 2.588138e-40

# DRIVE replication
sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52519188-52679188.txt "rs4784227_T=1" "" "" "" "" rs4784227.replication
# LRT p value = 1.808780e-50 on 1 degree(s) of freedom with OR = 1.229432e+00 (1.196436e+00 to 1.263338e+00), frequency = 2.644816e-01 (29276) = 2.822479e-01/2.433550e-01 (16971/12305), p = 1.646459e-48

# LD, h1 (reduced) vs GWAS hit
sh /home/wletsou/scripts/ukbb_haplotype_haplotype_LD2.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52519188-52679188.txt "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "rs4784227_T=1" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3
# Haplotype 0 vs. haplotype 0:
# 0.213987        0.23893 0.213987
# rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0       rs4784227_T=1   0.867186        1


mkdir phase2
cp ukbb_bca_cases.indiv phase2/ukbb_bca_cases.indiv
cp ukbb_bca_controls.indiv phase2/ukbb_bca_controls.indiv
cp dbgap28544_cases.indiv phase2/dbgap28544_cases.indiv
cp dbgap28544_controls.indiv phase2/dbgap28544_controls.indiv
cd phase2

# get h1 and list of SNPs in upstream region
HAPLOTYPE_SNPS=$(echo "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" | sed 's/_[A-Z]=[0-9]//g') # snps in h1, no alleles

awk 'BEGIN{OFS="\t"} ( ($2>=52100000 && $2<52519188) || ($2>52679188 && $2<=53100000) || "'$HAPLOTYPE_SNPS'" ~ $3){print $3}' /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/ukbb_genotyped_snps/ukbb_bca_cases+ukbb_bca_controls.chr16.52000000-54000000.haplotypes.vcf > ukbb_snp_list.chr16.52100000-53100000.txt

# check that SNPs have good "imputation" quality
awk 'NR==FNR{snp[$1]; next} ($8 in snp && $15==1){print $8}' ukbb_snp_list.chr16.52100000-53100000.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr16.qced.anno.info > ukbb_snp_list.chr16.52100000-53100000.tmp
test -f ukbb_snp_list.chr16.52100000-53100000.tmp && mv ukbb_snp_list.chr16.52100000-53100000.tmp ukbb_snp_list.chr16.52100000-53100000.txt

# vcf files of region (UKBB and dbGaP)
bsub -P SJLIFE -J ukbb_topmed_merge.chr16.52100000-53100000 -oo ukbb_topmed_merge.chr16.52100000-53100000.out -eo ukbb_topmed_merge.chr16.52100000-53100000.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_topmed_merge.sh ukbb_snp_list.chr16.52100000-53100000.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr16.qced.anno.info 16 52100000,53100000 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/raw/chr16"

bsub -P SJLIFE -J dbgap_haplotype_merge.chr16.52100000-53100000 -oo dbgap_haplotype_merge.chr16.52100000-53100000.out -eo dbgap_haplotype_merge.chr16.52100000-53100000.err -R "rusage[mem=1000]" -w "done(179596006)" "sh /home/wletsou/scripts/dbgap_haplotype_merge.sh dbgap28544_cases,dbgap28544_controls \"/research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt\" \"/research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/chr16.qced.info\" ukbb_snp_list.chr16.52100000-53100000.txt 16 52100000,53100000 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/vcfs/dbgap28544.topmedr2.cleaned.hg38.chr16.vcf.gz"

# extract haplotypes for cases and controls
bsub -P SJLIFE -J ukbb_bca_overlap.chr16.52100000-53100000 -oo ukbb_bca_overlap.chr16.52100000-53100000.out -eo ukbb_bca_overlap.chr16.52100000-53100000.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_submission.v2.sh ukbb_bca_cases,ukbb_bca_controls dbgap28544_cases,dbgap28544_controls \"\" 16 52100000,53100000 \"rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0\" 1 ukbb_snp_list.chr16.52100000-53100000.txt ukbb.topmed.hg19_chr16.52100000-53100000.hg38.vcf.gz,dbgap28544_cases+dbgap28544_controls.hg19_chr16.52100000-53100000.hg38.vcf.gz /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr16.qced.anno.info \"\" \"\" \"\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2 /home/wletsou/scripts"

awk 'NR==1{char=length(NF-1)} NR>1{printf "1\t"; for (i=2;i<=NF;i++) {printf "%s%0"char"d_%s",(i>2?",":""),i-1,$i}; printf "\n"}' haplotype_estimates.ukbb_bca_cases.chr16.52100000-53100000.subset.txt | awk 'BEGIN{OFS="\t"} {seen[$2]+=$1} END{for (i in seen) {print seen[i],i} }' > Pattern_combined.Iteration000.chr16.52100000-53100000_2,j.txt # starting with all chromosomes, not just pairs

# First round of overlaps:
bsub -P SJLIFE -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr16.52100000-53100000 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr16.52100000-53100000.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr16.52100000-53100000.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ChromosomeOverlap_iteration_sub_parallel.v3.sh chr16.52100000-53100000 2 2,j 50 0 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2 /home/wletsou/scripts"

# p-values for first round of overlaps
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr16.52100000-53100000.subset.txt,haplotype_estimates.ukbb_bca_controls.chr16.52100000-53100000.subset.txt Closed_patterns_stats.chr16.52100000-53100000_2,j.txt 50 \"1\" \"Iteration001.chr16.52100000-53100000\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2 /home/wletsou/scripts"

# remove older results if staring again
rm -r *Iteration00[2-9]*
awk '($1<2){print $0}' Closed_patterns_stats.chr16.52100000-53100000_2,j.txt > Closed_patterns_stats.chr16.52100000-53100000_2,j.tmp
mv Closed_patterns_stats.chr16.52100000-53100000_2,j.tmp Closed_patterns_stats.chr16.52100000-53100000_2,j.txt

test -f Pattern_combined_old.Iteration001.chr16.52100000-53100000_2,j.txt && mv Pattern_combined_old.Iteration001.chr16.52100000-53100000_2,j.txt Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.txt
awk 'NR==FNR{seen[$1]; next} ($2 in seen){print $0}' <(awk '(NR>1 && $5+0<1.5e-5 && $4+0>1){print $0}' fisher_exact.Iteration001.chr16.52100000-53100000.patterns_0000001-5954654.txt | awk 'BEGIN{array[0]} {delete array; n=split($1,array,","); if (n<snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]+0 || snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]+0==0) {snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]=n; seen[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]=$1} } END{for (i in seen) {print seen[i]} }') Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.txt > Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.tmp # take shortest pattern of a family with the same p-value
mv Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.txt Pattern_combined_old.Iteration001.chr16.52100000-53100000_2,j.txt
mv Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.tmp Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.txt

# Begin iterations:
bsub -P SJLIFE -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr16.52100000-53100000 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr16.52100000-53100000.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr16.52100000-53100000.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ChromosomeOverlap_iteration_sub_parallel.v3.sh chr16.52100000-53100000 2 2,j 50 1 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2 /home/wletsou/scripts"

# (filtered) patterns appearing for the first time at Iterations 001-009
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1){print $0}' Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.txt Closed_patterns_stats.chr16.52100000-53100000_2,j.txt > Closed_patterns_stats.Iteration001-009.chr16.52100000-53100000_2,j.txt

# p-values for Iterations 001-009
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-009 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-009.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-009.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr16.52100000-53100000.subset.txt,haplotype_estimates.ukbb_bca_controls.chr16.52100000-53100000.subset.txt Closed_patterns_stats.Iteration001-009.chr16.52100000-53100000_2,j.txt 50 \"\" \"Iteration001-009.chr16.52100000-53100000\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2 /home/wletsou/scripts"

# top closed patterns
awk '($5+0<1e-4 && $4+0>1){print $0}' fisher_exact.Iteration001-009.chr16.52100000-53100000.patterns_0000001-2517921.txt > fisher_exact.Iteration001-009.chr16.52100000-53100000.patterns_0000001-2517921.Results.txt

# tranlate haplotype to rsid
bsub -P SJLIFE -J ukbb_haplotype_translate2_sub -oo ukbb_haplotype_translate2_sub.out -eo ukbb_haplotype_translate2_sub.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr16.52100000-53100000.subset.txt fisher_exact.Iteration001-009.chr16.52100000-53100000.patterns_0000001-2517921.Results.txt 50"

# precompute haplotype counts
bsub -P SJLIFE -J ukbb_haplotype_model9_sub.ukbb_discovery -oo ukbb_haplotype_model9_sub.v2.ukbb_discovery.out -eo ukbb_haplotype_model9_sub.v2.ukbb_discovery.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52100000-53100000.txt fisher_exact.Iteration001-009.chr16.52100000-53100000.patterns_0000001-2517921.Results.translated.txt \"1\" \"rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0\" \"1\" 16 52100000,53100000 50 \"ukbb_discovery\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2 /home/wletsou/scripts"

bsub -P SJLIFE -J ukbb_haplotype_model9_sub.drive_replication -oo ukbb_haplotype_model9_sub.v2.drive_replication.out -eo ukbb_haplotype_model9_sub.v2.drive_replication.err -R "rusage[mem=256]" -w "done(179707180)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52100000-53100000.txt fisher_exact.Iteration001-009.chr16.52100000-53100000.patterns_0000001-2517921.Results.translated.txt \"1\" \"rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0\" \"1\" 16 52100000,53100000 50 \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2 /home/wletsou/scripts"

# forward selection
bsub -P SJLIFE -J ukbb_haplotype_model9_iterate.v2.discovery -oo ukbb_haplotype_model9_iterate.v2.discovery.out -eo ukbb_haplotype_model9_iterate.v2.discovery.err -R "rusage[mem=256]" -w "done(179707185)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.v2.sh ukbb_discovery.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt \"1\" \"1\" 1e-5 ukbb_discovery.Significant_patterns.txt \"ukbb_discovery\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2 /home/wletsou/scripts"

bsub -P SJLIFE -J dbgap28544_haplotype_model9_iterate.v2.replication -oo dbgap28544_haplotype_model9_iterate.v2.replication.out -eo dbgap28544_haplotype_model9_iterate.v2.replication.err -w "done(179707186)" -R "rusage[mem=256]" "sh /home/wletsou/scripts/dbgap28544_haplotype_model9_iterate.v2.sh drive_replication.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt \"1\" \"1\" 1e-5 drive_replication.Significant_patterns.txt \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2 /home/wletsou/scripts"

# Including contiguous haplotypes (Iteration000)

# (filtered) patterns appearing for the first time at Iterations 000,001-009
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1 || $1==0){print $0}' Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.txt Closed_patterns_stats.chr16.52100000-53100000_2,j.txt > Closed_patterns_stats.Iteration000+Iteration001-009.chr16.52100000-53100000_2,j.txt

# p-values for Iterations 000,001-009
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-009 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-009.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-009.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr16.52100000-53100000.subset.txt,haplotype_estimates.ukbb_bca_controls.chr16.52100000-53100000.subset.txt Closed_patterns_stats.Iteration000+Iteration001-009.chr16.52100000-53100000_2,j.txt 50 \"\" \"Iteration000+Iteration001-009.chr16.52100000-53100000\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2 /home/wletsou/scripts"

# top closed patterns
awk '($5+0<1e-4 && $4+0>1){print $0}' fisher_exact.Iteration000+Iteration001-009.chr16.52100000-53100000.patterns_0000001-2521373.txt > fisher_exact.Iteration000+Iteration001-009.chr16.52100000-53100000.patterns_0000001-2521373.Results.txt

# minimum p-value for contiguous patterns
awk 'BEGIN{min=1} NR==FNR{seen[$2]; next} ($1 in seen){ if ($5<min+0) {min=$5} } END{print min}' Pattern_combined.Iteration000.chr16.52100000-53100000_2,j.txt fisher_exact.Iteration000+Iteration001-009.chr16.52100000-53100000.patterns_0000001-2521373.txt
# 0.00141667698272622

# tranlate haplotype to rsid
bsub -P SJLIFE -J ukbb_haplotype_translate2_sub.Iteration000+Iteration001-009 -oo ukbb_haplotype_translate2_sub.Iteration000+Iteration001-009.out -eo ukbb_haplotype_translate2_sub.Iteration000+Iteration001-009.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr16.52100000-53100000.subset.txt fisher_exact.Iteration000+Iteration001-009.chr16.52100000-53100000.patterns_0000001-2521373.Results.txt 50"

# precompute haplotype counts
bsub -P SJLIFE -J ukbb_haplotype_model9_sub.ukbb_discovery.Iteration000+Iteration001-009 -oo ukbb_haplotype_model9_sub.v2.ukbb_discovery.Iteration000+Iteration001-009.out -eo ukbb_haplotype_model9_sub.v2.ukbb_discovery.Iteration000+Iteration001-009.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52100000-53100000.txt fisher_exact.Iteration000+Iteration001-009.chr16.52100000-53100000.patterns_0000001-2521373.Results.translated.txt \"1\" \"rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0\" \"1\" 16 52100000,53100000 50 \"ukbb_discovery.Iteration000+Iteration001-009\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2 /home/wletsou/scripts"

bsub -P SJLIFE -J ukbb_haplotype_model9_sub.drive_replication.Iteration000+Iteration001-009 -oo ukbb_haplotype_model9_sub.v2.drive_replication.Iteration000+Iteration001-009.out -eo ukbb_haplotype_model9_sub.v2.drive_replication.Iteration000+Iteration001-009.err -R "rusage[mem=256]" -w "done(179707180)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52100000-53100000.txt fisher_exact.Iteration000+Iteration001-009.chr16.52100000-53100000.patterns_0000001-2521373.Results.translated.txt \"1\" \"rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0\" \"1\" 16 52100000,53100000 50 \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2 /home/wletsou/scripts"

# forward selection
bsub -P SJLIFE -J ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-009 -oo ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-009.out -eo ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-009.err -R "rusage[mem=256]" -w "done(179707185)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.v2.sh ukbb_discovery.Iteration000+Iteration001-009.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt \"1\" \"1\" 1e-5 ukbb_discovery.Iteration000+Iteration001-009.Significant_patterns.txt \"ukbb_discovery.Iteration000+Iteration001-009\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2 /home/wletsou/scripts"

bsub -P SJLIFE -J dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-009 -oo dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-009.out -eo dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-009.err -w "done(179707186)" -R "rusage[mem=256]" "sh /home/wletsou/scripts/dbgap28544_haplotype_model9_iterate.v2.sh drive_replication.Iteration000+Iteration001-009.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt \"1\" \"1\" 1e-5 drive_replication.Iteration000+Iteration001-009.Significant_patterns.txt \"drive_replication.Iteration000+Iteration001-009\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2 /home/wletsou/scripts"

awk -v alpha=0.05 'BEGIN{OFS = "\t"} NR==FNR && $2 > 1 && $5 < 1e-5{seen[$1]=$0; n+=1; next} ($1 in seen && $2 > 1 && $5 < alpha){m+=1} END{printf("%s of %s replicated at p < %0.2e\n",m + 0,n + 0,alpha)}' ukbb_discovery.Conditional_haplotype_effects.h2.txt drive_replication.Conditional_haplotype_effects.h2.txt
# 1 of 123 replicated at p < 5.00e-02

awk 'BEGIN{OFS = "\t"} NR==FNR && $2 > 1 && $5 < 1e-5{seen[$1]=$0; next} ($1 in seen && $2 > 1 && $5 < 0.05){a = $1; $1 = ""; print seen[a],$0}' ukbb_discovery.Conditional_haplotype_effects.h2.txt drive_replication.Conditional_haplotype_effects.h2.txt | column -t
# rs117213474_A=0,rs62042828_T=0,rs56359341_G=0,rs117328487_C=0,rs111279855_C=0,rs1420281_C=0,rs60526343_G=0,rs72792181_C=0,rs116988362_C=0,rs4785118_A=0,rs17268400_C=0,rs117108822_A=0,rs72792197_C=0,rs117285516_A=0,rs116048640_A=0,rs113083456_T=0,rs1345316_A=1,rs2287071_G=0,rs112797109_T=0,rs72794126_T=0,rs74393050_A=0,rs17270950_C=0,rs117832727_C=0,rs113443627_G=0,rs117226581_T=0,rs116841226_A=0,rs117691720_C=0,rs34042425_T=0,rs194392_G=0,rs17535828_C=0,rs1362554_T=0,rs79687889_G=0,rs75191414_A=0,rs117773755_A=0,rs80072244_C=0,rs80179772_C=0,rs146046759_A=0,rs3743797_T=0,rs13333858_G=0,rs9935437_G=0,rs74017815_C=0,rs16951186_T=0,rs28617640_A=0,rs62043288_T=0,rs117175422_A=0,rs12324961_C=0,rs12927162_G=0,rs118136222_T=0,rs118155214_C=0,rs77700051_C=0,rs117106695_C=0,rs118186521_C=0,rs16951525_C=0,rs79501461_T=0,rs77362202_T=0,rs117155117_G=0,rs62043164_G=0,rs17361619_C=0,rs77960456_G=0,rs117544975_C=0,rs113555472_A=0,rs59793894_G=0,rs74020833_C=0,rs116849383_T=0,rs17296337_G=0,rs4377151_G=0,rs117947188_C=0,rs61131526_A=0,rs117553291_G=0,rs17296940_C=0,rs76028093_C=0,rs79291013_G=0,rs8044632_C=0,rs9939896_T=0,rs79993873_C=0,rs78320298_G=0,rs62041381_C=0,rs117222080_T=0,rs117551222_C=0,rs7500472_A=0,rs12446016_G=0,rs150704739_T=0,rs142367768_T=0,rs4783791_G=0,rs117911537_T=0,rs1344484_C=0,rs73597288_A=0,rs78322982_A=0,rs62043282_A=0,rs62043323_A=0,rs80356178_A=0,rs34261756_C=0,rs117400881_C=0,rs4784276_C=0,rs78724580_A=0,rs72812119_G=0,rs76365437_A=0,rs13335861_G=0,rs36099013_A=0,rs72812136_T=0,rs34382563_T=0,rs12598778_C=0,rs12922187_G=0,rs79088862_A=0,rs4784287_T=1,rs117297253_A=0,rs3931698_T=1,rs62049054_C=0,rs113465432_T=0  2.023716e+00  1.670281e+00  2.451939e+00  6.728349e-06  3.076770e-03  1114  5.826212e-03  2.932747e-03  105  1009  1.561212e+00  1.221255e+00  1.995802e+00  4.465622e-02  2.547610e-03  282  3.010245e-03  1.997469e-03  181  101

mkdir permutation1
cp dbgap28544_cases.indiv permutation1/dbgap28544_cases.indiv
cp dbgap28544_controls.indiv permutation1/dbgap28544_controls.indiv
cd permutation1

# permute UKBB case-and control-h1-carriers
bsub -P SJLIFE -J ukbb_haplotype_permute.v2 -oo ukbb_haplotype_permute.v2.out -eo ukbb_haplotype_permute.v2.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_permute.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52100000-53100000.txt rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0 20200116 ukbb_bca_cases,ukbb_bca_controls ukb.bca.pheno.txt /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/permutation1 /home/wletsou/scripts"

# extract haplotypes of permuted cases and controls
bsub -P SJLIFE -J ukbb_bca_overlap.chr16.52100000-53100000 -oo ukbb_bca_overlap.chr16.52100000-53100000.out -eo ukbb_bca_overlap.chr16.52100000-53100000.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_submission.v2.sh ukbb_bca_cases,ukbb_bca_controls dbgap28544_cases,dbgap28544_controls \"\" 16 52100000,53100000 rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0 1 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/ukbb_snp_list.chr16.52100000-53100000.txt /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/ukbb.topmed.hg19_chr16.52100000-53100000.hg38.vcf.gz,/scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/dbgap28544_cases+dbgap28544_controls.hg19_chr16.52100000-53100000.hg38.vcf.gz /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr16.qced.anno.info \"\" \"\" \"\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/permutation1 /home/wletsou/scripts"

awk 'NR==1{char=length(NF-1)} NR>1{printf "1\t"; for (i=2;i<=NF;i++) {printf "%s%0"char"d_%s",(i>2?",":""),i-1,$i}; printf "\n"}' haplotype_estimates.ukbb_bca_cases.chr16.52100000-53100000.subset.txt | awk 'BEGIN{OFS="\t"} {seen[$2]+=$1} END{for (i in seen) {print seen[i],i} }' > Pattern_combined.Iteration000.chr16.52100000-53100000_2,j.txt # starting with all chromosomes, not just pairs

# First round of overlaps:
bsub -P SJLIFE -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr16.52100000-53100000 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr16.52100000-53100000.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr16.52100000-53100000.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ChromosomeOverlap_iteration_sub_parallel.v3.sh chr16.52100000-53100000 2 2,j 50 0 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/permutation1 /home/wletsou/scripts"

# p-values for first round of overlaps
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr16.52100000-53100000.subset.txt,haplotype_estimates.ukbb_bca_controls.chr16.52100000-53100000.subset.txt Closed_patterns_stats.chr16.52100000-53100000_2,j.txt 50 \"1\" \"Iteration001.chr16.52100000-53100000\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/permutation1 /home/wletsou/scripts"

# remove older results if staring again
rm -r *Iteration00[2-9]*

awk '($1<2){print $0}' Closed_patterns_stats.chr16.52100000-53100000_2,j.txt > Closed_patterns_stats.chr16.52100000-53100000_2,j.tmp
mv Closed_patterns_stats.chr16.52100000-53100000_2,j.tmp Closed_patterns_stats.chr16.52100000-53100000_2,j.txt

test -f Pattern_combined_old.Iteration001.chr16.52100000-53100000_2,j.txt && mv Pattern_combined_old.Iteration001.chr16.52100000-53100000_2,j.txt Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.txt
awk 'NR==FNR{seen[$1]; next} ($2 in seen){print $0}' <(awk '(NR>1 && $5+0<8e-6 && $4+0>1){print $0}' fisher_exact.Iteration001.chr16.52100000-53100000.patterns_0000001-5840900.txt | awk 'BEGIN{array[0]} {delete array; n=split($1,array,","); if (n<snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]+0 || snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]+0==0) {snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]=n; seen[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]=$1} } END{for (i in seen) {print seen[i]} }') Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.txt > Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.tmp # take shortest pattern of a family with the same p-value
mv Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.txt Pattern_combined_old.Iteration001.chr16.52100000-53100000_2,j.txt
mv Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.tmp Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.txt

# Begin iterations:
bsub -P SJLIFE -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr16.52100000-53100000 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr16.52100000-53100000.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr16.52100000-53100000.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ChromosomeOverlap_iteration_sub_parallel.v3.sh chr16.52100000-53100000 2 2,j 50 1 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/permutation1 /home/wletsou/scripts"

# (filtered) patterns appearing for the first time at Iterations 001-006
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1){print $0}' Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.txt Closed_patterns_stats.chr16.52100000-53100000_2,j.txt > Closed_patterns_stats.Iteration001-006.chr16.52100000-53100000_2,j.txt

# p-values for Iterations 001-006
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-006 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-006.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-006.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr16.52100000-53100000.subset.txt,haplotype_estimates.ukbb_bca_controls.chr16.52100000-53100000.subset.txt Closed_patterns_stats.Iteration001-006.chr16.52100000-53100000_2,j.txt 50 \"\" \"Iteration001-006.chr16.52100000-53100000\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/permutation1 /home/wletsou/scripts"

# top closed patterns
awk '($5+0<1e-4 && $4+0>1){print $0}' fisher_exact.Iteration001-006.chr16.52100000-53100000.patterns_0000001-1569485.txt > fisher_exact.Iteration001-006.chr16.52100000-53100000.patterns_0000001-1569485.Results.txt

# tranlate haplotype to rsid
bsub -P SJLIFE -J ukbb_haplotype_translate2_sub -oo ukbb_haplotype_translate2_sub.out -eo ukbb_haplotype_translate2_sub.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr16.52100000-53100000.subset.txt fisher_exact.Iteration001-006.chr16.52100000-53100000.patterns_0000001-1569485.Results.txt 50"

# precompute haplotype counts
bsub -P SJLIFE -J ukbb_haplotype_model9_sub.v2.ukbb_discovery -oo ukbb_haplotype_model9_sub.v2.ukbb_discovery.out -eo ukbb_haplotype_model9_sub.v2.ukbb_discovery.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/permutation1/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52100000-53100000.txt fisher_exact.Iteration001-006.chr16.52100000-53100000.patterns_0000001-1569485.Results.translated.txt \"1\" \"rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0\" \"1\" 16 52100000,53100000 50 \"ukbb_discovery\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/permutation1 /home/wletsou/scripts"

bsub -P SJLIFE -J ukbb_haplotype_model9_sub.v2.drive_replication -oo ukbb_haplotype_model9_sub.v2.drive_replication.out -eo ukbb_haplotype_model9_sub.v2.drive_replication.err -R "rusage[mem=256]" -w "done(180226811)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52100000-53100000.txt fisher_exact.Iteration001-006.chr16.52100000-53100000.patterns_0000001-1569485.Results.translated.txt \"1\" \"rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0\" \"1\" 16 52100000,53100000 50 \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/permutation1 /home/wletsou/scripts"

# forward selection
bsub -P SJLIFE -J ukbb_haplotype_model9_iterate.v2.discovery -oo ukbb_haplotype_model9_iterate.v2.discovery.out -eo ukbb_haplotype_model9_iterate.v2.discovery.err -R "rusage[mem=256]" -w "done(180226935)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.v2.sh ukbb_discovery.allele_counts.txt /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/permutation1/ukb.bca.pheno.txt \"1\" \"1\" 1e-5 ukbb_discovery.Significant_patterns.txt \"ukbb_discovery\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/permutation1 /home/wletsou/scripts"

bsub -P SJLIFE -J dbgap28544_haplotype_model9_iterate.v2.replication -oo dbgap28544_haplotype_model9_iterate.v2.replication.out -eo dbgap28544_haplotype_model9_iterate.v2.replication.err -w "done(180227127)" -R "rusage[mem=256]" "sh /home/wletsou/scripts/dbgap28544_haplotype_model9_iterate.v2.sh drive_replication.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt \"1\" \"1\" 1e-5 drive_replication.Significant_patterns.txt \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/permutation1 /home/wletsou/scripts"

# replication
awk -v alpha=0.05 'BEGIN{OFS = "\t"} NR==FNR && $2 > 1 && $5 < 1e-5{seen[$1]=$0; n+=1; next} ($1 in seen && $2 > 1 && $5 < alpha){m+=1} END{printf("%s of %s replicated at p < %0.2e\n",m + 0,n + 0,alpha)}' ukbb_discovery.Conditional_haplotype_effects.h2.txt drive_replication.Conditional_haplotype_effects.h2.txt
# 6 of 60 replicated at p < 5.00e-02
awk -v alpha=0.05 'BEGIN{OFS = "\t"} NR==FNR && $2 > 1 && $5 < 1e-5{seen[$1]=$0; next} ($1 in seen && $2 > 1 && $5 < alpha){a = $1; $1 = ""; print seen[a],$0}' ukbb_discovery.Conditional_haplotype_effects.h2.txt drive_replication.Conditional_haplotype_effects.h2.txt | column -t

# (filtered) patterns appearing for the first time at Iterations 000,001-006
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1 || $1==0){print $0}' Pattern_combined.Iteration001.chr16.52100000-53100000_2,j.txt Closed_patterns_stats.chr16.52100000-53100000_2,j.txt > Closed_patterns_stats.Iteration000+Iteration001-006.chr16.52100000-53100000_2,j.txt

# p-values for Iterations 000,001-006
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-006 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-006.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-006.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr16.52100000-53100000.subset.txt,haplotype_estimates.ukbb_bca_controls.chr16.52100000-53100000.subset.txt Closed_patterns_stats.Iteration000+Iteration001-006.chr16.52100000-53100000_2,j.txt 50 \"\" \"Iteration000+Iteration001-006.chr16.52100000-53100000\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.3/phase2/permutation1 /home/wletsou/scripts"

# top closed patterns
awk '($5+0<1e-4 && $4+0>1){print $0}' fisher_exact.Iteration000+Iteration001-006.chr16.52100000-53100000.patterns_0000001-1572904.txt > fisher_exact.Iteration000+Iteration001-006.chr16.52100000-53100000.patterns_0000001-1572904.Results.txt

# minimum p-value for contiguous patterns
awk 'BEGIN{min=1} NR==FNR{seen[$2]; next} ($1 in seen){ if ($5<min+0) {min=$5} } END{print min}' Pattern_combined.Iteration000.chr16.52100000-53100000_2,j.txt fisher_exact.Iteration000+Iteration001-006.chr16.52100000-53100000.patterns_0000001-1572904.txt
# 0.00302786258461715
