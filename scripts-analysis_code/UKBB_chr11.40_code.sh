mkdir UKBB_chr11.40
cd UKBB_chr11.40

mkdir /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/ukbb_genotyped_snps
cd ukbb_genotyped_snps

awk 'BEGIN{OFS="\t"} ($3==1){print $1,$1}' /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_cases.indiv
awk 'BEGIN{OFS="\t"} ($3==0){print $1,$1}' /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_controls.indiv

bsub -P SJLIFE -J ukbb_hybrid_haplotype3.chr11.68850000-69500000 -oo ukbb_hybrid_haplotype3.chr11.68850000-69231641.out -eo ukbb_hybrid_haplotype3.chr11.68850000-69500000.err -R "rusage[mem=20000]" "sh /home/wletsou/scripts/ukbb_hybrid_haplotype3.sh ukbb_bca_cases.indiv,ukbb_bca_controls.indiv 11 68850000,69500000 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/biobank/gds/ukb.bca.hap/ukb.bca.hap.chr11.new.vcf.gz"

cd /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40

awk 'BEGIN{OFS="\t"} ($3==1){print $1,$1}' /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_cases.indiv
awk 'BEGIN{OFS="\t"} ($3==0){print $1,$1}' /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_controls.indiv

# get genotyped SNPs in region
awk 'BEGIN{OFS="\t"} ($2>=69231642 && $2<=69431642){print $3}' /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/ukbb_genotyped_snps/ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69500000.haplotypes.vcf > ukbb_snp_list.chr11.69231642-69431642.txt

# check that SNPs have good "imputation" quality
awk 'NR==FNR{snp[$1]; next} ($8 in snp && $15==1){print $8}' ukbb_snp_list.chr11.69231642-69431642.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr11.qced.anno.info > ukbb_snp_list.chr11.69231642-69431642.tmp
test -f ukbb_snp_list.chr11.69231642-69431642.tmp && mv ukbb_snp_list.chr11.69231642-69431642.tmp ukbb_snp_list.chr11.69231642-69431642.txt

# vcf file for UKBB
bsub -P SJLIFE -J ukbb_topmed_merge.chr11.69231642-69431642 -oo ukbb_topmed_merge.chr11.69231642-69431642.out -eo ukbb_topmed_merge.chr11.69231642-69431642.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_topmed_merge.sh ukbb_snp_list.chr11.69231642-69431642.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr11.qced.anno.info 11 69231642,69431642 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/raw/chr11"

# vcf file for dbGaP
bsub -P SJLIFE -J dbgap_haplotype_merge.chr11.69231642-69431642 -oo dbgap_haplotype_merge.chr11.69231642-69431642.out -eo dbgap_haplotype_merge.chr11.69231642-69431642.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/dbgap_haplotype_merge.sh dbgap28544_cases,dbgap28544_controls /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/chr11.qced.info ukbb_snp_list.chr11.69231642-69431642.txt 11 69231642,69431642 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/vcfs/dbgap28544.topmedr2.cleaned.hg38.chr11.vcf.gz"

# extract haplotypes for cases and controls
bsub -P SJLIFE -J ukbb_bca_overlap.chr11.69231642-69431642 -oo ukbb_bca_overlap.chr11.69231642-69431642.out -eo ukbb_bca_overlap.chr11.69231642-69431642.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_submission.v2.sh ukbb_bca_cases,ukbb_bca_controls dbgap28544_cases,dbgap28544_controls \"\" 11 69231642,69431642 \"\" \"\" ukbb_snp_list.chr11.69231642-69431642.txt ukbb.topmed.hg19_chr11.69231642-69431642.hg38.vcf.gz,dbgap28544_cases+dbgap28544_controls.hg19_chr11.69231642-69431642.hg38.vcf.gz /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr11.qced.anno.info \"\" \"\" \"\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40 /home/wletsou/scripts"

awk 'NR==1{char=length(NF-1)} NR>1{printf "1\t"; for (i=2;i<=NF;i++) {printf "%s%0"char"d_%s",(i>2?",":""),i-1,$i}; printf "\n"}' haplotype_estimates.ukbb_bca_cases.chr11.69231642-69431642.txt | awk 'BEGIN{OFS="\t"} {seen[$2]+=$1} END{for (i in seen) {print seen[i],i} }' > Pattern_combined.Iteration000.chr11.69231642-69431642_2,j.txt # starting with all chromosomes, not just pairs

# First round of overlaps:
bsub -P SJLIFE -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr11.69231642-69431642 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr11.69231642-69431642.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr11.69231642-69431642.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ChromosomeOverlap_iteration_sub_parallel.v3.sh chr11.69231642-69431642 2 2,j 50 0 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40 /home/wletsou/scripts"

# p-values for first round of overlaps
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr11.69231642-69431642.txt,haplotype_estimates.ukbb_bca_controls.chr11.69231642-69431642.txt Closed_patterns_stats.chr11.69231642-69431642_2,j.txt 50 \"1\" \"Iteration001.chr11.69231642-69431642\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40 /home/wletsou/scripts"

# remove older results if staring again
rm -r *Iteration00[2-9]*
awk '($1<2){print $0}' Closed_patterns_stats.chr11.69231642-69431642_2,j.txt > Closed_patterns_stats.chr11.69231642-69431642_2,j.tmp
mv Closed_patterns_stats.chr11.69231642-69431642_2,j.tmp Closed_patterns_stats.chr11.69231642-69431642_2,j.txt

test -f Pattern_combined_old.Iteration001.chr11.69231642-69431642_2,j.txt && mv Pattern_combined_old.Iteration001.chr11.69231642-69431642_2,j.txt Pattern_combined.Iteration001.chr11.69231642-69431642_2,j.txt
awk 'NR==FNR{seen[$1]; next} ($2 in seen){print $0}' <(awk '(NR>1 && $5+0<1e-9 && $4+0>1){print $1}' fisher_exact.Iteration001.chr11.69231642-69431642.patterns_0000001-2962894.txt) Pattern_combined.Iteration001.chr11.69231642-69431642_2,j.txt > Pattern_combined.Iteration001.chr11.69231642-69431642_2,j.tmp
mv Pattern_combined.Iteration001.chr11.69231642-69431642_2,j.txt Pattern_combined_old.Iteration001.chr11.69231642-69431642_2,j.txt
mv Pattern_combined.Iteration001.chr11.69231642-69431642_2,j.tmp Pattern_combined.Iteration001.chr11.69231642-69431642_2,j.txt

# Begin iterations:
bsub -P SJLIFE -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr11.69231642-69431642 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr11.69231642-69431642.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr11.69231642-69431642.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ChromosomeOverlap_iteration_sub_parallel.v3.sh chr11.69231642-69431642 2 2,j 50 1 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40 /home/wletsou/scripts"

# (filtered) patterns appearing for the first time at Iterations 001-005
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1){print $0}' Pattern_combined.Iteration001.chr11.69231642-69431642_2,j.txt Closed_patterns_stats.chr11.69231642-69431642_2,j.txt > Closed_patterns_stats.Iteration001-005.chr11.69231642-69431642_2,j.txt

# p-values for Iterations 001-005
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr11.69231642-69431642.txt,haplotype_estimates.ukbb_bca_controls.chr11.69231642-69431642.txt Closed_patterns_stats.Iteration001-005.chr11.69231642-69431642_2,j.txt 50 \"\" \"Iteration001-005.chr11.69231642-69431642\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40 /home/wletsou/scripts"

# top closed patterns
awk '($5+0<1e-15 && $4+0>1){print $0}' fisher_exact.Iteration001-005.chr11.69231642-69431642.patterns_000001-585850.txt > fisher_exact.Iteration001-005.chr11.69231642-69431642.patterns_000001-585850.Results.txt

awk 'BEGIN{min=1} NR==FNR{seen[$2]; next} ($1 in seen){if ($5 < min ) {min=$5} } END{print min}' Pattern_combined.Iteration000.chr11.69231642-69431642_2,j.txt fisher_exact.Iteration000+Iteration001-005.chr11.69231642-69431642.patterns_000001-588370.txt
# 2.3775821003027e-05

# tranlate haplotype to rsid
bsub -P SJLIFE -J ukbb_haplotype_translate2_sub -oo ukbb_haplotype_translate2_sub.out -eo ukbb_haplotype_translate2_sub.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr11.69231642-69431642.txt fisher_exact.Iteration001-005.chr11.69231642-69431642.patterns_000001-585850.Results.txt 50"

# precompute haplotype counts
bsub -P SJLIFE -J ukbb_haplotype_model9_sub.ukbb_discovery -oo ukbb_haplotype_model9_sub.v2.ukbb_discovery.out -eo ukbb_haplotype_model9_sub.v2.ukbb_discovery.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.69231642-69431642.txt fisher_exact.Iteration001-005.chr11.69231642-69431642.patterns_000001-585850.Results.translated.txt \"\" \"\" \"\" 11 69231642,69431642 50 \"ukbb_discovery\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40 /home/wletsou/scripts"

bsub -P SJLIFE -J ukbb_haplotype_model9_sub.drive_replication -oo ukbb_haplotype_model9_sub.v2.drive_replication.out -eo ukbb_haplotype_model9_sub.v2.drive_replication.err -R "rusage[mem=256]" -w "done(179279647)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr11.69231642-69431642.txt fisher_exact.Iteration001-005.chr11.69231642-69431642.patterns_000001-585850.Results.translated.txt \"\" \"\" \"\" 11 69231642,69431642 50 \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40 /home/wletsou/scripts"

# forward selection
bsub -P SJLIFE -J ukbb_haplotype_model9_iterate.v2.discovery -oo ukbb_haplotype_model9_iterate.v2.discovery.out -eo ukbb_haplotype_model9_iterate.v2.discovery.err -R "rusage[mem=256]" -w "done(179279664)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.v2.sh ukbb_discovery.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt \"\" \"\" 1e-5 ukbb_discovery.Significant_patterns.txt \"ukbb_discovery\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40 /home/wletsou/scripts"

bsub -P SJLIFE -J dbgap28544_haplotype_model9_iterate.v2.replication -oo dbgap28544_haplotype_model9_iterate.v2.replication.out -eo dbgap28544_haplotype_model9_iterate.v2.replication.err -w "done(179287492)" -R "rusage[mem=256]" "sh /home/wletsou/scripts/dbgap28544_haplotype_model9_iterate.v2.sh drive_replication.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt \"\" \"\" 1e-5 drive_replication.Significant_patterns.txt \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40 /home/wletsou/scripts"

# Include contiguous haplotypes (Iteration 000)

awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1 || $1 == 0){print $0}' Pattern_combined.Iteration001.chr11.69231642-69431642_2,j.txt Closed_patterns_stats.chr11.69231642-69431642_2,j.txt > Closed_patterns_stats.Iteration000+Iteration001-005.chr11.69231642-69431642_2,j.txt

# p-values for Iterations 000,001-005
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr11.69231642-69431642.txt,haplotype_estimates.ukbb_bca_controls.chr11.69231642-69431642.txt Closed_patterns_stats.Iteration000+Iteration001-005.chr11.69231642-69431642_2,j.txt 50 \"\" \"Iteration000+Iteration001-005.chr11.69231642-69431642\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40 /home/wletsou/scripts"

# top closed patterns
awk '($5+0<1e-15 && $4+0>1){print $0}' fisher_exact.Iteration000+Iteration001-005.chr11.69231642-69431642.patterns_000001-588370.txt > fisher_exact.Iteration000+Iteration001-005.chr11.69231642-69431642.patterns_000001-588370.Results.txt

awk 'BEGIN{min=1} NR==FNR{seen[$2]; next} ($1 in seen){ if ($5<min+0) {min=$5} } END{print min}' Pattern_combined.Iteration000.chr11.69231642-69431642_2,j.txt fisher_exact.Iteration000+Iteration001-005.chr11.69231642-69431642.patterns_000001-588370.txt
# 2.3775821003027e-05

awk 'NR==FNR && $1==0{seen[$3]; next} ($5+0<1e-15 && $4+0>1 || $1 in seen){print $0}' Closed_patterns_stats.Iteration000+Iteration001-005.chr11.69231642-69431642_2,j.txt fisher_exact.Iteration000+Iteration001-005.chr11.69231642-69431642.patterns_000001-588370.txt > fisher_exact.Iteration000+Iteration001-005.chr11.69231642-69431642.patterns_000001-588370.Results.txt # force to include Iteration000

# tranlate haplotype to rsid
bsub -P SJLIFE -J ukbb_haplotype_translate2_sub.Iteration000+Iteration001-005 -oo ukbb_haplotype_translate2_sub.Iteration000+Iteration001-005.out -eo ukbb_haplotype_translate2_sub.Iteration000+Iteration001-005.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr11.69231642-69431642.txt fisher_exact.Iteration000+Iteration001-005.chr11.69231642-69431642.patterns_000001-588370.Results.txt 50"

# precompute haplotype counts
bsub -P SJLIFE -J ukbb_haplotype_model9_sub.ukbb_discovery.Iteration000+Iteration001-005 -oo ukbb_haplotype_model9_sub.v2.ukbb_discovery.Iteration000+Iteration001-005.out -eo ukbb_haplotype_model9_sub.v2.ukbb_discovery.Iteration000+Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.69231642-69431642.txt fisher_exact.Iteration000+Iteration001-005.chr11.69231642-69431642.patterns_000001-588370.Results.translated.txt \"\" \"\" \"\" 11 69231642,69431642 50 \"ukbb_discovery.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40 /home/wletsou/scripts"

bsub -P SJLIFE -J ukbb_haplotype_model9_sub.drive_replication.Iteration000+Iteration001-005 -oo ukbb_haplotype_model9_sub.v2.drive_replication.Iteration000+Iteration001-005.out -eo ukbb_haplotype_model9_sub.v2.drive_replication.Iteration000+Iteration001-005.err -R "rusage[mem=256]" -w "done(179279647)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr11.69231642-69431642.txt fisher_exact.Iteration000+Iteration001-005.chr11.69231642-69431642.patterns_000001-588370.Results.translated.txt \"\" \"\" \"\" 11 69231642,69431642 50 \"drive_replication.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40 /home/wletsou/scripts"

# forward selection
bsub -P SJLIFE -J ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-005 -oo ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-005.out -eo ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-005.err -R "rusage[mem=256]" -w "done(179279664)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.v2.sh ukbb_discovery.Iteration000+Iteration001-005.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt \"\" \"\" 1e-5 ukbb_discovery.Iteration000+Iteration001-005.Significant_patterns.txt \"ukbb_discovery.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40 /home/wletsou/scripts"

bsub -P SJLIFE -J dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-005 -oo dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-005.out -eo dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-005.err -w "done(179287492)" -R "rusage[mem=256]" "sh /home/wletsou/scripts/dbgap28544_haplotype_model9_iterate.v2.sh drive_replication.Iteration000+Iteration001-005.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt \"\" \"\" 1e-5 drive_replication.Iteration000+Iteration001-005.Significant_patterns.txt \"drive_replication.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40 /home/wletsou/scripts"

# reduction
sh /home/wletsou/scripts/haplotype_rpart.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.69231642-69431642.txt "rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs2298764_C=0,rs117222887_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "" 1 # h1
#  reduced = rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0

# discovery
sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.69231642-69431642.txt "rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "" "" "" "" h1.discovery
# LRT p value = 1.542142e-21 on 1 degree(s) of freedom with HR = 1.255279e+00 (1.199565e+00 to 1.313580e+00), frequency = 9.579968e-02 (34686) = 1.165242e-01/9.471408e-02 (2100/32586), p = 3.960361e-21

# DRIVE replication
sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr11.69231642-69431642.txt "rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "" "" "" "" h1.replication
# LRT p value = 3.109647e-21 on 1 degree(s) of freedom with OR = 1.205913e+00 (1.159833e+00 to 1.253823e+00), frequency = 1.058703e-01 (11719) = 1.135078e-01/9.678823e-02 (6825/4894), p = 1.793186e-19

sh /home/wletsou/scripts/vcf2bed.v3.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.69231642-69431642.haplotypes.hg38.vcf "rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "" "h1.bed" "1" # h1, with alleles and header

awk 'BEGIN{OFS="\t"} ("rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" ~ $3){print $1,$2,$3,$4,$5}' ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.69231642-69431642.haplotypes.hg38.vcf
# h1 rows of vcf file

# chr11   69424973        rs7105934       G       A
# chr11   69432853        rs79373485      C       T
# chr11   69458108        rs4444099       C       T
# chr11   69475035        rs117752342     C       T
# chr11   69485307        rs1122316       G       A
# chr11   69508496        rs76809977      C       T
# chr11   69517492        rs559664        A       G
# chr11   69517944        rs657315        C       T
# chr11   69524849        rs498931        A       G
# chr11   69528387        rs71465432      C       T
# chr11   69539114        rs183782062     G       T
# chr11   69562741        rs79442425      C       T
# chr11   69564519        rs111929748     G       T
# chr11   69572513        rs117490805     G       A
# chr11   69581315        rs79241527      C       T
# chr11   69603668        rs72932500      A       G
# chr11   69609602        rs57162717      A       T

mkdir phase2
cp ukbb_bca_cases.indiv phase2/ukbb_bca_cases.indiv
cp ukbb_bca_controls.indiv phase2/ukbb_bca_controls.indiv
cp dbgap28544_cases.indiv phase2/dbgap28544_cases.indiv
cp dbgap28544_controls.indiv phase2/dbgap28544_controls.indiv
cd phase2

# get h1 and list of SNPs in upstream region
HAPLOTYPE_SNPS=$(echo "rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" | sed 's/_[A-Z]=[0-9]//g') # snps in h1, no alleles

awk 'BEGIN{OFS="\t"} ( ($2>=68850000 && $2<69231642) || "'$HAPLOTYPE_SNPS'" ~ $3){print $3}' /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/ukbb_genotyped_snps/ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69500000.haplotypes.vcf > ukbb_snp_list.chr11.68850000-69231641.txt

# check that SNPs have good "imputation" quality
awk 'NR==FNR{snp[$1]; next} ($8 in snp && $15==1){print $8}' ukbb_snp_list.chr11.68850000-69231641.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr11.qced.anno.info > ukbb_snp_list.chr11.68850000-69231641.tmp
test -f ukbb_snp_list.chr11.68850000-69231641.tmp && mv ukbb_snp_list.chr11.68850000-69231641.tmp ukbb_snp_list.chr11.68850000-69231641.txt

# vcf files of region (UKBB and dbGaP)
bsub -P SJLIFE -J ukbb_topmed_merge.chr11.68850000-69231641 -oo ukbb_topmed_merge.chr11.68850000-69231641.out -eo ukbb_topmed_merge.chr11.68850000-69231641.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_topmed_merge.sh ukbb_snp_list.chr11.68850000-69231641.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr11.qced.anno.info 11 68850000,69231641 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/raw/chr11"

bsub -P SJLIFE -J dbgap_haplotype_merge.chr11.68850000-69231641 -oo dbgap_haplotype_merge.chr11.68850000-69231641.out -eo dbgap_haplotype_merge.chr11.68850000-69231641.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/dbgap_haplotype_merge.sh dbgap28544_cases,dbgap28544_controls \"/research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt\" \"/research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/chr11.qced.info\" ukbb_snp_list.chr11.68850000-69231641.txt 11 68850000,69231641 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/vcfs/dbgap28544.topmedr2.cleaned.hg38.chr11.vcf.gz"

# extract haplotypes for cases and controls
bsub -P SJLIFE -J ukbb_bca_overlap.chr11.68850000-69231641 -oo ukbb_bca_overlap.chr11.68850000-69231641.out -eo ukbb_bca_overlap.chr11.68850000-69231641.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_submission.v2.sh ukbb_bca_cases,ukbb_bca_controls dbgap28544_cases,dbgap28544_controls \"\" 11 68850000,69231641 rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0 1 ukbb_snp_list.chr11.68850000-69231641.txt ukbb.topmed.hg19_chr11.68850000-69231641.hg38.vcf.gz,dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.hg38.vcf.gz /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr11.qced.anno.info \"\" \"\" \"\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

awk 'NR==1{char=length(NF-1)} NR>1{printf "1\t"; for (i=2;i<=NF;i++) {printf "%s%0"char"d_%s",(i>2?",":""),i-1,$i}; printf "\n"}' haplotype_estimates.ukbb_bca_cases.chr11.68850000-69231641.subset.txt | awk 'BEGIN{OFS="\t"} {seen[$2]+=$1} END{for (i in seen) {print seen[i],i} }' > Pattern_combined.Iteration000.chr11.68850000-69231641_2,j.txt # starting with all chromosomes, not just pairs

# First round of overlaps:
bsub -P SJLIFE -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr11.68850000-69231641 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr11.68850000-69231641.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr11.68850000-69231641.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ChromosomeOverlap_iteration_sub_parallel.v3.sh chr11.68850000-69231641 2 2,j 50 0 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

# p-values for first round of overlaps
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr11.68850000-69231641.subset.txt,haplotype_estimates.ukbb_bca_controls.chr11.68850000-69231641.subset.txt Closed_patterns_stats.chr11.68850000-69231641_2,j.txt 50 \"1\" \"Iteration001.chr11.68850000-69231641\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

# remove older results if staring again
rm -r *Iteration00[2-9]*

awk '($1<2){print $0}' Closed_patterns_stats.chr11.68850000-69231641_2,j.txt > Closed_patterns_stats.chr11.68850000-69231641_2,j.tmp
mv Closed_patterns_stats.chr11.68850000-69231641_2,j.tmp Closed_patterns_stats.chr11.68850000-69231641_2,j.txt

test -f Pattern_combined_old.Iteration001.chr11.68850000-69231641_2,j.txt && mv Pattern_combined_old.Iteration001.chr11.68850000-69231641_2,j.txt Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.txt
awk 'NR==FNR{seen[$1]; next} ($2 in seen){print $0}' <(awk '(NR>1 && $5+0<1e-4 && $4+0>1){print $0}' fisher_exact.Iteration001.chr11.68850000-69231641.patterns_000001-241811.txt | awk 'BEGIN{array[0]} {delete array; n=split($1,array,","); if (n<snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]+0 || snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]+0==0) {snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]=n; seen[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]=$1} } END{for (i in seen) {print seen[i]} }') Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.txt > Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.tmp # take shortest pattern of a family with the same p-value
mv Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.txt Pattern_combined_old.Iteration001.chr11.68850000-69231641_2,j.txt
mv Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.tmp Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.txt

# Begin iterations:
bsub -P SJLIFE -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr11.68850000-69231641 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr11.68850000-69231641.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr11.68850000-69231641.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ChromosomeOverlap_iteration_sub_parallel.v3.sh chr11.68850000-69231641 2 2,j 50 1 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

# (filtered) patterns appearing for the first time at Iterations 001-005
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1){print $0}' Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.txt Closed_patterns_stats.chr11.68850000-69231641_2,j.txt > Closed_patterns_stats.Iteration001-005.chr11.68850000-69231641_2,j.txt

# p-values for Iterations 001-005
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr11.68850000-69231641.subset.txt,haplotype_estimates.ukbb_bca_controls.chr11.68850000-69231641.subset.txt Closed_patterns_stats.Iteration001-005.chr11.68850000-69231641_2,j.txt 50 \"\" \"Iteration001-005.chr11.68850000-69231641\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

# top closed patterns
awk '($5+0<1e-4 && $4+0>1){print $0}' fisher_exact.Iteration001-005.chr11.68850000-69231641.patterns_000001-167998.txt > fisher_exact.Iteration001-005.chr11.68850000-69231641.patterns_000001-167998.Results.txt

# tranlate haplotype to rsid
bsub -P SJLIFE -J ukbb_haplotype_translate2_sub -oo ukbb_haplotype_translate2_sub.out -eo ukbb_haplotype_translate2_sub.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr11.68850000-69231641.subset.txt fisher_exact.Iteration001-005.chr11.68850000-69231641.patterns_000001-167998.Results.txt 50"

# precompute haplotype counts
bsub -P SJLIFE -J ukbb_haplotype_model9_sub.v2.ukbb_discovery -oo ukbb_haplotype_model9_sub.v2.ukbb_discovery.out -eo ukbb_haplotype_model9_sub.v2.ukbb_discovery.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69231641.txt fisher_exact.Iteration001-005.chr11.68850000-69231641.patterns_000001-167998.Results.translated.txt \"1\" \"rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0\" \"1\" 11 68850000,69231641 50 \"ukbb_discovery\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

bsub -P SJLIFE -J ukbb_haplotype_model9_sub.v2.drive_replication -oo ukbb_haplotype_model9_sub.v2.drive_replication.out -eo ukbb_haplotype_model9_sub.v2.drive_replication.err -R "rusage[mem=256]" -w "done(179326475)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr11.68850000-69231641.txt fisher_exact.Iteration001-005.chr11.68850000-69231641.patterns_000001-167998.Results.translated.txt \"1\" \"rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0\" \"1\" 11 68850000,69231641 50 \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

# forward selection
bsub -P SJLIFE -J ukbb_haplotype_model9_iterate.v2.discovery -oo ukbb_haplotype_model9_iterate.v2.discovery.out -eo ukbb_haplotype_model9_iterate.v2.discovery.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.v2.sh ukbb_discovery.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt \"1\" \"1\" 1e-5 ukbb_discovery.Significant_patterns.txt \"ukbb_discovery\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

bsub -P SJLIFE -J dbgap28544_haplotype_model9_iterate.v2.replication -oo dbgap28544_haplotype_model9_iterate.v2.replication.out -eo dbgap28544_haplotype_model9_iterate.v2.replication.err -w "done(179327022)" -R "rusage[mem=256]" "sh /home/wletsou/scripts/dbgap28544_haplotype_model9_iterate.v2.sh drive_replication.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt \"1\" \"1\" 1e-5 drive_replication.Significant_patterns.txt \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

# Including contiguous haplotypes (Iteration000)

# (filtered) patterns appearing for the first time at Iterations 000,001-005
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1 || $1 == 0){print $0}' Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.txt Closed_patterns_stats.chr11.68850000-69231641_2,j.txt > Closed_patterns_stats.Iteration000+Iteration001-005.chr11.68850000-69231641_2,j.txt

# p-values for Iterations 000,001-005
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr11.68850000-69231641.subset.txt,haplotype_estimates.ukbb_bca_controls.chr11.68850000-69231641.subset.txt Closed_patterns_stats.Iteration000+Iteration001-005.chr11.68850000-69231641_2,j.txt 50 \"\" \"Iteration000+Iteration001-005.chr11.68850000-69231641\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

# top closed patterns
awk '($5+0<1e-4 && $4+0>1){print $0}' fisher_exact.Iteration000+Iteration001-005.chr11.68850000-69231641.patterns_000001-168697.txt > fisher_exact.Iteration000+Iteration001-005.chr11.68850000-69231641.patterns_000001-168697.Results.txt

# minimum p-value for contiguous patterns
awk 'BEGIN{min=1} NR==FNR{seen[$2]; next} ($1 in seen){ if ($5<min+0) {min=$5} } END{print min}' Pattern_combined.Iteration000.chr11.68850000-69231641_2,j.txt fisher_exact.Iteration000+Iteration001-005.chr11.68850000-69231641.patterns_000001-168697.txt
# 0.00313117114689362

awk 'NR==FNR && $1==0{seen[$3]; next} ($5+0<1e-4 && $4+0>1 || $1 in seen){print $0}' Closed_patterns_stats.Iteration000+Iteration001-005.chr11.68850000-69231641_2,j.txt fisher_exact.Iteration000+Iteration001-005.chr11.68850000-69231641.patterns_000001-168697.txt > fisher_exact.Iteration000+Iteration001-005.chr11.68850000-69231641.patterns_000001-168697.Results.txt # force to include Iteration000

# tranlate haplotype to rsid
bsub -P SJLIFE -J ukbb_haplotype_translate2_sub -oo ukbb_haplotype_translate2_sub.out -eo ukbb_haplotype_translate2_sub.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr11.68850000-69231641.subset.txt fisher_exact.Iteration000+Iteration001-005.chr11.68850000-69231641.patterns_000001-168697.Results.txt 50"

# precompute haplotype counts
bsub -P SJLIFE -J ukbb_haplotype_model9_sub.v2.ukbb_discovery.Iteration000+Iteration001-005 -oo ukbb_haplotype_model9_sub.v2.ukbb_discovery.Iteration000+Iteration001-005.out -eo ukbb_haplotype_model9_sub.v2.ukbb_discovery.Iteration000+Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69231641.txt fisher_exact.Iteration000+Iteration001-005.chr11.68850000-69231641.patterns_000001-168697.Results.translated.txt \"1\" \"rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0\" \"1\" 11 68850000,69231641 50 \"ukbb_discovery.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

bsub -P SJLIFE -J ukbb_haplotype_model9_sub.v2.drive_replication.Iteration000+Iteration001-005 -oo ukbb_haplotype_model9_sub.v2.drive_replication.Iteration000+Iteration001-005.out -eo ukbb_haplotype_model9_sub.v2.drive_replication.Iteration000+Iteration001-005.err -R "rusage[mem=256]" -w "done(179326475)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr11.68850000-69231641.txt fisher_exact.Iteration000+Iteration001-005.chr11.68850000-69231641.patterns_000001-168697.Results.translated.txt \"1\" \"rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0\" \"1\" 11 68850000,69231641 50 \"drive_replication.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

# forward selection
bsub -P SJLIFE -J ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-005 -oo ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-005.out -eo ukbb_haplotype_model9_iterate.v2.discovery.Iteration000+Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.v2.sh ukbb_discovery.Iteration000+Iteration001-005.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt \"1\" \"1\" 1e-5 ukbb_discovery.Iteration000+Iteration001-005.Significant_patterns.txt \"ukbb_discovery.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

bsub -P SJLIFE -J dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-005 -oo dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-005.out -eo dbgap28544_haplotype_model9_iterate.v2.replication.Iteration000+Iteration001-005.err -w "done(179327022)" -R "rusage[mem=256]" "sh /home/wletsou/scripts/dbgap28544_haplotype_model9_iterate.v2.sh drive_replication.Iteration000+Iteration001-005.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt \"1\" \"1\" 1e-5 drive_replication.Iteration000+Iteration001-005.Significant_patterns.txt \"drive_replication.Iteration000+Iteration001-005\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

# replication
awk -v alpha=0.05 'BEGIN{OFS = "\t"} NR==FNR && $2 > 1 && $5 < 1e-5{seen[$1]=$0; n+=1; next} ($1 in seen && $2 > 1 && $5 < alpha){m+=1} END{printf("%s of %s replicated at p < %0.2e\n",m + 0,n + 0,alpha)}' ukbb_discovery.Conditional_haplotype_effects.h2.txt drive_replication.Conditional_haplotype_effects.h2.txt
# 14 of 106 replicated at p < 5.00e-02
awk 'BEGIN{OFS = "\t"} NR==FNR && $2 > 1 && $5 < 1e-5{seen[$1]=$0; next} ($1 in seen && $2 > 1 && $5 < 0.05){a = $1; $1 = ""; print seen[a],$0}' ukbb_discovery.Conditional_haplotype_effects.h2.txt drive_replication.Conditional_haplotype_effects.h2.txt | column -t

awk 'BEGIN{OFS = "\t"} NR==FNR && $2 > 1 && $5 < 1e-5{seen[$1]=$0; next} ($1 in seen && $2 > 1 && $5 < 0.05){a = $1; $1 = ""; print seen[a],$0}' ukbb_discovery.Conditional_haplotype_effects.h2.txt drive_replication.Conditional_haplotype_effects.h2.txt | cut -f1 > replicated_patterns.chr11.68850000-69231641.txt

# counts for replicated haplotypes
bsub -P SJLIFE -J ukbb_haplotype_model9_sub.v2.replicated_patterns -oo ukbb_haplotype_model9_sub.v2.replicated_patterns.out -eo ukbb_haplotype_model9_sub.v2.replicated_patterns.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69231641.txt replicated_patterns.chr11.68850000-69231641.txt \"1\" \"rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0\" \"1\" 11 68850000,69231641 50 \"replicated_patterns\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

# forward selection for replicated patterns
bsub -P SJLIFE -J ukbb_haplotype_model9_iterate.v2.replicated_patterns -oo ukbb_haplotype_model9_iterate.v2.replicated_patterns.out -eo ukbb_haplotype_model9_iterate.v2.replicated_patterns.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.v2.sh replicated_patterns.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt \"1\" \"1\" 1e-5 replicated_patterns.Significant_patterns.txt \"replicated_patterns\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts"

# final model
awk '($6<1e-5){printf "%s:",$2} END{printf "\n"}' replicated_patterns.Significant_patterns.txt

sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69231641.txt "rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs74471298_A=0:rs2376558_C=1,rs61881030_A=0,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs12577213_A=0,rs73520512_T=0,rs3892895_G=1,rs7116054_T=0,rs74897859_T=0,rs17308715_T=0,rs67039008_A=1,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs74551015_A=0,rs117821797_T=0,rs12275849_C=0,rs56103266_C=0,rs117131950_A=0,rs116926312_T=0,rs11539762_A=0,rs142581206_T=0,rs11228610_C=0,rs12274095_A=0,rs12224376_A=0,rs7124547_T=1,rs11603814_G=1,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs79487139_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs76855139_T=0,rs74471298_A=0" "1" "rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "1" "" h2+h3.ukbb_discovery /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts

Rscript /home/wletsou/scripts/ukbb_haplotype_model5.R file=/scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/h2+h3.ukbb_discovery.new_allele_counts.txt groups=1,2,3 verbose=1 n=2 n_covs=10
# [1] "X000" "X100" "X101" "X110"
# [1] "DT[,h0:=X000]"
# [1] "DT[,h1.3:=X100+X101]"
# [1] "DT[,h1:=X100]"
# [1] "DT[,h2:=X110]"
# [1] "DT[,h3:=X101]"

#    affected  h1.3    h1  h2 h3
# 1:        0 32470 32465 116  5
# 2:        1  2074  2068  26  6

#    affected       h1.3         h1           h2           h3
# 1:        0 0.09437691 0.09436238 0.0003371642 1.453294e-05
# 2:        1 0.11508157 0.11474864 0.0014426812 3.329264e-04

#            h1.3           h1           h2           h3
# 1: 2.665607e-19 9.573352e-19 8.900272e-09 5.648169e-06

#     h1.3    h1  h2 h3
# 1: 34544 34533 142 11

#          h1.3         h1           h2           h3
# 1: 0.09540749 0.09537711 0.0003921915 3.038103e-05

# Surv(age,affected) ~ h1.3 + h2 + pc01 + pc02 + pc03 + pc04 + pc05 + pc06 + pc07 + pc08 + pc09 + pc10
# Call:
# coxph(formula = eval(parse(text = Y)), data = X)

#   n= 181034, number of events= 9011

#            coef  exp(coef)   se(coef)      z Pr(>|z|)
# h1.3  0.2181018  1.2437136  0.0232964  9.362  < 2e-16 ***
# h2    1.4374647  4.2100086  0.1964154  7.318 2.51e-13 ***
# pc01 -0.0010036  0.9989969  0.0069136 -0.145   0.8846
# pc02 -0.0146110  0.9854952  0.0071569 -2.042   0.0412 *
# pc03  0.0074527  1.0074806  0.0069030  1.080   0.2803
# pc04  0.0036003  1.0036068  0.0051789  0.695   0.4869
# pc05 -0.0006872  0.9993130  0.0022688 -0.303   0.7620
# pc06  0.0067864  1.0068095  0.0065716  1.033   0.3017
# pc07  0.0054262  1.0054410  0.0059056  0.919   0.3582
# pc08 -0.0013405  0.9986603  0.0058786 -0.228   0.8196
# pc09 -0.0026178  0.9973856  0.0022827 -1.147   0.2515
# pc10  0.0057683  1.0057850  0.0050951  1.132   0.2576
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#      exp(coef) exp(-coef) lower .95 upper .95
# h1.3    1.2437     0.8040    1.1882    1.3018
# h2      4.2100     0.2375    2.8648    6.1869
# pc01    0.9990     1.0010    0.9856    1.0126
# pc02    0.9855     1.0147    0.9718    0.9994
# pc03    1.0075     0.9926    0.9939    1.0212
# pc04    1.0036     0.9964    0.9935    1.0138
# pc05    0.9993     1.0007    0.9949    1.0038
# pc06    1.0068     0.9932    0.9939    1.0199
# pc07    1.0054     0.9946    0.9939    1.0171
# pc08    0.9987     1.0013    0.9872    1.0102
# pc09    0.9974     1.0026    0.9929    1.0019
# pc10    1.0058     0.9942    0.9958    1.0159

# Concordance= 0.526  (se = 0.003 )
# Likelihood ratio test= 129.9  on 12 df,   p=<2e-16
# Wald test            = 152.7  on 12 df,   p=<2e-16
# Score (logrank) test = 162.7  on 12 df,   p=<2e-16

# Surv(age,affected) ~ h1 + h2 + h3 + pc01 + pc02 + pc03 + pc04 + pc05 + pc06 + pc07 + pc08 + pc09 + pc10
# Call:
# coxph(formula = eval(parse(text = Y)), data = X)

#   n= 181034, number of events= 9011

#            coef  exp(coef)   se(coef)      z Pr(>|z|)
# h1    0.2146142  1.2393837  0.0233306  9.199  < 2e-16 ***
# h2    1.4378414  4.2115949  0.1964156  7.320 2.47e-13 ***
# h3    2.8182842 16.7480892  0.4084815  6.899 5.22e-12 ***
# pc01 -0.0009650  0.9990355  0.0069144 -0.140   0.8890
# pc02 -0.0144919  0.9856126  0.0071579 -2.025   0.0429 *
# pc03  0.0076005  1.0076294  0.0069013  1.101   0.2708
# pc04  0.0036878  1.0036946  0.0051798  0.712   0.4765
# pc05 -0.0007007  0.9992996  0.0022692 -0.309   0.7575
# pc06  0.0068818  1.0069055  0.0065716  1.047   0.2950
# pc07  0.0053823  1.0053968  0.0059044  0.912   0.3620
# pc08 -0.0012551  0.9987457  0.0058790 -0.213   0.8309
# pc09 -0.0026244  0.9973790  0.0022832 -1.149   0.2504
# pc10  0.0057719  1.0057886  0.0050959  1.133   0.2574
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#      exp(coef) exp(-coef) lower .95 upper .95
# h1      1.2394    0.80685    1.1840    1.2974
# h2      4.2116    0.23744    2.8659    6.1892
# h3     16.7481    0.05971    7.5208   37.2963
# pc01    0.9990    1.00097    0.9856    1.0127
# pc02    0.9856    1.01460    0.9719    0.9995
# pc03    1.0076    0.99243    0.9941    1.0214
# pc04    1.0037    0.99632    0.9936    1.0139
# pc05    0.9993    1.00070    0.9949    1.0038
# pc06    1.0069    0.99314    0.9940    1.0200
# pc07    1.0054    0.99463    0.9938    1.0171
# pc08    0.9987    1.00126    0.9873    1.0103
# pc09    0.9974    1.00263    0.9929    1.0019
# pc10    1.0058    0.99424    0.9958    1.0159

# Concordance= 0.526  (se = 0.003 )
# Likelihood ratio test= 150  on 13 df,   p=<2e-16
# Wald test            = 198.3  on 13 df,   p=<2e-16
# Score (logrank) test = 251.4  on 13 df,   p=<2e-16

# LRT p value = 7.333569e-06 on 1 degree(s) of freedom with HR = 1.674809e+01 (7.520821e+00 to 3.729626e+01), frequency = 3.038103e-05 (11) = 3.329264e-04/1.453294e-05 (6/5), p = 5.648169e-06


Rscript /home/wletsou/scripts/ukbb_backward_selection2.R file=h2+h3.ukbb_discovery.new_allele_counts.txt n=2 groups=1,2,3 threshold=0.05 n_covs=10
#    remaining_vars         pval
# 1:             h2 2.362604e-07
# 2:             h3 7.333569e-06

# final model
awk '($6<1e-5){printf "%s:",$2} END{printf "\n"}' ukbb_discovery.Significant_patterns.txt

sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69231641.txt "rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs74471298_A=0:rs2376558_C=1,rs61881030_A=0,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs12577213_A=0,rs73520512_T=0,rs3892895_G=1,rs7116054_T=0,rs74897859_T=0,rs17308715_T=0,rs67039008_A=1,rs17149775_G=1,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs12418948_T=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=1,rs10896449_G=1,rs7130881_G=0,rs111762835_A=0,rs4930672_A=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs78656034_T=0,rs74551015_A=0,rs117821797_T=0,rs12275849_C=0,rs117131950_A=0,rs7103126_C=0,rs11539762_A=0,rs142581206_T=0,rs11228610_C=0,rs12274095_A=0,rs73512137_T=0,rs12224376_A=0,rs11603814_G=0,rs72932105_T=0,rs12789955_G=0,rs117186144_T=0,rs79487139_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12790064_G=1,rs76855139_T=0,rs74471298_A=0:rs2376558_C=1,rs61881030_A=0,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs12577213_A=0,rs73520512_T=0,rs3892895_G=1,rs7116054_T=0,rs74897859_T=0,rs17308715_T=0,rs67039008_A=1,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs74551015_A=0,rs117821797_T=0,rs12275849_C=0,rs56103266_C=0,rs117131950_A=0,rs116926312_T=0,rs11539762_A=0,rs142581206_T=0,rs11228610_C=0,rs12274095_A=0,rs12224376_A=0,rs7124547_T=1,rs11603814_G=1,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs79487139_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs76855139_T=0,rs74471298_A=0" "1" "rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "1" "" h2-h4.ukbb_discovery /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts

Rscript /home/wletsou/scripts/ukbb_haplotype_model5.R file=/scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/h2-h4.ukbb_discovery.new_allele_counts.txt groups=1,2,3,4 verbose=1 n=3 n_covs=10
# [1] "X0000" "X1000" "X1001" "X1010" "X1100"
# [1] "DT[,h0:=X0000]"
# [1] "DT[,h1.4:=X1000+X1001]"
# [1] "DT[,h1:=X1000]"
# [1] "DT[,h2:=X1100]"
# [1] "DT[,h3:=X1010]"
# [1] "DT[,h4:=X1001]"

#    affected  h1.4    h1  h2 h3 h4
# 1:        0 32421 32416 116 49  5
# 2:        1  2057  2051  26 17  6

#    affected       h1.4         h1           h2           h3           h4
# 1:        0 0.09423449 0.09421996 0.0003371642 0.0001424228 1.453294e-05
# 2:        1 0.11413828 0.11380535 0.0014426812 0.0009432915 3.329264e-04

#            h1.4           h1           h2           h3           h4
# 1: 5.144621e-18 1.766029e-17 8.900272e-09 1.687137e-08 5.648169e-06

#     h1.4    h1  h2 h3 h4
# 1: 34478 34467 142 66 11

#          h1.4         h1           h2           h3           h4
# 1: 0.09522521 0.09519483 0.0003921915 0.0001822862 3.038103e-05

# Surv(age,affected) ~ h1.4 + h2 + h3 + pc01 + pc02 + pc03 + pc04 + pc05 + pc06 + pc07 + pc08 + pc09 + pc10
# Call:
# coxph(formula = eval(parse(text = Y)), data = X)

#   n= 181034, number of events= 9011

#            coef  exp(coef)   se(coef)      z Pr(>|z|)
# h1.4  0.2111493  1.2350967  0.0233848  9.029  < 2e-16 ***
# h2    1.4385348  4.2145163  0.1964159  7.324 2.41e-13 ***
# h3    1.7695637  5.8682925  0.2428330  7.287 3.17e-13 ***
# pc01 -0.0010891  0.9989115  0.0069129 -0.158   0.8748
# pc02 -0.0147107  0.9853970  0.0071568 -2.055   0.0398 *
# pc03  0.0073840  1.0074114  0.0069036  1.070   0.2848
# pc04  0.0035029  1.0035091  0.0051790  0.676   0.4988
# pc05 -0.0005938  0.9994064  0.0022690 -0.262   0.7936
# pc06  0.0067676  1.0067906  0.0065713  1.030   0.3031
# pc07  0.0052896  1.0053036  0.0059059  0.896   0.3704
# pc08 -0.0014112  0.9985898  0.0058785 -0.240   0.8103
# pc09 -0.0026333  0.9973702  0.0022830 -1.153   0.2487
# pc10  0.0055516  1.0055670  0.0050950  1.090   0.2759
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#      exp(coef) exp(-coef) lower .95 upper .95
# h1.4    1.2351     0.8097    1.1798    1.2930
# h2      4.2145     0.2373    2.8679    6.1935
# h3      5.8683     0.1704    3.6460    9.4452
# pc01    0.9989     1.0011    0.9855    1.0125
# pc02    0.9854     1.0148    0.9717    0.9993
# pc03    1.0074     0.9926    0.9939    1.0211
# pc04    1.0035     0.9965    0.9934    1.0137
# pc05    0.9994     1.0006    0.9950    1.0039
# pc06    1.0068     0.9933    0.9939    1.0198
# pc07    1.0053     0.9947    0.9937    1.0170
# pc08    0.9986     1.0014    0.9872    1.0102
# pc09    0.9974     1.0026    0.9929    1.0018
# pc10    1.0056     0.9945    0.9956    1.0157

# Concordance= 0.527  (se = 0.003 )
# Likelihood ratio test= 155.9  on 13 df,   p=<2e-16
# Wald test            = 199.2  on 13 df,   p=<2e-16
# Score (logrank) test = 224.3  on 13 df,   p=<2e-16

# Surv(age,affected) ~ h1 + h2 + h3 + h4 + pc01 + pc02 + pc03 + pc04 + pc05 + pc06 + pc07 + pc08 + pc09 + pc10
# Call:
# coxph(formula = eval(parse(text = Y)), data = X)

#   n= 181034, number of events= 9011

#            coef  exp(coef)   se(coef)      z Pr(>|z|)
# h1    0.2076314  1.2307595  0.0234194  8.866  < 2e-16 ***
# h2    1.4389142  4.2161156  0.1964161  7.326 2.37e-13 ***
# h3    1.7700066  5.8708918  0.2428332  7.289 3.12e-13 ***
# h4    2.8208394 16.7909390  0.4084818  6.906 5.00e-12 ***
# pc01 -0.0010510  0.9989495  0.0069137 -0.152   0.8792
# pc02 -0.0145911  0.9855149  0.0071579 -2.038   0.0415 *
# pc03  0.0075323  1.0075608  0.0069019  1.091   0.2751
# pc04  0.0035903  1.0035968  0.0051799  0.693   0.4882
# pc05 -0.0006071  0.9993931  0.0022693 -0.268   0.7891
# pc06  0.0068624  1.0068860  0.0065714  1.044   0.2964
# pc07  0.0052451  1.0052589  0.0059047  0.888   0.3744
# pc08 -0.0013264  0.9986745  0.0058789 -0.226   0.8215
# pc09 -0.0026400  0.9973635  0.0022836 -1.156   0.2476
# pc10  0.0055542  1.0055696  0.0050958  1.090   0.2757
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#      exp(coef) exp(-coef) lower .95 upper .95
# h1      1.2308    0.81251    1.1755    1.2886
# h2      4.2161    0.23719    2.8689    6.1959
# h3      5.8709    0.17033    3.6476    9.4494
# h4     16.7909    0.05956    7.5401   37.3917
# pc01    0.9989    1.00105    0.9855    1.0126
# pc02    0.9855    1.01470    0.9718    0.9994
# pc03    1.0076    0.99250    0.9940    1.0213
# pc04    1.0036    0.99642    0.9935    1.0138
# pc05    0.9994    1.00061    0.9950    1.0038
# pc06    1.0069    0.99316    0.9940    1.0199
# pc07    1.0053    0.99477    0.9937    1.0170
# pc08    0.9987    1.00133    0.9872    1.0102
# pc09    0.9974    1.00264    0.9929    1.0018
# pc10    1.0056    0.99446    0.9956    1.0157

# Concordance= 0.527  (se = 0.003 )
# Likelihood ratio test= 176.2  on 14 df,   p=<2e-16
# Wald test            = 244.8  on 14 df,   p=<2e-16
# Score (logrank) test = 313.2  on 14 df,   p=<2e-16

# LRT p value = 6.938959e-06 on 1 degree(s) of freedom with HR = 1.679094e+01 (7.540058e+00 to 3.739171e+01), frequency = 3.038103e-05 (11) = 3.329264e-04/1.453294e-05 (6/5), p = 5.648169e-06

Rscript /home/wletsou/scripts/ukbb_backward_selection2.R file=h2-h4.ukbb_discovery.new_allele_counts.txt n=3 groups=1,2,3,4 threshold=0.05 n_covs=10

#    remaining_vars         pval
# 1:             h2 2.030050e-07
# 2:             h3 3.186593e-07
# 3:             h4 6.938959e-06

sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69231641.txt "rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs111762835_A=0,rs4930672_A=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs142581206_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs76855139_T=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs111762835_A=0,rs4930672_A=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs76855139_T=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs111762835_A=0,rs4930672_A=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs74471298_A=0:rs2376558_C=1,rs61881030_A=0,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs12577213_A=0,rs73520512_T=0,rs3892895_G=1,rs7116054_T=0,rs74897859_T=0,rs17308715_T=0,rs67039008_A=1,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs74551015_A=0,rs117821797_T=0,rs12275849_C=0,rs56103266_C=0,rs117131950_A=0,rs116926312_T=0,rs11539762_A=0,rs142581206_T=0,rs11228610_C=0,rs12274095_A=0,rs12224376_A=0,rs7124547_T=1,rs11603814_G=1,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs79487139_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs76855139_T=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs76855139_T=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs142581206_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs111762835_A=0,rs4930672_A=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs142581206_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs142581206_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs76855139_T=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs76855139_T=0,rs74471298_A=0" "1" "rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "1" "" h2-h15.ukbb_discovery /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts
# Rscript /home/wletsou/scripts/ukbb_haplotype_model5.R file=/scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/h2-h15.ukbb_discovery.new_allele_counts.txt groups=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 verbose=1 n=14 n_covs=10

#    affected h1.15    h1  h2  h3  h4  h5  h6 h7  h8  h9 h10 h11 h12 h13 h14 h15
# 1:        0 32566 32417 108 108 107 155 148  5 100 108 101 115 116 109 101 149
# 2:        1  2090  2063  22  23  22  29  28  6  21  24  21  25  26  23  21  27

#    affected      h1.15         h1           h2           h3           h4           h5          h6           h7           h8           h9          h10          h11          h12          h13          h14          h15
# 1:        0 0.09465595 0.09422287 0.0003139115 0.0003139115 0.0003110049 0.0004505212 0.000430175 1.453294e-05 0.0002906588 0.0003139115 0.0002935654 0.0003342576 0.0003371642 0.0003168181 0.0002935654 0.0004330816
# 2:        1 0.11596937 0.11447120 0.0012207302 0.0012762180 0.0012207302 0.0016091444 0.001553657 3.329264e-04 0.0011652425 0.0013317057 0.0011652425 0.0013871934 0.0014426812 0.0012762180 0.0011652425 0.0014981689

#           h1.15           h1           h2           h3           h4           h5           h6           h7           h8           h9          h10          h11          h12          h13          h14          h15
# 1: 2.736497e-20 1.448091e-18 4.952813e-07 1.385016e-07 4.317433e-07 4.243274e-08 5.857782e-08 5.648169e-06 5.794485e-07 3.744068e-08 6.671242e-07 2.860703e-08 8.900272e-09 1.598263e-07 6.671242e-07 2.135436e-07

#    h1.15    h1  h2  h3  h4  h5  h6 h7  h8  h9 h10 h11 h12 h13 h14 h15
# 1: 34656 34480 130 131 129 184 176 11 121 132 122 140 142 132 122 176

#         h1.15         h1           h2           h3           h4           h5           h6           h7           h8           h9          h10          h11          h12          h13          h14          h15
# 1: 0.09571683 0.09523073 0.0003590486 0.0003618105 0.0003562867 0.0005081918 0.0004860965 3.038103e-05 0.0003341914 0.0003645724 0.0003369533 0.0003866677 0.0003921915 0.0003645724 0.0003369533 0.0004860965

sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr11.68850000-69231641.txt "rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs111762835_A=0,rs4930672_A=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs142581206_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs76855139_T=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs111762835_A=0,rs4930672_A=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs76855139_T=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs111762835_A=0,rs4930672_A=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs74471298_A=0:rs2376558_C=1,rs61881030_A=0,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs12577213_A=0,rs73520512_T=0,rs3892895_G=1,rs7116054_T=0,rs74897859_T=0,rs17308715_T=0,rs67039008_A=1,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs74551015_A=0,rs117821797_T=0,rs12275849_C=0,rs56103266_C=0,rs117131950_A=0,rs116926312_T=0,rs11539762_A=0,rs142581206_T=0,rs11228610_C=0,rs12274095_A=0,rs12224376_A=0,rs7124547_T=1,rs11603814_G=1,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs79487139_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs76855139_T=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs76855139_T=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs142581206_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs111762835_A=0,rs4930672_A=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs72932105_T=0,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs7102705_G=1,rs117186144_T=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs142581206_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs142581206_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs76855139_T=0,rs74471298_A=0:rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs76855139_T=0,rs74471298_A=0" "1" "rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "1" "" h2-h15.drive_replication /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2 /home/wletsou/scripts

# Rscript /home/wletsou/scripts/dbgap_haplotype_model3.R file=/scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/h2-h15.drive_replication.new_allele_counts.txt groups=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 verbose=1 n=14 n_covs=10

#    affected h1.15   h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 h11 h12 h13 h14 h15
# 1:        0  4888 4862 25 22 25 30 25  1 20 20  20  25  22  22  22  26
# 2:        1  6804 6748 55 51 55 63 55 13 49 50  49  56  52  51  50  56

#    affected      h1.15         h1           h2           h3           h4           h5           h6           h7           h8           h9          h10          h11          h12          h13          h14          h15
# 1:        0 0.09666957 0.09615537 0.0004944229 0.0004350922 0.0004944229 0.0005933075 0.0004944229 1.977692e-05 0.0003955383 0.0003955383 0.0003955383 0.0004944229 0.0004350922 0.0004350922 0.0004350922 0.0005141998
# 2:        1 0.11315859 0.11222725 0.0009147153 0.0008481905 0.0009147153 0.0010477648 0.0009147153 2.162054e-04 0.0008149282 0.0008315593 0.0008149282 0.0009313465 0.0008648217 0.0008481905 0.0008315593 0.0009313465

#           h1.15          h1          h2          h3          h4         h5          h6          h7          h8          h9         h10         h11         h12         h13        h14        h15
# 1: 5.252388e-19 3.09829e-18 0.009679642 0.009262831 0.009679642 0.00919264 0.009679642 0.005005311 0.005277603 0.003859875 0.005277603 0.007314411 0.006874466 0.009262831 0.01241295 0.01078046

#    h1.15    h1 h2 h3 h4 h5 h6 h7 h8 h9 h10 h11 h12 h13 h14 h15
# 1: 11692 11610 80 73 80 93 80 14 69 70  69  81  74  73  72  82

#        h1.15        h1           h2           h3           h4           h5           h6           h7           h8           h9          h10          h11          h12          h13          h14          h15
# 1: 0.1056264 0.1048856 0.0007227261 0.0006594876 0.0007227261 0.0008401691 0.0007227261 0.0001264771 0.0006233513 0.0006323854 0.0006233513 0.0007317602 0.0006685217 0.0006594876 0.0006504535 0.0007407943


sh /home/wletsou/scripts/vcf2bed2.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.haplotypes.hg38.vcf "rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs74471298_A=0,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "" "h2.bed" # h2, with alleles

sh /home/wletsou/scripts/vcf2bed2.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.haplotypes.hg38.vcf "rs2376558_C=1,rs61881030_A=0,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs12577213_A=0,rs73520512_T=0,rs3892895_G=1,rs7116054_T=0,rs74897859_T=0,rs17308715_T=0,rs67039008_A=1,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs74551015_A=0,rs117821797_T=0,rs12275849_C=0,rs56103266_C=0,rs117131950_A=0,rs116926312_T=0,rs11539762_A=0,rs142581206_T=0,rs11228610_C=0,rs12274095_A=0,rs12224376_A=0,rs7124547_T=1,rs11603814_G=1,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs79487139_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs76855139_T=0,rs74471298_A=0,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "" "h3.bed" # h3, with alleles

sh /home/wletsou/scripts/ukbb_haplotype_haplotype_LD2.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69231641.txt "rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs74471298_A=0,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "rs2376558_C=1,rs61881030_A=0,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs12577213_A=0,rs73520512_T=0,rs3892895_G=1,rs7116054_T=0,rs74897859_T=0,rs17308715_T=0,rs67039008_A=1,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs74551015_A=0,rs117821797_T=0,rs12275849_C=0,rs56103266_C=0,rs117131950_A=0,rs116926312_T=0,rs11539762_A=0,rs142581206_T=0,rs11228610_C=0,rs12274095_A=0,rs12224376_A=0,rs7124547_T=1,rs11603814_G=1,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs79487139_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs76855139_T=0,rs74471298_A=0,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0"

# replication

sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr11.68850000-69231641.txt "rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs74471298_A=0" "1" "rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "1" "" h2.replication
# LRT p value = 2.471105e-02 on 1 degree(s) of freedom with OR = 2.099783e+00 (1.271988e+00 to 3.466298e+00), frequency = 6.685217e-04 (74) = 8.648217e-04/4.350922e-04 (52/22), p = 6.874466e-03

sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr11.68850000-69231641.txt "rs2376558_C=1,rs61881030_A=0,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs12577213_A=0,rs73520512_T=0,rs3892895_G=1,rs7116054_T=0,rs74897859_T=0,rs17308715_T=0,rs67039008_A=1,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs74551015_A=0,rs117821797_T=0,rs12275849_C=0,rs56103266_C=0,rs117131950_A=0,rs116926312_T=0,rs11539762_A=0,rs142581206_T=0,rs11228610_C=0,rs12274095_A=0,rs12224376_A=0,rs7124547_T=1,rs11603814_G=1,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs79487139_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs76855139_T=0,rs74471298_A=0" "1" "rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "1" "" h3.replication
# LRT p value = 2.805067e-03 on 1 degree(s) of freedom with OR = 1.169116e+01 (1.530030e+00 to 8.933363e+01), frequency = 1.264771e-04 (14) = 2.162054e-04/1.977692e-05 (13/1), p = 5.005311e-03

# impute GWAS hit
mkdir ukbb_imputed_snps
cp dbgap28544_cases.indiv ukbb_imputed_snps/dbgap28544_cases.indiv
cp dbgap28544_controls.indiv ukbb_imputed_snps/dbgap28544_controls.indiv
cp ukbb_bca_cases.indiv ukbb_imputed_snps/ukbb_bca_cases.indiv
cp ukbb_bca_controls.indiv ukbb_imputed_snps/ukbb_bca_controls.indiv

cd ukbb_imputed_snps

# vcf file for UKBB
bsub -P SJLIFE -J ukbb_topmed_merge.chr11.68850000-69231641 -oo ukbb_topmed_merge.chr11.68850000-69231641.out -eo ukbb_topmed_merge.chr11.68850000-69231641.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_topmed_merge.sh \"rs2376558,rs896973,rs61881030,rs149949063,rs144905179,rs1123665,rs3829241,rs74674490,rs72930627,rs72930631,rs11228502,rs12577213,rs73520512,rs3892895,rs10896431,rs7116054,rs74897859,rs17308715,rs67039008,rs55703863,rs17149775,rs117236867,rs11228530,rs17309046,rs56134379,rs72932540,rs12418948,rs78033785,rs74342245,rs35435383,rs4495899,rs117694794,rs10792029,rs10896445,rs78208050,rs61881109,rs75590802,rs11228565,rs7937094,rs7931342,rs10896449,rs7130881,rs7119988,rs111762835,rs4930672,rs11228599,rs7940107,rs72930267,rs7946255,rs116951063,rs118004919,rs35637432,rs34737133,rs11825015,rs78656034,rs74551015,rs117821797,rs12275849,rs56103266,rs117131950,rs7103126,rs116926312,rs11539762,rs142581206,rs11228610,rs12274095,rs73512137,rs11601693,rs10896461,rs12224376,rs78919175,rs7124547,rs11603814,rs72932105,rs12789955,rs7481709,rs115223540,rs28378931,rs7102705,rs117186144,rs79487139,rs11600497,rs72932198,rs12802601,rs12796465,rs4275647,rs61882193,rs12790064,rs11263638,rs11263641,rs76855139,rs74471298,rs7105934,rs79373485,rs4444099,rs117752342,rs1122316,rs76809977,rs559664,rs657315,rs498931,rs71465432,rs183782062,rs79442425,rs111929748,rs117490805,rs79241527,rs72932500,rs57162717,rs554219\" /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr11.qced.anno.info 11 68850000,69231641 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/raw/chr11 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/ukbb_imputed_snps /home/wletsou/scripts"

# vcf file for dbGaP
bsub -P SJLIFE -J dbgap_haplotype_merge.chr11.68850000-69231641 -oo dbgap_haplotype_merge.chr11.68850000-69231641.out -eo dbgap_haplotype_merge.chr11.68850000-69231641.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/dbgap_haplotype_merge.sh dbgap28544_cases,dbgap28544_controls /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/chr11.qced.info \"rs2376558,rs896973,rs61881030,rs149949063,rs144905179,rs1123665,rs3829241,rs74674490,rs72930627,rs72930631,rs11228502,rs12577213,rs73520512,rs3892895,rs10896431,rs7116054,rs74897859,rs17308715,rs67039008,rs55703863,rs17149775,rs117236867,rs11228530,rs17309046,rs56134379,rs72932540,rs12418948,rs78033785,rs74342245,rs35435383,rs4495899,rs117694794,rs10792029,rs10896445,rs78208050,rs61881109,rs75590802,rs11228565,rs7937094,rs7931342,rs10896449,rs7130881,rs7119988,rs111762835,rs4930672,rs11228599,rs7940107,rs72930267,rs7946255,rs116951063,rs118004919,rs35637432,rs34737133,rs11825015,rs78656034,rs74551015,rs117821797,rs12275849,rs56103266,rs117131950,rs7103126,rs116926312,rs11539762,rs142581206,rs11228610,rs12274095,rs73512137,rs11601693,rs10896461,rs12224376,rs78919175,rs7124547,rs11603814,rs72932105,rs12789955,rs7481709,rs115223540,rs28378931,rs7102705,rs117186144,rs79487139,rs11600497,rs72932198,rs12802601,rs12796465,rs4275647,rs61882193,rs12790064,rs11263638,rs11263641,rs76855139,rs74471298,rs7105934,rs79373485,rs4444099,rs117752342,rs1122316,rs76809977,rs559664,rs657315,rs498931,rs71465432,rs183782062,rs79442425,rs111929748,rs117490805,rs79241527,rs72932500,rs57162717,rs554219\" 11 68850000,69231641 /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/vcfs/dbgap28544.topmedr2.cleaned.hg38.chr11.vcf.gz"

# extract imputed SNPs in range
bsub -P SJLIFE -J ukbb_haplotype_extract3.chr11.68850000-69231641 -oo ukbb_haplotype_extract3.chr11.68850000-69231641.out -eo ukbb_haplotype_extract3.chr11.68850000-69231641.err -R "rusage[mem=20000]" -q large_mem "sh /home/wletsou/scripts/ukbb_haplotype_extract3.sh ukbb_bca_cases.indiv,ukbb_bca_controls.indiv,dbgap28544_cases.indiv,dbgap28544_controls.indiv \"rs2376558,rs896973,rs61881030,rs149949063,rs144905179,rs1123665,rs3829241,rs74674490,rs72930627,rs72930631,rs11228502,rs12577213,rs73520512,rs3892895,rs10896431,rs7116054,rs74897859,rs17308715,rs67039008,rs55703863,rs17149775,rs117236867,rs11228530,rs17309046,rs56134379,rs72932540,rs12418948,rs78033785,rs74342245,rs35435383,rs4495899,rs117694794,rs10792029,rs10896445,rs78208050,rs61881109,rs75590802,rs11228565,rs7937094,rs7931342,rs10896449,rs7130881,rs7119988,rs111762835,rs4930672,rs11228599,rs7940107,rs72930267,rs7946255,rs116951063,rs118004919,rs35637432,rs34737133,rs11825015,rs78656034,rs74551015,rs117821797,rs12275849,rs56103266,rs117131950,rs7103126,rs116926312,rs11539762,rs142581206,rs11228610,rs12274095,rs73512137,rs11601693,rs10896461,rs12224376,rs78919175,rs7124547,rs11603814,rs72932105,rs12789955,rs7481709,rs115223540,rs28378931,rs7102705,rs117186144,rs79487139,rs11600497,rs72932198,rs12802601,rs12796465,rs4275647,rs61882193,rs12790064,rs11263638,rs11263641,rs76855139,rs74471298,rs7105934,rs79373485,rs4444099,rs117752342,rs1122316,rs76809977,rs559664,rs657315,rs498931,rs71465432,rs183782062,rs79442425,rs111929748,rs117490805,rs79241527,rs72932500,rs57162717,rs554219\" /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr11.qced.anno.info 11 68850000,69231641 ukbb.topmed.hg19_chr11.68850000-69231641.hg38.vcf.gz,dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.hg38.vcf.gz /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/ukbb_imputed_snps /home/wletsou/scripts"

awk 'FILENAME==ARGV[1]{seen[$1]; next} FILENAME==ARGV[2]{seen[$1]; next} ($1 in seen || FNR==1){print $0}' ukbb_bca_cases.indiv ukbb_bca_controls.indiv haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.chr11.68850000-69231641.txt > haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69231641.txt

awk 'FILENAME==ARGV[1]{seen[$1]; next} FILENAME==ARGV[2]{seen[$1]; next} ($1 in seen || FNR==1){print $0}' dbgap28544_cases.indiv dbgap28544_controls.indiv haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.chr11.68850000-69231641.txt > haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr11.68850000-69231641.txt

# discovery
sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69231641.txt "rs554219_G=1" "" "" "" "" rs554219.discovery
# LRT p value = 3.009969e-22 on 1 degree(s) of freedom with HR = 1.235763e+00 (1.185297e+00 to 1.288379e+00), frequency = 1.194444e-01 (43247) = 1.424925e-01/1.182371e-01 (2568/40679), p = 1.286710e-21

# DRIVE replication
sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr11.68850000-69231641.txt "rs554219_G=1" "" "" "" "" rs554219.replication
# LLRT p value = 7.858252e-27 on 1 degree(s) of freedom with OR = 1.210425e+00 (1.168739e+00 to 1.253598e+00), frequency = 1.367307e-01 (15135) = 1.463877e-01/1.252472e-01 (8802/6333), p = 1.663833e-24

# LD, h1 (reduced) vs GWAS hit
sh /home/wletsou/scripts/ukbb_haplotype_haplotype_LD2.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69231641.txt "rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "rs554219_G=1" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/ukbb_imputed_snps
# Haplotype 0 vs. haplotype 0:
# 0.0957997       0.119444        0.0956505
# rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0       rs554219_G=1    0.778311        0.998232

sh /home/wletsou/scripts/vcf2bed.v3.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.haplotypes.hg38.vcf "rs554219_G=1" "" "rs554219.bed" "0" # rs554219, with alleles, no header

awk 'BEGIN{OFS="\t"} {$9="0,0,0"}1' rs554219.bed > rs554219.tmp && mv rs554219.tmp rs554219.bed # blue to black

sh /home/wletsou/scripts/vcf2bed.v3.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.haplotypes.hg38.vcf "rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs74471298_A=0,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "" "h2.bed" "0" # h2, with alleles, no header

awk 'BEGIN{OFS="\t"} ("rs2376558_C=1,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs72930627_A=0,rs72930631_C=0,rs11228502_C=0,rs12577213_A=0,rs73520512_T=0,rs10896431_T=0,rs74897859_T=0,rs55703863_C=0,rs17149775_G=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=0,rs10896449_G=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs117821797_T=0,rs56103266_C=0,rs117131950_A=0,rs7103126_C=1,rs116926312_T=0,rs12274095_A=0,rs10896461_G=1,rs78919175_T=0,rs7124547_T=1,rs12789955_G=0,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs74471298_A=0,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" ~ $3){print $1,$2,$3,$4,$5}' ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.haplotypes.hg38.vcf
# h2 rows of vcf file

# chr11   69083946        rs2376558       T       C
# chr11   69085685        rs149949063     A       G
# chr11   69085729        rs144905179     G       A
# chr11   69092678        rs74674490      A       G
# chr11   69099653        rs72930627      C       A
# chr11   69100367        rs72930631      T       C
# chr11   69110120        rs11228502      T       C
# chr11   69114683        rs12577213      C       A
# chr11   69114776        rs73520512      C       T
# chr11   69119549        rs10896431      C       T
# chr11   69126381        rs74897859      C       T
# chr11   69133728        rs55703863      G       C
# chr11   69142111        rs17149775      A       G
# chr11   69145367        rs117236867     A       G
# chr11   69151517        rs11228530      A       C
# chr11   69151962        rs17309046      T       C
# chr11   69153340        rs56134379      G       A
# chr11   69154575        rs72932540      A       G
# chr11   69168018        rs78033785      C       T
# chr11   69171308        rs74342245      G       A
# chr11   69180342        rs35435383      G       A
# chr11   69195706        rs117694794     C       T
# chr11   69199415        rs10792029      A       G
# chr11   69201946        rs78208050      A       G
# chr11   69209346        rs61881109      G       A
# chr11   69210118        rs75590802      G       A
# chr11   69212239        rs7937094       C       T
# chr11   69227030        rs7931342       T       G
# chr11   69227200        rs10896449      A       G
# chr11   69228491        rs7130881       A       G
# chr11   69248404        rs7119988       G       A
# chr11   69250616        rs111762835     G       A
# chr11   69251764        rs4930672       G       A
# chr11   69258819        rs11228599      G       T
# chr11   69260303        rs7940107       G       A
# chr11   69263523        rs72930267      A       T
# chr11   69263730        rs7946255       T       C
# chr11   69265264        rs116951063     C       T
# chr11   69267799        rs118004919     G       A
# chr11   69269130        rs35637432      G       A
# chr11   69275141        rs34737133      C       T
# chr11   69277004        rs11825015      C       T
# chr11   69285720        rs78656034      C       T
# chr11   69287175        rs117821797     C       T
# chr11   69289127        rs56103266      G       C
# chr11   69295686        rs117131950     G       A
# chr11   69295926        rs7103126       T       C
# chr11   69296009        rs116926312     C       T
# chr11   69296300        rs12274095      C       A
# chr11   69306180        rs10896461      A       G
# chr11   69309047        rs78919175      C       T
# chr11   69309919        rs7124547       C       T
# chr11   69311667        rs12789955      T       G
# chr11   69321285        rs115223540     C       T
# chr11   69321751        rs28378931      A       G
# chr11   69338687        rs117186144     C       T
# chr11   69360358        rs72932198      G       A
# chr11   69361135        rs12802601      T       C
# chr11   69367952        rs12796465      A       G
# chr11   69377704        rs61882193      A       G
# chr11   69382816        rs12790064      A       G
# chr11   69386345        rs11263638      G       A
# chr11   69414699        rs74471298      G       A
# chr11   69424973        rs7105934       G       A
# chr11   69432853        rs79373485      C       T
# chr11   69458108        rs4444099       C       T
# chr11   69475035        rs117752342     C       T
# chr11   69485307        rs1122316       G       A
# chr11   69508496        rs76809977      C       T
# chr11   69517492        rs559664        A       G
# chr11   69517944        rs657315        C       T
# chr11   69524849        rs498931        A       G
# chr11   69528387        rs71465432      C       T
# chr11   69539114        rs183782062     G       T
# chr11   69562741        rs79442425      C       T
# chr11   69564519        rs111929748     G       T
# chr11   69572513        rs117490805     G       A
# chr11   69581315        rs79241527      C       T
# chr11   69603668        rs72932500      A       G
# chr11   69609602        rs57162717      A       T

sh /home/wletsou/scripts/vcf2bed.v3.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.haplotypes.hg38.vcf "rs2376558_C=1,rs61881030_A=0,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs12577213_A=0,rs73520512_T=0,rs3892895_G=1,rs7116054_T=0,rs74897859_T=0,rs17308715_T=0,rs67039008_A=1,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs74551015_A=0,rs117821797_T=0,rs12275849_C=0,rs56103266_C=0,rs117131950_A=0,rs116926312_T=0,rs11539762_A=0,rs142581206_T=0,rs11228610_C=0,rs12274095_A=0,rs12224376_A=0,rs7124547_T=1,rs11603814_G=1,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs79487139_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs76855139_T=0,rs74471298_A=0,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "" "h3.bed" "0" # h3, with alleles, no header

awk 'BEGIN{OFS="\t"} ("rs2376558_C=1,rs61881030_A=0,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs12577213_A=0,rs73520512_T=0,rs3892895_G=1,rs7116054_T=0,rs74897859_T=0,rs17308715_T=0,rs67039008_A=1,rs55703863_C=0,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs11228565_A=0,rs7937094_T=0,rs7130881_G=0,rs7119988_A=0,rs111762835_A=0,rs4930672_A=0,rs11228599_T=0,rs7940107_A=0,rs72930267_T=0,rs7946255_C=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs11825015_T=0,rs78656034_T=0,rs74551015_A=0,rs117821797_T=0,rs12275849_C=0,rs56103266_C=0,rs117131950_A=0,rs116926312_T=0,rs11539762_A=0,rs142581206_T=0,rs11228610_C=0,rs12274095_A=0,rs12224376_A=0,rs7124547_T=1,rs11603814_G=1,rs115223540_T=0,rs28378931_G=0,rs117186144_T=0,rs79487139_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12796465_G=0,rs4275647_T=0,rs61882193_G=0,rs12790064_G=1,rs11263638_A=0,rs11263641_C=0,rs76855139_T=0,rs74471298_A=0,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" ~ $3){print $1,$2,$3,$4,$5}' ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.haplotypes.hg38.vcf # h3 rows of vcf file

# chr11   69083946        rs2376558       T       C
# chr11   69085616        rs61881030      C       A
# chr11   69085685        rs149949063     A       G
# chr11   69085729        rs144905179     G       A
# chr11   69092678        rs74674490      A       G
# chr11   69114683        rs12577213      C       A
# chr11   69114776        rs73520512      C       T
# chr11   69117287        rs3892895       A       G
# chr11   69124227        rs7116054       C       T
# chr11   69126381        rs74897859      C       T
# chr11   69127338        rs17308715      C       T
# chr11   69130541        rs67039008      G       A
# chr11   69133728        rs55703863      G       C
# chr11   69145367        rs117236867     A       G
# chr11   69151517        rs11228530      A       C
# chr11   69151962        rs17309046      T       C
# chr11   69153340        rs56134379      G       A
# chr11   69168018        rs78033785      C       T
# chr11   69171308        rs74342245      G       A
# chr11   69180342        rs35435383      G       A
# chr11   69195706        rs117694794     C       T
# chr11   69199415        rs10792029      A       G
# chr11   69201946        rs78208050      A       G
# chr11   69209346        rs61881109      G       A
# chr11   69210118        rs75590802      G       A
# chr11   69211113        rs11228565      G       A
# chr11   69212239        rs7937094       C       T
# chr11   69228491        rs7130881       A       G
# chr11   69248404        rs7119988       G       A
# chr11   69250616        rs111762835     G       A
# chr11   69251764        rs4930672       G       A
# chr11   69258819        rs11228599      G       T
# chr11   69260303        rs7940107       G       A
# chr11   69263523        rs72930267      A       T
# chr11   69263730        rs7946255       T       C
# chr11   69265264        rs116951063     C       T
# chr11   69267799        rs118004919     G       A
# chr11   69269130        rs35637432      G       A
# chr11   69275141        rs34737133      C       T
# chr11   69277004        rs11825015      C       T
# chr11   69285720        rs78656034      C       T
# chr11   69286582        rs74551015      G       A
# chr11   69287175        rs117821797     C       T
# chr11   69287948        rs12275849      T       C
# chr11   69289127        rs56103266      G       C
# chr11   69295686        rs117131950     G       A
# chr11   69296009        rs116926312     C       T
# chr11   69296043        rs11539762      G       A
# chr11   69296259        rs142581206     G       T
# chr11   69296261        rs11228610      G       C
# chr11   69296300        rs12274095      C       A
# chr11   69306623        rs12224376      G       A
# chr11   69309919        rs7124547       C       T
# chr11   69311104        rs11603814      C       G
# chr11   69321285        rs115223540     C       T
# chr11   69321751        rs28378931      A       G
# chr11   69338687        rs117186144     C       T
# chr11   69352233        rs79487139      C       T
# chr11   69358550        rs11600497      C       A
# chr11   69360358        rs72932198      G       A
# chr11   69361135        rs12802601      T       C
# chr11   69367952        rs12796465      A       G
# chr11   69371038        rs4275647       C       T
# chr11   69377704        rs61882193      A       G
# chr11   69382816        rs12790064      A       G
# chr11   69386345        rs11263638      G       A
# chr11   69392822        rs11263641      T       C
# chr11   69409182        rs76855139      C       T
# chr11   69414699        rs74471298      G       A
# chr11   69424973        rs7105934       G       A
# chr11   69432853        rs79373485      C       T
# chr11   69458108        rs4444099       C       T
# chr11   69475035        rs117752342     C       T
# chr11   69485307        rs1122316       G       A
# chr11   69508496        rs76809977      C       T
# chr11   69517492        rs559664        A       G
# chr11   69517944        rs657315        C       T
# chr11   69524849        rs498931        A       G
# chr11   69528387        rs71465432      C       T
# chr11   69539114        rs183782062     G       T
# chr11   69562741        rs79442425      C       T
# chr11   69564519        rs111929748     G       T
# chr11   69572513        rs117490805     G       A
# chr11   69581315        rs79241527      C       T
# chr11   69603668        rs72932500      A       G
# chr11   69609602        rs57162717      A       T

sh /home/wletsou/scripts/vcf2bed.v3.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.haplotypes.hg38.vcf "rs2376558_C=1,rs61881030_A=0,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs12577213_A=0,rs73520512_T=0,rs3892895_G=1,rs7116054_T=0,rs74897859_T=0,rs17308715_T=0,rs67039008_A=1,rs17149775_G=1,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs12418948_T=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=1,rs10896449_G=1,rs7130881_G=0,rs111762835_A=0,rs4930672_A=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs78656034_T=0,rs74551015_A=0,rs117821797_T=0,rs12275849_C=0,rs117131950_A=0,rs7103126_C=0,rs11539762_A=0,rs142581206_T=0,rs11228610_C=0,rs12274095_A=0,rs73512137_T=0,rs12224376_A=0,rs11603814_G=0,rs72932105_T=0,rs12789955_G=0,rs117186144_T=0,rs79487139_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12790064_G=1,rs76855139_T=0,rs74471298_A=0,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "" "g4.bed" "0" # non-replicated haplotype g4, with alleles, no header

awk 'BEGIN{OFS="\t"} ("rs2376558_C=1,rs61881030_A=0,rs149949063_G=0,rs144905179_A=0,rs74674490_G=0,rs12577213_A=0,rs73520512_T=0,rs3892895_G=1,rs7116054_T=0,rs74897859_T=0,rs17308715_T=0,rs67039008_A=1,rs17149775_G=1,rs117236867_G=0,rs11228530_C=0,rs17309046_C=0,rs56134379_A=0,rs72932540_G=0,rs12418948_T=0,rs78033785_T=0,rs74342245_A=0,rs35435383_A=0,rs117694794_T=0,rs10792029_G=1,rs78208050_G=0,rs61881109_A=0,rs75590802_A=0,rs7937094_T=0,rs7931342_G=1,rs10896449_G=1,rs7130881_G=0,rs111762835_A=0,rs4930672_A=0,rs116951063_T=0,rs118004919_A=0,rs35637432_A=0,rs34737133_T=0,rs78656034_T=0,rs74551015_A=0,rs117821797_T=0,rs12275849_C=0,rs117131950_A=0,rs7103126_C=0,rs11539762_A=0,rs142581206_T=0,rs11228610_C=0,rs12274095_A=0,rs73512137_T=0,rs12224376_A=0,rs11603814_G=0,rs72932105_T=0,rs12789955_G=0,rs117186144_T=0,rs79487139_T=0,rs11600497_A=0,rs72932198_A=0,rs12802601_C=0,rs12790064_G=1,rs76855139_T=0,rs74471298_A=0,rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" ~ $3){print $1,$2,$3,$4,$5}' ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.haplotypes.hg38.vcf # g4 rows of vcf file

# chr11   69083946        rs2376558       T       C
# chr11   69085616        rs61881030      C       A
# chr11   69085685        rs149949063     A       G
# chr11   69085729        rs144905179     G       A
# chr11   69092678        rs74674490      A       G
# chr11   69114683        rs12577213      C       A
# chr11   69114776        rs73520512      C       T
# chr11   69117287        rs3892895       A       G
# chr11   69124227        rs7116054       C       T
# chr11   69126381        rs74897859      C       T
# chr11   69127338        rs17308715      C       T
# chr11   69130541        rs67039008      G       A
# chr11   69142111        rs17149775      A       G
# chr11   69145367        rs117236867     A       G
# chr11   69151517        rs11228530      A       C
# chr11   69151962        rs17309046      T       C
# chr11   69153340        rs56134379      G       A
# chr11   69154575        rs72932540      A       G
# chr11   69159054        rs12418948      C       T
# chr11   69168018        rs78033785      C       T
# chr11   69171308        rs74342245      G       A
# chr11   69180342        rs35435383      G       A
# chr11   69195706        rs117694794     C       T
# chr11   69199415        rs10792029      A       G
# chr11   69201946        rs78208050      A       G
# chr11   69209346        rs61881109      G       A
# chr11   69210118        rs75590802      G       A
# chr11   69212239        rs7937094       C       T
# chr11   69227030        rs7931342       T       G
# chr11   69227200        rs10896449      A       G
# chr11   69228491        rs7130881       A       G
# chr11   69250616        rs111762835     G       A
# chr11   69251764        rs4930672       G       A
# chr11   69265264        rs116951063     C       T
# chr11   69267799        rs118004919     G       A
# chr11   69269130        rs35637432      G       A
# chr11   69275141        rs34737133      C       T
# chr11   69285720        rs78656034      C       T
# chr11   69286582        rs74551015      G       A
# chr11   69287175        rs117821797     C       T
# chr11   69287948        rs12275849      T       C
# chr11   69295686        rs117131950     G       A
# chr11   69295926        rs7103126       T       C
# chr11   69296043        rs11539762      G       A
# chr11   69296259        rs142581206     G       T
# chr11   69296261        rs11228610      G       C
# chr11   69296300        rs12274095      C       A
# chr11   69296517        rs73512137      C       T
# chr11   69306623        rs12224376      G       A
# chr11   69311104        rs11603814      C       G
# chr11   69311184        rs72932105      C       T
# chr11   69311667        rs12789955      T       G
# chr11   69338687        rs117186144     C       T
# chr11   69352233        rs79487139      C       T
# chr11   69358550        rs11600497      C       A
# chr11   69360358        rs72932198      G       A
# chr11   69361135        rs12802601      T       C
# chr11   69382816        rs12790064      A       G
# chr11   69409182        rs76855139      C       T
# chr11   69414699        rs74471298      G       A
# chr11   69424973        rs7105934       G       A
# chr11   69432853        rs79373485      C       T
# chr11   69458108        rs4444099       C       T
# chr11   69475035        rs117752342     C       T
# chr11   69485307        rs1122316       G       A
# chr11   69508496        rs76809977      C       T
# chr11   69517492        rs559664        A       G
# chr11   69517944        rs657315        C       T
# chr11   69524849        rs498931        A       G
# chr11   69528387        rs71465432      C       T
# chr11   69539114        rs183782062     G       T
# chr11   69562741        rs79442425      C       T
# chr11   69564519        rs111929748     G       T
# chr11   69572513        rs117490805     G       A
# chr11   69581315        rs79241527      C       T
# chr11   69603668        rs72932500      A       G
# chr11   69609602        rs57162717      A       T


sh /home/wletsou/scripts/vcf2bed.v3.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.haplotypes.hg38.vcf "rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" "" "h1.bed" "0" # h1, with alleles, no header

awk 'BEGIN{OFS="\t"} ("rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0" ~ $3){print $1,$2,$3,$4,$5}' ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.haplotypes.hg38.vcf # h1 rows of vcf file

# chr11   69424973        rs7105934       G       A
# chr11   69432853        rs79373485      C       T
# chr11   69458108        rs4444099       C       T
# chr11   69475035        rs117752342     C       T
# chr11   69485307        rs1122316       G       A
# chr11   69508496        rs76809977      C       T
# chr11   69517492        rs559664        A       G
# chr11   69517944        rs657315        C       T
# chr11   69524849        rs498931        A       G
# chr11   69528387        rs71465432      C       T
# chr11   69539114        rs183782062     G       T
# chr11   69562741        rs79442425      C       T
# chr11   69564519        rs111929748     G       T
# chr11   69572513        rs117490805     G       A
# chr11   69581315        rs79241527      C       T
# chr11   69603668        rs72932500      A       G
# chr11   69609602        rs57162717      A       T

module load ucsc/112922

fetchChromSizes hg38 > hg38.sizes # chrom sizes for bigBed file

# convert bed to bigBed
bedToBigBed rs554219.bed hg38.sizes rs554219.bb
bedToBigBed h1.bed hg38.sizes h1.bb
bedToBigBed h2.bed hg38.sizes h2.bb
bedToBigBed h3.bed hg38.sizes h3.bb
bedToBigBed g4.bed hg38.sizes g4.bb

# move to ClusterHome
cp rs554219.bed /home/wletsou/Chromosome_Overlap_results/UKBB_chr11.40/phase2/ukbb_imputed_snps/rs554219.bed
cp h1.bed /home/wletsou/Chromosome_Overlap_results/UKBB_chr11.40/phase2/ukbb_imputed_snps/h1.bed
cp h2.bed /home/wletsou/Chromosome_Overlap_results/UKBB_chr11.40/phase2/ukbb_imputed_snps/h2.bed
cp h3.bed /home/wletsou/Chromosome_Overlap_results/UKBB_chr11.40/phase2/ukbb_imputed_snps/h3.bed
cp g4.bed /home/wletsou/Chromosome_Overlap_results/UKBB_chr11.40/phase2/ukbb_imputed_snps/g4.bed

cp rs554219.bb /home/wletsou/Chromosome_Overlap_results/UKBB_chr11.40/phase2/ukbb_imputed_snps/rs554219.bb
cp h1.bb /home/wletsou/Chromosome_Overlap_results/UKBB_chr11.40/phase2/ukbb_imputed_snps/h1.bb
cp h2.bb /home/wletsou/Chromosome_Overlap_results/UKBB_chr11.40/phase2/ukbb_imputed_snps/h2.bb
cp h3.bb /home/wletsou/Chromosome_Overlap_results/UKBB_chr11.40/phase2/ukbb_imputed_snps/h3.bb
cp g4.bb /home/wletsou/Chromosome_Overlap_results/UKBB_chr11.40/phase2/ukbb_imputed_snps/g4.bb


module load samtools/1.2
module load tabix/0.2.6

bgzip -c rs554219.bed > rs554219.gz
tabix -p bed rs554219.gz

cp rs554219.gz /home/wletsou/Chromosome_Overlap_results/UKBB_chr11.40/rs554219.gz
cp rs554219.gz.tbi /home/wletsou/Chromosome_Overlap_results/UKBB_chr11.40/rs554219.gz.tbi

cd ..
# Permutation analysis

mkdir permutation1
cp dbgap28544_cases.indiv permutation1/dbgap28544_cases.indiv
cp dbgap28544_controls.indiv permutation1/dbgap28544_controls.indiv
cd permutation1

# permute UKBB case-and control-h1-carriers
bsub -P SJLIFE -J ukbb_haplotype_permute.v2 -oo ukbb_haplotype_permute.v2.out -eo ukbb_haplotype_permute.v2.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_permute.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69231641.txt rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0 20200116 ukbb_bca_cases,ukbb_bca_controls ukb.bca.pheno.txt /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/permutation1 /home/wletsou/scripts"

# extract haplotypes of permuted cases and controls
bsub -P SJLIFE -J ukbb_bca_overlap.chr11.68850000-69231641 -oo ukbb_bca_overlap.chr11.68850000-69231641.out -eo ukbb_bca_overlap.chr11.68850000-69231641.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_submission.v2.sh ukbb_bca_cases,ukbb_bca_controls dbgap28544_cases,dbgap28544_controls \"\" 11 68850000,69231641 rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0 1 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/ukbb_snp_list.chr11.68850000-69231641.txt /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/ukbb.topmed.hg19_chr11.68850000-69231641.hg38.vcf.gz,/scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/dbgap28544_cases+dbgap28544_controls.hg19_chr11.68850000-69231641.hg38.vcf.gz /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/topmed/qc/chr11.qced.anno.info \"\" \"\" \"\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/permutation1 /home/wletsou/scripts"

awk 'NR==1{char=length(NF-1)} NR>1{printf "1\t"; for (i=2;i<=NF;i++) {printf "%s%0"char"d_%s",(i>2?",":""),i-1,$i}; printf "\n"}' haplotype_estimates.ukbb_bca_cases.chr11.68850000-69231641.subset.txt | awk 'BEGIN{OFS="\t"} {seen[$2]+=$1} END{for (i in seen) {print seen[i],i} }' > Pattern_combined.Iteration000.chr11.68850000-69231641_2,j.txt # starting with all chromosomes, not just pairs

# First round of overlaps:
bsub -P SJLIFE -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr11.68850000-69231641 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr11.68850000-69231641.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration000.chr11.68850000-69231641.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ChromosomeOverlap_iteration_sub_parallel.v3.sh chr11.68850000-69231641 2 2,j 50 0 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/permutation1 /home/wletsou/scripts"

# p-values for first round of overlaps
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr11.68850000-69231641.subset.txt,haplotype_estimates.ukbb_bca_controls.chr11.68850000-69231641.subset.txt Closed_patterns_stats.chr11.68850000-69231641_2,j.txt 50 \"1\" \"Iteration001.chr11.68850000-69231641\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/permutation1 /home/wletsou/scripts"

# remove older results if staring again
rm -r *Iteration00[2-9]*

awk '($1<2){print $0}' Closed_patterns_stats.chr11.68850000-69231641_2,j.txt > Closed_patterns_stats.chr11.68850000-69231641_2,j.tmp
mv Closed_patterns_stats.chr11.68850000-69231641_2,j.tmp Closed_patterns_stats.chr11.68850000-69231641_2,j.txt

test -f Pattern_combined_old.Iteration001.chr11.68850000-69231641_2,j.txt && mv Pattern_combined_old.Iteration001.chr11.68850000-69231641_2,j.txt Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.txt
awk 'NR==FNR{seen[$1]; next} ($2 in seen){print $0}' <(awk '(NR>1 && $5+0<2.9e-4 && $4+0>1){print $0}' fisher_exact.Iteration001.chr11.68850000-69231641.patterns_000001-246392.txt | awk 'BEGIN{array[0]} {delete array; n=split($1,array,","); if (n<snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]+0 || snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]+0==0) {snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]=n; seen[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]=$1} } END{for (i in seen) {print seen[i]} }') Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.txt > Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.tmp # take shortest pattern of a family with the same p-value
mv Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.txt Pattern_combined_old.Iteration001.chr11.68850000-69231641_2,j.txt
mv Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.tmp Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.txt

# Begin iterations:
bsub -P SJLIFE -J ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr11.68850000-69231641 -oo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr11.68850000-69231641.out -eo ChromosomeOverlap_iteration_sub_parallel.v3.Iteration001.chr11.68850000-69231641.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ChromosomeOverlap_iteration_sub_parallel.v3.sh chr11.68850000-69231641 2 2,j 50 1 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/permutation1 /home/wletsou/scripts"

# (filtered) patterns appearing for the first time at Iterations 001-005
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1){print $0}' Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.txt Closed_patterns_stats.chr11.68850000-69231641_2,j.txt > Closed_patterns_stats.Iteration001-005.chr11.68850000-69231641_2,j.txt

# p-values for Iterations 001-005
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr11.68850000-69231641.subset.txt,haplotype_estimates.ukbb_bca_controls.chr11.68850000-69231641.subset.txt Closed_patterns_stats.Iteration001-005.chr11.68850000-69231641_2,j.txt 50 \"\" \"Iteration001-005.chr11.68850000-69231641\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/permutation1 /home/wletsou/scripts"

# top closed patterns
awk '($5+0<1e-4 && $4+0>1){print $0}' fisher_exact.Iteration001-005.chr11.68850000-69231641.patterns_00001-80462.txt > fisher_exact.Iteration001-005.chr11.68850000-69231641.patterns_00001-80462.Results.txt

# tranlate haplotype to rsid
bsub -P SJLIFE -J ukbb_haplotype_translate2_sub -oo ukbb_haplotype_translate2_sub.out -eo ukbb_haplotype_translate2_sub.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr11.68850000-69231641.subset.txt fisher_exact.Iteration001-005.chr11.68850000-69231641.patterns_00001-80462.Results.txt 50"

# precompute haplotype counts
bsub -P SJLIFE -J ukbb_haplotype_model9_sub.v2.ukbb_discovery -oo ukbb_haplotype_model9_sub.v2.ukbb_discovery.out -eo ukbb_haplotype_model9_sub.v2.ukbb_discovery.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/permutation1/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr11.68850000-69231641.txt fisher_exact.Iteration001-005.chr11.68850000-69231641.patterns_00001-80462.Results.translated.txt \"1\" \"rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0\" \"1\" 11 68850000,69231641 50 \"ukbb_discovery\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/permutation1 /home/wletsou/scripts"

bsub -P SJLIFE -J ukbb_haplotype_model9_sub.v2.drive_replication -oo ukbb_haplotype_model9_sub.v2.drive_replication.out -eo ukbb_haplotype_model9_sub.v2.drive_replication.err -R "rusage[mem=256]" -w "done(179338829)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.v2.sh /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr11.68850000-69231641.txt fisher_exact.Iteration001-005.chr11.68850000-69231641.patterns_00001-80462.Results.translated.txt \"1\" \"rs7105934_A=0,rs79373485_T=0,rs4444099_T=0,rs117752342_T=0,rs1122316_A=0,rs76809977_T=0,rs559664_G=1,rs657315_T=0,rs498931_G=0,rs71465432_T=0,rs183782062_T=0,rs79442425_T=0,rs111929748_T=0,rs117490805_A=0,rs79241527_T=0,rs72932500_G=0,rs57162717_T=0\" \"1\" 11 68850000,69231641 50 \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/permutation1 /home/wletsou/scripts"

# forward selection
bsub -P SJLIFE -J ukbb_haplotype_model9_iterate.v2.discovery -oo ukbb_haplotype_model9_iterate.v2.discovery.out -eo ukbb_haplotype_model9_iterate.v2.discovery.err -R "rusage[mem=256]" -w "done(179338831)" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.v2.sh ukbb_discovery.allele_counts.txt /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/permutation1/ukb.bca.pheno.txt \"1\" \"1\" 1e-5 ukbb_discovery.Significant_patterns.txt \"ukbb_discovery\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/permutation1 /home/wletsou/scripts"

bsub -P SJLIFE -J dbgap28544_haplotype_model9_iterate.v2.replication -oo dbgap28544_haplotype_model9_iterate.v2.replication.out -eo dbgap28544_haplotype_model9_iterate.v2.replication.err -w "done(179338832)" -R "rusage[mem=256]" "sh /home/wletsou/scripts/dbgap28544_haplotype_model9_iterate.v2.sh drive_replication.allele_counts.txt /research_jude/rgs01_jude/groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt \"1\" \"1\" 1e-5 drive_replication.Significant_patterns.txt \"drive_replication\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/permutation1 /home/wletsou/scripts"

# replication
awk -v alpha=0.05 'BEGIN{OFS = "\t"} NR==FNR && $2 > 1 && $5 < 1e-5{seen[$1]=$0; n+=1; next} ($1 in seen && $2 > 1 && $5 < alpha){m+=1} END{printf("%s of %s replicated at p < %0.2e\n",m + 0,n + 0,alpha)}' ukbb_discovery.Conditional_haplotype_effects.h2.txt drive_replication.Conditional_haplotype_effects.h2.txt
# 0 of 0 replicated at p < 5.00e-02

# Including contiguous haplotypes (Iteration000)

# (filtered) patterns appearing for the first time at Iterations 000,001-005
awk 'BEGIN{OFS="\t"} NR==FNR{seen[$2]; next} ($3 in seen || $1>1 || $1 == 0){print $0}' Pattern_combined.Iteration001.chr11.68850000-69231641_2,j.txt Closed_patterns_stats.chr11.68850000-69231641_2,j.txt > Closed_patterns_stats.Iteration000+Iteration001-005.chr11.68850000-69231641_2,j.txt

# p-values for Iterations 000,001-005
bsub -P SJLIFE -J ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005 -oo ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005.out -eo ChromosomeOverlap_haplotype_count_sub.v5.Iteration000+Iteration001-005.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ChromosomeOverlap_haplotype_count_sub.v5.sh haplotype_estimates.ukbb_bca_cases.chr11.68850000-69231641.subset.txt,haplotype_estimates.ukbb_bca_controls.chr11.68850000-69231641.subset.txt Closed_patterns_stats.Iteration000+Iteration001-005.chr11.68850000-69231641_2,j.txt 50 \"\" \"Iteration000+Iteration001-005.chr11.68850000-69231641\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr11.40/phase2/permutation1 /home/wletsou/scripts"

# top closed patterns
awk '($5+0<1e-4 && $4+0>1){print $0}' fisher_exact.Iteration000+Iteration001-005.chr11.68850000-69231641.patterns_00001-81167.txt > fisher_exact.Iteration000+Iteration001-005.chr11.68850000-69231641.patterns_00001-81167.Results.txt

# minimum p-value for contiguous patterns
awk 'BEGIN{min=1} NR==FNR{seen[$2]; next} ($1 in seen){ if ($5<min+0) {min=$5} } END{print min}' Pattern_combined.Iteration000.chr11.68850000-69231641_2,j.txt fisher_exact.Iteration000+Iteration001-005.chr11.68850000-69231641.patterns_00001-81167.txt
# 0.000846287480800234