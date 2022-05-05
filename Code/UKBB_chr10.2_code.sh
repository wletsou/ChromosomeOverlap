# mkdir UKBB_chr10.1 # create project folder in your chosen directory

# set your directories
# DIRECTORY=UKBB_chr10.1
# WORKING_DIR=/research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/ChromosomeOverlapScripts

# cd $DIRECTORY

# extract UKBB-genotyped SNPs in entire TAD
# mkdir ukbb_genotyped_snps
# cd ${DIRECTORY}/ukbb_genotyped_snps # move down one level

# subjects lists
# awk 'BEGIN{OFS="\t"} ($3==1){print $1,$1}' /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_cases.indiv
# awk 'BEGIN{OFS="\t"} ($3==0){print $1,$1}' /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_controls.indiv

# bsub -P SJLIFE -J ukbb_hybrid_haplotype2.chr10.123000000-124000000 -oo ukbb_hybrid_haplotype2.chr10.123000000-124000000.out -eo ukbb_hybrid_haplotype2.chr10.123000000-124000000.err -R "rusage[mem=20000]" "sh /home/wletsou/scripts/ukbb_hybrid_haplotype2.sh ukbb_bca_cases.indiv,ukbb_bca_controls.indiv 10 123000000,124000000 /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/biobank/gds/ukb.bca.hap/ukb.bca.hap.chr10.vcf.gz" # requires gds to be converted into vcf, see ukbb_gds2vcf.R if not already done; if running for the first time, use ukb.bca.hap.chr10.vcf.gz as the vcf file, otherwise ukb.bca.hap.chr10.new.vcf.gz

# bsub -P SJLIFE -J ukbb_hybrid_haplotype2.chr10.122000000-124000000 -oo ukbb_hybrid_haplotype2.chr10.122000000-124000000.out -eo ukbb_hybrid_haplotype2.chr10.122000000-124000000.err -R "rusage[mem=30000]" -q large_mem "sh /home/wletsou/scripts/ukbb_hybrid_haplotype2.sh ukbb_bca_cases.indiv,ukbb_bca_controls.indiv 10 122000000,124000000 /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/biobank/gds/ukb.bca.hap/ukb.bca.hap.chr10.new.vcf.gz" # requires gds to be converted into vcf, see ukbb_gds2vcf.R if not already done; if running for the first time, use ukb.bca.hap.chr10.vcf.gz as the vcf file, otherwise ukb.bca.hap.chr10.new.vcf.gz

# cd ${DIRECTORY} # move back up

# vcf and haplotypes for Phase 1

# get genotyped SNPs in region
# awk 'BEGIN{OFS="\t"} ($2>=123240311 && $2<=123440311){print $3}' ${DIRECTORY}/ukbb_genotyped_snps/ukbb_bca_cases+ukbb_bca_controls.chr10.123000000-124000000.haplotypes.vcf > ukbb_snp_list.chr10.123240311-123440311.txt

# check that SNPs have good imputation quality
# awk 'NR==FNR{snp[$1]; next} ($8 in snp && $15==1){print $8}' ukbb_snp_list.chr10.123240311-123440311.txt /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr10.120580000.122580000.qced.anno.info > ukbb_snp_list.chr10.123240311-123440311.tmp
# test -f ukbb_snp_list.chr10.123240311-123440311.tmp && mv ukbb_snp_list.chr10.123240311-123440311.tmp ukbb_snp_list.chr10.123240311-123440311.txt

# subjects lists
# awk 'BEGIN{OFS="\t"} ($3==1){print $1,$1}' /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_cases.indiv
# awk 'BEGIN{OFS="\t"} ($3==0){print $1,$1}' /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_controls.indiv

# vcf file for UKBB
# bsub -P SJLIFE -J ukbb_topmed_extract.chr10.123240311-123440311 -oo ukbb_topmed_extract.chr10.123240311-123440311.out -eo ukbb_topmed_extract.chr10.123240311-123440311.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_topmed_extract.sh ukbb_snp_list.chr10.123240311-123440311.txt /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr10.120580000.122580000.qced.anno.info 10 123240311,123440311 /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/ukb.bca.topmedr2.chr10.120580000.122580000.qced.hg38.vcf.gz"

# vcf file for dbGaP
# bsub -P SJLIFE -J dbgap_haplotype_merge.chr10.123240311-123440311 -oo dbgap_haplotype_merge.chr10.123240311-123440311.out -eo dbgap_haplotype_merge.chr10.123240311-123440311.err -R "rusage[mem=1000]" "sh ${WORKING_DIR}/dbgap_haplotype_merge.sh dbgap28544_cases,dbgap28544_controls "/research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt" "/research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/chr10.qced.info" ukbb_snp_list.chr10.123240311-123440311.txt 10 123240311,123440311 /research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/vcfs/dbgap28544.topmedr2.cleaned.hg38.chr10.vcf.gz" # will take 2-3 hrs to index chrXX vcf if running for the first time

# check that there are no overlapping samples
# module load bcftools/1.10.2

# awk 'NR==FNR{seen[$1]; next} ($1 in seen){print $0}' <(bcftools query -l dbgap28544_cases+dbgap28544_controls.hg19_chr10.123240311-123440311.hg38.vcf.gz) <(bcftools query -l ukbb.topmed.hg19_chr10.123240311-123440311.hg38.vcf.gz) # dbGaP samples in UKBB, if any

# awk 'NR==FNR{seen[$1]; next} ($1 in seen){print $0}' <(bcftools query -l ukbb.topmed.hg19_chr10.123240311-123440311.hg38.vcf.gz) <(bcftools query -l dbgap28544_cases+dbgap28544_controls.hg19_chr10.123240311-123440311.hg38.vcf.gz) # UKBB samples in dbGaP, if any

# check for missing SNPs
# awk 'NR==FNR{seen[$3]; next} ($3 in seen==0){print $3}' <(bcftools view --no-header dbgap28544_cases+dbgap28544_controls.hg19_chr10.123240311-123440311.hg38.vcf.gz) <(bcftools view --no-header ukbb.topmed.hg19_chr10.123240311-123440311.hg38.vcf.gz) # SNPs in UKBB not in DRIVE, if any

# awk 'NR==FNR{seen[$3]; next} ($3 in seen==0){print $3}' <(bcftools view --no-header ukbb.topmed.hg19_chr10.123240311-123440311.hg38.vcf.gz) <(bcftools view --no-header dbgap28544_cases+dbgap28544_controls.hg19_chr10.123240311-123440311.hg38.vcf.gz) # SNPs in DRIVE not in UKBB, if any

# first round of overlaps
# bsub -P SJLIFE -J ukbb_bca_overlap.chr10.123240311-123440311 -oo ukbb_bca_overlap.chr10.123240311-123440311.out -eo ukbb_bca_overlap.chr10.123240311-123440311.err -R "rusage[mem=256]" "sh ${WORKING_DIR}/ukbb_submission13.sh ukbb_bca_cases,ukbb_bca_controls,dbgap28544_cases,dbgap28544_controls 10 123240311,123440311 \"\" ukbb_snp_list.chr10.123240311-123440311.txt ukbb.topmed.hg19_chr10.123240311-123440311.hg38.vcf.gz,dbgap28544_cases+dbgap28544_controls.hg19_chr10.123240311-123440311.hg38.vcf.gz /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr10.120580000.122580000.qced.anno.info 0.00000001 \"\" \"\" \"\" ${DIRECTORY} ${WORKING_DIR}"

# sample of all cases and controls (for LD calculations)
# awk 'BEGIN{OFS="\t"} {print $0}' ukbb_bca_cases.indiv ukbb_bca_controls.indiv > ukbb_bca_cases+ukbb_bca_controls.indiv
# awk 'NR==FNR{id[$1]; next} (NR!=FNR && FNR==1){print $0} (NR!=FNR && FNR>1 && $1 in id){print $0}' ukbb_bca_cases+ukbb_bca_controls.indiv haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.chr10.123240311-123440311.txt > haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr10.123240311-123440311.txt

# set a threshold for remaining overlaps
# mv Pattern_combined_old.Iteration000.ukbb_bca_cases.chr10.123240311-123440311_2,j.txt Pattern_combined.Iteration000.ukbb_bca_cases.chr10.123240311-123440311_2,j.txt

# awk 'BEGIN{OFS="\t"} NR==FNR{count[$2]=$1; next} ($1 in count){print count[$1],$1}' Pattern_combined.Iteration000.ukbb_bca_cases.chr10.123240311-123440311_2,j.txt <(awk '($5+0<1e-19){print $1}' Pattern_differences.Iteration000.ukbb_bca_cases.chr10.123240311-123440311_2,j.p_vals.txt) > Pattern_combined.Iteration000.ukbb_bca_cases.chr10.123240311-123440311_2,j.tmp # change p-value in $5+0<p-value to get patterns passing a new threshold

# mv Pattern_combined.Iteration000.ukbb_bca_cases.chr10.123240311-123440311_2,j.txt Pattern_combined_old.Iteration000.ukbb_bca_cases.chr10.123240311-123440311_2,j.txt

# mv Pattern_combined.Iteration000.ukbb_bca_cases.chr10.123240311-123440311_2,j.tmp Pattern_combined.Iteration000.ukbb_bca_cases.chr10.123240311-123440311_2,j.txt

# remove old files from previous run (can skip if running iterations for the first time)
# for file in *Iteration00[1-9]*; do rm -r $file; done # keeps Iteration000 and removes Iterations 001-009; if Iterations 010 and above exist, need to rewrite code to remove them

# for file in fisher_exact*; do rm -r $file; done

# for file in Core*; do rm $file; done

# for file in Closed*; do rm $file; done

# Begin iterations:
# bsub -P SJLIFE -J pattern_overlap_iterate3.ukbb_bca_cases.chr10.123240311-123440311 -oo ${DIRECTORY}/pattern_overlap_iterate3.ukbb_bca_cases.chr10.123240311-123440311.out -eo ${DIRECTORY}/pattern_overlap_iterate3.ukbb_bca_cases.chr10.123240311-123440311.err -R "rusage[mem=512]" "sh ${WORKING_DIR}/pattern_overlap_iterate3.sh haplotype_estimates.ukbb_bca_cases.chr10.123240311-123440311.txt,haplotype_estimates.ukbb_bca_controls.chr10.123240311-123440311.txt 10 123240311,123440311 100 50 2 2,j \"\" \"ukbb_bca_cases\" ${DIRECTORY} ${WORKING_DIR}"

# details on iterations, make plots of total and closed patterns
# sh ${WORKING_DIR}/ukbb_pattern_count_summary.sh ukbb_bca_cases 10 123240311,123440311 2,j

# Haplotypes with p < 1e-XX
# test -f fisher_exact.ukbb_bca_cases.Results.txt && rm fisher_exact.ukbb_bca_cases.Results.txt; touch fisher_exact.ukbb_bca_cases.Results.txt; for file in fisher_exact.ukbb_bca_cases.Iteration*.txt; do awk '($6+0<1e-55){print $0}' $file >> fisher_exact.ukbb_bca_cases.Results.txt; done;

# translated significant haplotypes from column numbers to rsids
# bsub -P SJLIFE -J ukbb_haplotype_translate2_sub -oo ukbb_haplotype_translate2_sub.out -eo ukbb_haplotype_translate2_sub.err -R "rusage[mem=256]" -q standard "sh ${WORKING_DIR}/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr10.123240311-123440311.txt fisher_exact.ukbb_bca_cases.Results.txt 50"

# counts of translated haplotypes on each chromosome
# bsub -P SJLIFE -J ukbb_haplotype_model9_sub -oo ukbb_haplotype_model9_sub.out -eo ukbb_haplotype_model9_sub.err -R "rusage[mem=256]" "sh ${WORKING_DIR}/ukbb_haplotype_model9_sub.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr10.123240311-123440311.txt fisher_exact.ukbb_bca_cases.Results.translated.txt \"\" \"\" \"\" 10 123240311,123440311 50 \"\" \"\" /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1 /home/wletsou/scripts"

# add patterns to model
# bsub -P SJLIFE -J ukbb_haplotype_model9_iterate -oo ukbb_haplotype_model9_iterate.out -eo ukbb_haplotype_model9_iterate.err -R "rusage[mem=1000]" "sh ${WORKING_DIR}/ukbb_haplotype_model9_iterate.sh allele_counts.txt /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt \"\" \"\" 1e-5 Significant_patterns.txt ${DIRECTORY} ${WORKING_DIR}"

# awk 'BEGIN{n=split("rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs17542768_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)_[A-Z]=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } ($3 in haplotype_snps){found+=1; printf "%s[%s]%s",$3,(haplotype_snps[$3]==1?$5:$4),(found<n?",":"\n")}' ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr10.123240311-123440311.haplotypes.hg38.vcf # h1 alleles
# rs3135795[G],rs116941692[C],rs755793[A],rs56226109[G],rs2981579[A],rs17542768[A],rs45631549[C],rs45631614[G],rs45631630[G],rs17542858[G],rs118181559[G]

# reduction
# sh ${WORKING_DIR}/haplotype_rpart.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr10.123240311-123440311.txt rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs17542768_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0 "" 1 # h1
#  reduced = rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0

# awk 'BEGIN{n=split("rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)_[A-Z]=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } ($3 in haplotype_snps){found+=1; printf "%s[%s]%s",$3,(haplotype_snps[$3]==1?$5:$4),(found<n?",":"\n")}' ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr10.123240311-123440311.haplotypes.hg38.vcf # h1.reduced alleles
# rs3135795[G],rs116941692[C],rs755793[A],rs56226109[G],rs2981579[A],rs45631549[C],rs45631614[G],rs45631630[G],rs17542858[G],rs118181559[G]

# check that h1 and h1.reduced have the same effect sizes
# sh ${WORKING_DIR}/ukbb_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr10.123240311-123440311.txt "rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs17542768_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0" "" "" "" "" h1 ${DIRECTORY} ${WORKING_DIR}
# LRT p value = 1.020510e-65 on 1 degree(s) of freedom with HR = 1.295197e+00 (1.257649e+00 to 1.333865e+00), frequency = 3.797905e-01 (137510) = 4.402397e-01/3.766241e-01 (7934/129576)

# sh ${WORKING_DIR}/ukbb_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr10.123240311-123440311.txt "rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0" "" "" "" "" h1.reduced ${DIRECTORY} ${WORKING_DIR}
# LRT p value = 1.020510e-65 on 1 degree(s) of freedom with HR = 1.295197e+00 (1.257649e+00 to 1.333865e+00), frequency = 3.797905e-01 (137510) = 4.402397e-01/3.766241e-01 (7934/129576)

# effect of GWAS hit

# mkdir ukbb_imputed_snps
# cp ukbb_bca_cases.indiv ukbb_imputed_snps/ukbb_bca_cases.indiv
# cp ukbb_bca_controls.indiv ukbb_imputed_snps/ukbb_bca_controls.indiv
# cd ukbb_imputed_snps

# extract imputed SNPs in range (h1.full-length plus GWAS hit)
# bsub -P SJLIFE -J ukbb_haplotype_extract3.chr10.123240311-123440311 -oo ukbb_haplotype_extract3.chr10.123240311-123440311.out -eo ukbb_haplotype_extract3.chr10.123240311-123440311.err -R "rusage[mem=10000]" -q standard "sh ${WORKING_DIR}/ukbb_haplotype_extract3.sh ukbb_bca_cases.indiv,ukbb_bca_controls.indiv \"rs3135795,rs116941692,rs755793,rs56226109,rs2981579,rs17542768,rs45631549,rs45631614,rs45631630,rs17542858,rs118181559,rs2981578\" /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr10.120580000.122580000.qced.anno.info 10 123240311,123440311 /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/ukb.bca.topmedr2.chr10.120580000.122580000.qced.hg38.vcf.gz /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/ukbb_imputed_snps /home/wletsou/scripts"

# sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr10.123240311-123440311.txt "rs2981578_T=0" "" "" "" "" rs2981578 /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/ukbb_imputed_snps /home/wletsou/scripts
# LRT p value = 1.678709e-51 on 1 degree(s) of freedom with HR = 1.252367e+00 (1.216299e+00 to 1.289504e+00), frequency = 4.636781e-01 (167883) = 5.178116e-01/4.608424e-01 (9332/158551)

# LD with GWAS hit

# bsub -P SJLIFE -J ukbb_haplotype_haplotype_LD2 -oo ukbb_haplotype_haplotype_LD2.out -eo ukbb_haplotype_haplotype_LD2.err -R "rusage[mem=512]" "sh ${WORKING_DIR}/ukbb_haplotype_haplotype_LD2.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr10.123240311-123440311.txt \"rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0\" \"rs2981578_T=0\" ${DIRECTORY}" # h1 (reduced) vs GWAS hit
# 0.379791	0.463678	0.377954
# rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0	rs2981578_T=0	0.69558	0.990983

# plot SNPs

# echo "rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0" | sed 's/_[A-Z]=[0-9]//g' # snps in h1, no alleles
# convert .vcf to .bed at SNPs in h1 (reduced)
# sh /home/wletsou/scripts/vcf2bed.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr10.123240311-123440311.haplotypes.hg38.vcf "rs3135795,rs116941692,rs755793,rs56226109,rs2981579,rs45631549,rs45631614,rs45631630,rs17542858,rs118181559" "" h1.bed /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1 /home/wletsou/scripts # GRCh38 bed file

# cd ${DIRECTORY}
# mkdir phase2.1
# cd phase2.1

# get h1 and list of SNPs in remainder of TAD
# HAPLOTYPE_SNPS=$(echo "rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0" | sed 's/_[A-Z]=[0-9]//g') # snps in h1, no alleles

# awk 'BEGIN{OFS="\t"} ( ($2<123240311 && $2>=122850000) || "'$HAPLOTYPE_SNPS'" ~ $3){print $3}' ${DIRECTORY}/ukbb_genotyped_snps/ukbb_bca_cases+ukbb_bca_controls.chr10.123000000-124000000.haplotypes.vcf > TAD_snps.chr10.123440312-123675000.txt # http://3dgenome.fsm.northwestern.edu/view.php?method=Hi-C&species=human&assembly=hg19&source=inside&tissue=HMEC&type=Lieberman-raw&resolution=5&c_url=&transfer=&chr=chr10&start=122500000&end=124000000&sessionID=&browser=ucsc

# check that SNPs have good "imputation" quality
# awk 'NR==FNR{snp[$1]; next} ($8 in snp && $15==1){print $8}' TAD_snps.chr10.123440312-123675000.txt /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr10.120580000.122580000.qced.anno.info > TAD_snps.chr10.123440312-123675000.tmp
# test -f TAD_snps.chr10.123440312-123675000.tmp && mv TAD_snps.chr10.123440312-123675000.tmp TAD_snps.chr10.123440312-123675000.txt

# awk 'BEGIN{OFS="\t"} ($3==1){print $1,$1}' /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_cases.indiv
# awk 'BEGIN{OFS="\t"} ($3==0){print $1,$1}' /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_controls.indiv

# vcf file for UKBB
# bsub -P SJLIFE -J ukbb_topmed_extract.TAD_snps.chr10.123440312-123675000 -oo ukbb_topmed_extract.TAD_snps.chr10.123440312-123675000.out -eo ukbb_topmed_extract.TAD_snps.chr10.123440312-123675000.err -R "rusage[mem=1000]" "sh ${WORKING_DIR}/ukbb_topmed_extract.sh TAD_snps.chr10.123440312-123675000.txt /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr10.120580000.122580000.qced.anno.info 10 123440312,123675000 /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/ukb.bca.topmedr2.chr10.120580000.122580000.qced.hg38.vcf.gz ${DIRECTORY}/phase2 ${WORKING_DIR}"

# vcf file for dbGaP
# bsub -P SJLIFE -J dbgap_haplotype_merge.chr10.123440312-123675000 -oo dbgap_haplotype_merge.chr10.123440312-123675000.out -eo dbgap_haplotype_merge.chr10.123440312-123675000.err -R "rusage[mem=1000]" "sh ${WORKING_DIR}/dbgap_haplotype_merge.sh dbgap28544_cases,dbgap28544_controls "/research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt" "/research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/chr10.qced.info" TAD_snps.chr10.123440312-123675000.txt 10 123440312,123675000 /research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/vcfs/dbgap28544.topmedr2.cleaned.hg38.chr10.vcf.gz ${DIRECTORY}/phase2 ${WORKING_DIR}"

# check that there are no overlapping samples
# module load bcftools/1.10.2

# awk 'NR==FNR{seen[$1]; next} ($1 in seen){print $0}' <(bcftools query -l dbgap28544_cases+dbgap28544_controls.hg19_chr10.123440312-123675000.hg38.vcf.gz) <(bcftools query -l ukbb.topmed.hg19_chr10.123440312-123675000.hg38.vcf.gz) # dbGaP samples in UKBB, if any

# awk 'NR==FNR{seen[$1]; next} ($1 in seen){print $0}' <(bcftools query -l ukbb.topmed.hg19_chr10.123440312-123675000.hg38.vcf.gz) <(bcftools query -l dbgap28544_cases+dbgap28544_controls.hg19_chr10.123440312-123675000.hg38.vcf.gz) # UKBB samples in dbGaP, if any

# check for missing SNPs
# awk 'NR==FNR{seen[$3]; next} ($3 in seen==0){print $3}' <(bcftools view --no-header dbgap28544_cases+dbgap28544_controls.hg19_chr10.123440312-123675000.hg38.vcf.gz) <(bcftools view --no-header ukbb.topmed.hg19_chr10.123440312-123675000.hg38.vcf.gz) # SNPs in UKBB not in DRIVE, if any

# awk 'NR==FNR{seen[$3]; next} ($3 in seen==0){print $3}' <(bcftools view --no-header ukbb.topmed.hg19_chr10.123440312-123675000.hg38.vcf.gz) <(bcftools view --no-header dbgap28544_cases+dbgap28544_controls.hg19_chr10.123440312-123675000.hg38.vcf.gz) # SNPs in DRIVE not in UKBB, if any

# first round of overlaps (of h1-carrying chromosomes among cases)
# bsub -P SJLIFE -J ukbb_bca_overlap.chr10.123440312-123675000 -oo ukbb_bca_overlap.chr10.123440312-123675000.out -eo ukbb_bca_overlap.chr10.123440312-123675000.err -R "rusage[mem=256]" "sh ${WORKING_DIR}/ukbb_submission13.sh ukbb_bca_cases,ukbb_bca_controls,dbgap28544_cases,dbgap28544_controls 10 123440312,123675000 "rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0" TAD_snps.chr10.123440312-123675000.txt ukbb.topmed.hg19_chr10.123440312-123675000.hg38.vcf.gz,dbgap28544_cases+dbgap28544_controls.hg19_chr10.123440312-123675000.hg38.vcf.gz /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr10.120580000.122580000.qced.anno.info 0.0001 \"\" \"\" \"\" ${DIRECTORY}/phase2 ${WORKING_DIR}"
