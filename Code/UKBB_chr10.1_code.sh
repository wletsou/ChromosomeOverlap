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
# mkdir phase2
# cd phase2

# get h1 and list of SNPs in remainder of TAD
# HAPLOTYPE_SNPS=$(echo "rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0" | sed 's/_[A-Z]=[0-9]//g') # snps in h1, no alleles

# awk 'BEGIN{OFS="\t"} ( ($2>123440311 && $2<=123675000) || "'$HAPLOTYPE_SNPS'" ~ $3){print $3}' ${DIRECTORY}/ukbb_genotyped_snps/ukbb_bca_cases+ukbb_bca_controls.chr10.123000000-124000000.haplotypes.vcf > TAD_snps.chr10.123440312-123675000.txt # http://3dgenome.fsm.northwestern.edu/view.php?method=Hi-C&species=human&assembly=hg19&source=inside&tissue=HMEC&type=Lieberman-raw&resolution=5&c_url=&transfer=&chr=chr10&start=122500000&end=124000000&sessionID=&browser=ucsc

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

# set a threshold for remaining overlaps
# mv Pattern_combined_old.Iteration000.ukbb_bca_cases.chr10.123440312-123675000_2,j.txt Pattern_combined.Iteration000.ukbb_bca_cases.chr10.123440312-123675000_2,j.txt

# pick patterns with p<XX by changing awk '($5+0)<XX{print $0}...'
# awk 'BEGIN{OFS="\t"} NR==FNR{count[$2]=$1; next} ($1 in count){print count[$1],$1}' Pattern_combined.Iteration000.ukbb_bca_cases.chr10.123440312-123675000_2,j.txt <(awk '($5+0<1.7e-4){print $0}' Pattern_differences.Iteration000.ukbb_bca_cases.chr10.123440312-123675000_2,j.p_vals.txt | awk 'BEGIN{array[0]} {delete array; n=split($1,array,","); if (n<snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]+0 || snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]+0==0) {snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]=n; seen[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]=$1} } END{for (i in seen) {print seen[i]} }') > Pattern_combined.Iteration000.ukbb_bca_cases.chr10.123440312-123675000_2,j.tmp

# mv Pattern_combined.Iteration000.ukbb_bca_cases.chr10.123440312-123675000_2,j.txt Pattern_combined_old.Iteration000.ukbb_bca_cases.chr10.123440312-123675000_2,j.txt

# mv Pattern_combined.Iteration000.ukbb_bca_cases.chr10.123440312-123675000_2,j.tmp Pattern_combined.Iteration000.ukbb_bca_cases.chr10.123440312-123675000_2,j.txt

# remove old files from previous run (can skip if running iterations for the first time)
# for file in *Iteration00[1-9]*; do rm -r $file; done # keeps Iteration000 and removes Iterations 001-009; if Iterations 010 and above exist, need to rewrite code to remove them

# for file in fisher_exact*; do rm -r $file; done

# for file in Core*; do rm $file; done

# for file in Closed*; do rm $file; done

# Begin iterations:
# bsub -P SJLIFE -J pattern_overlap_iterate3.ukbb_bca_cases.chr10.123440312-123675000 -oo pattern_overlap_iterate3.ukbb_bca_cases.chr10.123440312-123675000.out -eo pattern_overlap_iterate3.ukbb_bca_cases.chr10.123440312-123675000.err -R "rusage[mem=512]" -R "status=='ok'" "sh /home/wletsou/scripts/pattern_overlap_iterate3.sh haplotype_estimates.ukbb_bca_cases.chr10.123440312-123675000.subset.txt,haplotype_estimates.ukbb_bca_controls.chr10.123440312-123675000.subset.txt 10 123440312,123675000 100 50 2 2,j \"\" \"ukbb_bca_cases\" /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/phase2 /home/wletsou/scripts"

# sh /home/wletsou/scripts/ukbb_pattern_count_summary.sh ukbb_bca_cases 10 123440312,123675000 2,j

# combined haplotype file for cases and controls
# awk 'BEGIN{OFS="\t"} {print $0}' ukbb_bca_cases.indiv ukbb_bca_controls.indiv > ukbb_bca_cases+ukbb_bca_controls.indiv
# awk 'NR==FNR{id[$1]; next} (NR!=FNR && FNR==1){print $0} (NR!=FNR && FNR>1 && $1 in id){print $0}' ukbb_bca_cases+ukbb_bca_controls.indiv haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.chr10.123440312-123675000.txt > haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr10.123440312-123675000.txt

# forward selection
# translated in parallel, p <
# test -f fisher_exact.ukbb_bca_cases.Results.txt && rm fisher_exact.ukbb_bca_cases.Results.txt; touch fisher_exact.ukbb_bca_cases.Results.txt; for file in fisher_exact.ukbb_bca_cases.Iteration*.txt; do awk '($6+0<1e-2){print $0}' $file >> fisher_exact.ukbb_bca_cases.Results.txt; done;

# remove longer patterns with the same frequency, OR, and p-value as a shorter pattern
# awk 'BEGIN{array[0]} {delete array; n=split($1,array,","); if (n<snps[sprintf("%s\t%s\t%s\t%s",$3,$4,$5,$6)]+0 || snps[sprintf("%s\t%s\t%s\t%s",$3,$4,$5,$6)]+0==0) {snps[sprintf("%s\t%s\t%s\t%s",$3,$4,$5,$6)]=n; seen[sprintf("%s\t%s\t%s\t%s",$3,$4,$5,$6)]=$0} }  END{for (i in seen) {print seen[i]} }' fisher_exact.ukbb_bca_cases.Results.txt > fisher_exact.ukbb_bca_cases.Results.tmp; test -f fisher_exact.ukbb_bca_cases.Results.tmp && mv fisher_exact.ukbb_bca_cases.Results.tmp fisher_exact.ukbb_bca_cases.Results.txt

# bsub -P SJLIFE -J ukbb_haplotype_translate2_sub -oo ukbb_haplotype_translate2_sub.out -eo ukbb_haplotype_translate2_sub.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr10.123440312-123675000.subset.txt fisher_exact.ukbb_bca_cases.Results.txt 50"

# bsub -P SJLIFE -J ukbb_haplotype_model9_sub -oo ukbb_haplotype_model9_sub.out -eo ukbb_haplotype_model9_sub.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr10.123440312-123675000.txt fisher_exact.ukbb_bca_cases.Results.translated.txt 1 \"rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0\" 1 10 123440312,123675000 50 \"\" \"\" /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/phase2_1.8e-4 /home/wletsou/scripts"

# bsub -P SJLIFE -J ukbb_haplotype_model9_iterate -oo ukbb_haplotype_model9_iterate.out -eo ukbb_haplotype_model9_iterate.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.sh allele_counts.txt /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt 1 1 1e-5 Significant_patterns.txt /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/phase2_1.8e-4 /home/wletsou/scripts"

# HAPLOTYPE_LIST=$(awk '(NR>1 && $6+0<1e-5 && $3 > 1){printf "%s%s",(found>0?":":""),$2; found+=1} END{printf "\n"}' Significant_patterns.txt)

# sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr10.123440312-123675000.txt "rs113940757_A=0,rs7924234_C=0,rs117256623_C=0,rs1696820_T=0,rs12780622_A=0,rs11200087_G=0,rs17102476_C=0,rs78516665_A=0,rs11200090_C=0,rs117867446_A=0,rs74158424_C=0,rs61874153_A=0,rs2935706_C=0,rs1693669_A=1,rs12413911_A=0,rs7096566_T=0,rs74158435_G=0,rs72836300_G=0,rs1693675_G=1,rs2420957_G=1,rs2420956_G=1,rs7893651_C=1,rs80195214_C=0,rs1874662_A=1,rs113793136_G=0,rs10749431_G=1,rs79652084_A=0,rs1909009_T=1,rs11594831_A=0,rs116938619_T=0,rs4751859_G=1,rs11200194_A=1,rs4444014_C=1,rs10510102_C=0,rs112376407_C=0,rs11200204_A=0,rs72838067_C=0,rs10732824_G=1,rs11200231_C=0,rs148095496_A=0" "1" "rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0" "1" "" h2 /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/phase2_1.8e-4 /home/wletsou/scripts
# Rscript /home/wletsou/scripts/ukbb_haplotype_model5.R file=/research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/phase2_1.8e-4/h2.new_allele_counts.txt groups=1,2 verbose=1 n=1 n_covs=10
# [1] "X00" "X10" "X11"
# [1] "DT[,h0:=X00]"
# [1] "DT[,h1.2:=X10+X11]"
# [1] "DT[,h1:=X10]"
# [1] "DT[,h2:=X11]"
#
#    affected   h1.2     h1 h2
# 1:        0 129576 129495 81
# 2:        1   7934   7914 20
#
#    affected      h1.2        h1           h2
# 1:        0 0.3766241 0.3763886 0.0002354336
# 2:        1 0.4402397 0.4391300 0.0011097547
#
#      h1.2     h1  h2
# 1: 137510 137409 101
#
#         h1.2        h1           h2
# 1: 0.3797905 0.3795116 0.0002789531
#
# Surv(age,affected) ~ h1 + h2 + pc01 + pc02 + pc03 + pc04 + pc05 + pc06 + pc07 + pc08 + pc09 + pc10
# Call:
# coxph(formula = eval(parse(text = Y)), data = X)
#
#   n= 181034, number of events= 9011
#
#            coef  exp(coef)   se(coef)      z Pr(>|z|)
# h1    0.2566032  1.2925322  0.0150207 17.083  < 2e-16 ***
# h2    1.6410258  5.1604603  0.2239832  7.327 2.36e-13 ***
# pc01 -0.0009421  0.9990583  0.0069138 -0.136   0.8916
# pc02 -0.0142253  0.9858754  0.0071597 -1.987   0.0469 *
# pc03  0.0068355  1.0068589  0.0068966  0.991   0.3216
# pc04  0.0038131  1.0038203  0.0051764  0.737   0.4613
# pc05 -0.0013079  0.9986929  0.0022694 -0.576   0.5644
# pc06  0.0067009  1.0067234  0.0065758  1.019   0.3082
# pc07  0.0044960  1.0045061  0.0059061  0.761   0.4465
# pc08 -0.0012275  0.9987733  0.0058811 -0.209   0.8347
# pc09 -0.0029074  0.9970968  0.0022830 -1.274   0.2028
# pc10  0.0062988  1.0063187  0.0050961  1.236   0.2165
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#      exp(coef) exp(-coef) lower .95 upper .95
# h1      1.2925     0.7737    1.2550    1.3311
# h2      5.1605     0.1938    3.3269    8.0047
# pc01    0.9991     1.0009    0.9856    1.0127
# pc02    0.9859     1.0143    0.9721    0.9998
# pc03    1.0069     0.9932    0.9933    1.0206
# pc04    1.0038     0.9962    0.9937    1.0141
# pc05    0.9987     1.0013    0.9943    1.0031
# pc06    1.0067     0.9933    0.9938    1.0198
# pc07    1.0045     0.9955    0.9929    1.0162
# pc08    0.9988     1.0012    0.9873    1.0104
# pc09    0.9971     1.0029    0.9926    1.0016
# pc10    1.0063     0.9937    0.9963    1.0164
#
# Concordance= 0.552  (se = 0.003 )
# Likelihood ratio test= 331  on 12 df,   p=<2e-16
# Wald test            = 351.6  on 12 df,   p=<2e-16
# Score (logrank) test = 363.6  on 12 df,   p=<2e-16
#
# LRT p value = 4.789725e-07 on 1 degree(s) of freedom with HR = 5.160460e+00 (3.326853e+00 to 8.004668e+00), frequency = 2.789531e-04 (101) = 1.109755e-03/2.354336e-04 (20/81)

# backward selection
# module load R/3.6.1
# Rscript /home/wletsou/scripts/ukbb_backward_selection2.R file=h2.new_allele_counts.txt n=1 groups=1,2 threshold=0.05 n_covs=10
# remaining_vars         pval
# 1:             h2 4.789725e-07

# reduction
# sh /home/wletsou/scripts/haplotype_rpart.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr10.123440312-123675000.txt "rs113940757_A=0,rs7924234_C=0,rs117256623_C=0,rs1696820_T=0,rs12780622_A=0,rs11200087_G=0,rs17102476_C=0,rs78516665_A=0,rs11200090_C=0,rs117867446_A=0,rs74158424_C=0,rs61874153_A=0,rs2935706_C=0,rs1693669_A=1,rs12413911_A=0,rs7096566_T=0,rs74158435_G=0,rs72836300_G=0,rs1693675_G=1,rs2420957_G=1,rs2420956_G=1,rs7893651_C=1,rs80195214_C=0,rs1874662_A=1,rs113793136_G=0,rs10749431_G=1,rs79652084_A=0,rs1909009_T=1,rs11594831_A=0,rs116938619_T=0,rs4751859_G=1,rs11200194_A=1,rs4444014_C=1,rs10510102_C=0,rs112376407_C=0,rs11200204_A=0,rs72838067_C=0,rs10732824_G=1,rs11200231_C=0,rs148095496_A=0" "rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0" 1 # h2
#  reduced = rs113940757_A=0,rs7924234_C=0,rs1696820_T=0,rs11200087_G=0,rs17102476_C=0,rs11200090_C=0,rs61874153_A=0,rs2935706_C=0,rs1693675_G=1,rs2420957_G=1,rs80195214_C=0,rs113793136_G=0,rs11200194_A=1,rs11200231_C=0

# dbGaP replication

# mkdir dbgap28544_imputed_snps
# cd dbgap28544_imputed_snps

# extract imputed SNPs in range (h1.full-length plus GWAS hit)
# bsub -P SJLIFE -J ukbb_haplotype_extract3.chr10.123240311-123440311 -oo ukbb_haplotype_extract3.chr10.123240311-123440311.out -eo ukbb_haplotype_extract3.chr10.123240311-123440311.err -R "rusage[mem=10000]" -q large_mem "sh /home/wletsou/scripts/ukbb_haplotype_extract3.sh dbgap28544_cases.indiv,dbgap28544_controls.indiv \"rs3135795,rs116941692,rs755793,rs56226109,rs2981579,rs17542768,rs45631549,rs45631614,rs45631630,rs17542858,rs118181559,rs2981578\" /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr10.120580000.122580000.qced.anno.info 10 123240311,123440311 /research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/vcfs/dbgap28544.topmedr2.cleaned.hg38.chr10.vcf.gz /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/dbgap28544_imputed_snps /home/wletsou/scripts"

# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr10.123240311-123440311.txt "rs2981578_T=0" "" "" "" "" rs2981578 /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/dbgap28544_imputed_snps /home/wletsou/scripts
# LRT p value = 1.655436e-58 on 1 degree(s) of freedom with OR = 1.217399e+00 (1.188571e+00 to 1.246927e+00), frequency = 5.018339e-01 (55549) = 5.239489e-01/4.755360e-01 (31504/24045)

# LD with GWAS hit

# bsub -P SJLIFE -J ukbb_haplotype_haplotype_LD2 -oo ukbb_haplotype_haplotype_LD2.out -eo ukbb_haplotype_haplotype_LD2.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ukbb_haplotype_haplotype_LD2.sh haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr10.123240311-123440311.txt \"rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0\" \"rs2981578_T=0\" /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/dbgap28544_imputed_snps" # h1 (reduced) vs GWAS hit
# 0.410626	0.501834	0.407789
# rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0	rs2981578_T=0	0.672574	0.986133

# cd /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/phase2_1.8e-4

# combined haplotype file for cases and controls
# awk 'BEGIN{OFS="\t"} {print $0}' dbgap28544_cases.indiv dbgap28544_controls.indiv > dbgap28544_cases+dbgap28544_controls.indiv
# awk 'NR==FNR{id[$1]; next} (NR!=FNR && FNR==1){print $0} (NR!=FNR && FNR>1 && $1 in id){print $0}' dbgap28544_cases+dbgap28544_controls.indiv haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.chr10.123440312-123675000.txt > haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr10.123440312-123675000.txt

# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr10.123440312-123675000.txt "rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0" "" "" "" "" dbgap28544.h1.reduced /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/phase2_1.8e-4 /home/wletsou/scripts # h1 reduced
# LRT p value = 1.465795e-74 on 1 degree(s) of freedom with OR = 1.253929e+00 (1.223743e+00 to 1.284859e+00), frequency = 4.106259e-01 (45453) = 4.354045e-01/3.811605e-01 (26180/19273)

# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr10.123440312-123675000.txt "rs113940757_A=0,rs7924234_C=0,rs117256623_C=0,rs1696820_T=0,rs12780622_A=0,rs11200087_G=0,rs17102476_C=0,rs78516665_A=0,rs11200090_C=0,rs117867446_A=0,rs74158424_C=0,rs61874153_A=0,rs2935706_C=0,rs1693669_A=1,rs12413911_A=0,rs7096566_T=0,rs74158435_G=0,rs72836300_G=0,rs1693675_G=1,rs2420957_G=1,rs2420956_G=1,rs7893651_C=1,rs80195214_C=0,rs1874662_A=1,rs113793136_G=0,rs10749431_G=1,rs79652084_A=0,rs1909009_T=1,rs11594831_A=0,rs116938619_T=0,rs4751859_G=1,rs11200194_A=1,rs4444014_C=1,rs10510102_C=0,rs112376407_C=0,rs11200204_A=0,rs72838067_C=0,rs10732824_G=1,rs11200231_C=0,rs148095496_A=0" "1" "rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0" "1" "" dbgap28544.h2.full-length /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/phase2_1.8e-4 /home/wletsou/scripts # h2 full-length
# LRT p value = 3.767744e-01 on 1 degree(s) of freedom with OR = 1.655345e+00 (8.881505e-01 to 3.085252e+00), frequency = 3.884653e-04 (43) = 4.490420e-04/3.164307e-04 (27/16)

# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr10.123440312-123675000.txt "rs113940757_A=0,rs7924234_C=0,rs1696820_T=0,rs11200087_G=0,rs17102476_C=0,rs11200090_C=0,rs61874153_A=0,rs2935706_C=0,rs1693675_G=1,rs2420957_G=1,rs80195214_C=0,rs113793136_G=0,rs11200194_A=1,rs11200231_C=0" "1" "rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0" "1" "" dbgap28544.h2.reduced /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/phase2_1.8e-4 /home/wletsou/scripts # h2 reduced
# LRT p value = 6.535891e-01 on 1 degree(s) of freedom with OR = 1.433375e+00 (7.970481e-01 to 2.577717e+00), frequency = 4.246016e-04 (47) = 4.656732e-04/3.757614e-04 (28/19)

# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr10.123440312-123675000.txt "rs113940757_A=0,rs1909004_A=1,rs117256623_C=0,rs7900009_T=0,rs73358541_T=0,rs1696820_T=0,rs12780622_A=0,rs11200087_G=0,rs17102476_C=0,rs78516665_A=0,rs117867446_A=0,rs74158424_C=0,rs2935706_C=0,rs4751849_G=0,rs1693669_A=1,rs11200120_C=0,rs12413911_A=0,rs7897273_T=0,rs7096566_T=0,rs80279826_G=0,rs74158435_G=0,rs72836300_G=0,rs1693675_G=1,rs2420957_G=1,rs2420956_G=1,rs7893651_C=1,rs80195214_C=0,rs1874662_A=1,rs113793136_G=0,rs10749431_G=1,rs79652084_A=0,rs1909009_T=1,rs9651480_G=1,rs11594831_A=0,rs116938619_T=0,rs4751859_G=1,rs11200194_A=0,rs4444014_C=0,rs10510102_C=1,rs112376407_C=0,rs11200204_A=0,rs72838067_C=0,rs10732824_G=1,rs11200231_C=0,rs148095496_A=0" "1" "rs3135795_T=0,rs116941692_T=0,rs755793_G=0,rs56226109_A=0,rs2981579_G=0,rs45631549_T=0,rs45631614_C=0,rs45631630_A=0,rs17542858_A=0,rs118181559_T=0" "1" "" dbgap28544.h3.full-length /research/projects/yasuigrp/EpiGenetics/common/drive4william/UKBB_chr10.1/phase2_1.8e-4 /home/wletsou/scripts # h3 full-length
# LRT p value = 8.650630e-01 on 1 degree(s) of freedom with OR = 1.350209e+00 (5.740069e-01 to 3.176034e+00), frequency = 1.987497e-04 (22) = 2.162054e-04/1.779922e-04 (13/9)
