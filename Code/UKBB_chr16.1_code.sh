# mkdir /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1
# cd /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1

# extract UKBB-genotyped SNPs in entire TAD

# mkdir /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/ukbb_genotyped_snps
# cd ukbb_genotyped_snps

# awk 'BEGIN{OFS="\t"} ($3==1){print $1,$1}' /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_cases.indiv
# awk 'BEGIN{OFS="\t"} ($3==0){print $1,$1}' /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_controls.indiv

# bsub -P SJLIFE -J ukbb_hybrid_haplotype2.chr16.52100000-53100000 -oo ukbb_hybrid_haplotype2.chr16.52100000-53100000.out -eo ukbb_hybrid_haplotype2.chr16.52100000-53100000.err -R "rusage[mem=20000]" "sh /home/wletsou/scripts/ukbb_hybrid_haplotype2.sh ukbb_bca_cases.indiv,ukbb_bca_controls.indiv 16 52100000,53100000 /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/biobank/gds/ukb.bca.hap/ukb.bca.hap.chr16.vcf.gz" # if running for the first time, use ukb.bca.hap.chr16.vcf.gz as the vcf file, otherwise ukb.bca.hap.chr16.new.vcf.gz

# bsub -P SJLIFE -J ukbb_hybrid_haplotype2.chr16.52000000-53200000 -oo ukbb_hybrid_haplotype2.chr16.52000000-53200000.out -eo ukbb_hybrid_haplotype2.chr16.52000000-53200000.err -R "rusage[mem=30000]" "sh /home/wletsou/scripts/ukbb_hybrid_haplotype2.sh ukbb_bca_cases.indiv,ukbb_bca_controls.indiv 16 52000000,53200000 /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/biobank/gds/ukb.bca.hap/ukb.bca.hap.chr16.new.vcf.gz"

# cd /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1

# vcf and haplotypes for Phase 1

# get genotyped SNPs in region
# awk 'BEGIN{OFS="\t"} ($2>=52519188 && $2<=52679188){print $3}' /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/ukbb_genotyped_snps/ukbb_bca_cases+ukbb_bca_controls.chr16.52100000-53100000.haplotypes.vcf > ukbb_snp_list.chr16.52519188-52679188.txt

# check that SNPs have good "imputation" quality
# awk 'NR==FNR{snp[$1]; next} ($8 in snp && $15==1){print $8}' ukbb_snp_list.chr16.52519188-52679188.txt /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr16.51560000.53560000.qced.anno.info > ukbb_snp_list.chr16.52519188-52679188.tmp
# test -f ukbb_snp_list.chr16.52519188-52679188.tmp && mv ukbb_snp_list.chr16.52519188-52679188.tmp ukbb_snp_list.chr16.52519188-52679188.txt

# awk 'BEGIN{OFS="\t"} ($3==1){print $1,$1}' /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_cases.indiv
# awk 'BEGIN{OFS="\t"} ($3==0){print $1,$1}' /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_controls.indiv

# vcf file for UKBB
# bsub -P SJLIFE -J ukbb_topmed_extract.chr16.52519188-52679188 -oo ukbb_topmed_extract.chr16.52519188-52679188.out -eo ukbb_topmed_extract.chr16.52519188-52679188.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_topmed_extract.sh ukbb_snp_list.chr16.52519188-52679188.txt /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr16.51560000.53560000.qced.anno.info 16 52519188,52679188 /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/ukb.bca.topmedr2.chr16.51560000.53560000.qced.hg38.vcf.gz"

# vcf file for dbGaP
# bsub -P SJLIFE -J dbgap_haplotype_merge.chr16.52519188-52679188 -oo dbgap_haplotype_merge.chr16.52519188-52679188.out -eo dbgap_haplotype_merge.chr16.52519188-52679188.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/dbgap_haplotype_merge.sh dbgap28544_cases,dbgap28544_controls "/research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt" "/research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/chr16.qced.info" ukbb_snp_list.chr16.52519188-52679188.txt 16 52519188,52679188 /research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/vcfs/dbgap28544.topmedr2.cleaned.hg38.chr16.vcf.gz"

# check that there are no overlapping samples
# module load bcftools/1.10.2

# awk 'NR==FNR{seen[$1]; next} ($1 in seen){print $0}' <(bcftools query -l dbgap28544_cases+dbgap28544_controls.hg19_chr16.52519188-52679188.hg38.vcf.gz) <(bcftools query -l ukbb.topmed.hg19_chr16.52519188-52679188.hg38.vcf.gz) # dbGaP samples in UKBB, if any

# awk 'NR==FNR{seen[$1]; next} ($1 in seen){print $0}' <(bcftools query -l ukbb.topmed.hg19_chr16.52519188-52679188.hg38.vcf.gz) <(bcftools query -l dbgap28544_cases+dbgap28544_controls.hg19_chr16.52519188-52679188.hg38.vcf.gz) # UKBB samples in dbGaP, if any

# check for missing SNPs
# awk 'NR==FNR{seen[$3]; next} ($3 in seen==0){print $3}' <(bcftools view --no-header dbgap28544_cases+dbgap28544_controls.hg19_chr16.52519188-52679188.hg38.vcf.gz) <(bcftools view --no-header ukbb.topmed.hg19_chr16.52519188-52679188.hg38.vcf.gz) # SNPs in UKBB not in DRIVE, if any

# awk 'NR==FNR{seen[$3]; next} ($3 in seen==0){print $3}' <(bcftools view --no-header ukbb.topmed.hg19_chr16.52519188-52679188.hg38.vcf.gz) <(bcftools view --no-header dbgap28544_cases+dbgap28544_controls.hg19_chr16.52519188-52679188.hg38.vcf.gz) # SNPs in DRIVE not in UKBB, if any

# first round of overlaps
# bsub -P SJLIFE -J ukbb_bca_overlap.chr16.52519188-52679188 -oo ukbb_bca_overlap.chr16.52519188-52679188.out -eo ukbb_bca_overlap.chr16.52519188-52679188.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_submission13.sh ukbb_bca_cases,ukbb_bca_controls,dbgap28544_cases,dbgap28544_controls 16 52519188,52679188 \"\" ukbb_snp_list.chr16.52519188-52679188.txt ukbb.topmed.hg19_chr16.52519188-52679188.hg38.vcf.gz,dbgap28544_cases+dbgap28544_controls.hg19_chr16.52519188-52679188.hg38.vcf.gz /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr16.51560000.53560000.qced.anno.info 0.00000001 \"\" \"\" \"\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1 /home/wletsou/scripts"

# sample of all cases and controls (for LD calculations)
# awk 'BEGIN{OFS="\t"} {print $0}' ukbb_bca_cases.indiv ukbb_bca_controls.indiv > ukbb_bca_cases+ukbb_bca_controls.indiv
# awk 'NR==FNR{id[$1]; next} (NR!=FNR && FNR==1){print $0} (NR!=FNR && FNR>1 && $1 in id){print $0}' ukbb_bca_cases+ukbb_bca_controls.indiv haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.chr16.52519188-52679188.txt > haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52519188-52679188.txt

# set a threshold for remaining overlaps
# mv Pattern_combined_old.Iteration000.ukbb_bca_cases.chr16.52519188-52679188_2,j.txt Pattern_combined.Iteration000.ukbb_bca_cases.chr16.52519188-52679188_2,j.txt

# awk 'BEGIN{OFS="\t"} NR==FNR{count[$2]=$1; next} ($1 in count){print count[$1],$1}' Pattern_combined.Iteration000.ukbb_bca_cases.chr16.52519188-52679188_2,j.txt <(awk '($5+0<1e-16){print $1}' Pattern_differences.Iteration000.ukbb_bca_cases.chr16.52519188-52679188_2,j.p_vals.txt) > Pattern_combined.Iteration000.ukbb_bca_cases.chr16.52519188-52679188_2,j.tmp

# mv Pattern_combined.Iteration000.ukbb_bca_cases.chr16.52519188-52679188_2,j.txt Pattern_combined_old.Iteration000.ukbb_bca_cases.chr16.52519188-52679188_2,j.txt

# mv Pattern_combined.Iteration000.ukbb_bca_cases.chr16.52519188-52679188_2,j.tmp Pattern_combined.Iteration000.ukbb_bca_cases.chr16.52519188-52679188_2,j.txt

# remove old files from previous run (can skip if running iterations for the first time)
# for file in *Iteration00[1-9]*; do rm -r $file; done # keeps Iteration000 and removes Iterations 001-009; if Iterations 010 and above exist, need to rewrite code to remove them

# for file in fisher_exact*; do rm -r $file; done

# for file in Core*; do rm $file; done

# for file in Closed*; do rm $file; done

# Begin iterations:
# bsub -P SJLIFE -J pattern_overlap_iterate3.ukbb_bca_cases.chr16.52519188-52679188 -oo /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/pattern_overlap_iterate3.ukbb_bca_cases.chr16.52519188-52679188.out -eo /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/pattern_overlap_iterate3.ukbb_bca_cases.chr16.52519188-52679188.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/pattern_overlap_iterate3.sh haplotype_estimates.ukbb_bca_cases.chr16.52519188-52679188.txt,haplotype_estimates.ukbb_bca_controls.chr16.52519188-52679188.txt 16 52519188,52679188 100 50 2 2,j \"\" \"ukbb_bca_cases\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1 /home/wletsou/scripts"

# details on iterations, make plots of total and closed patterns
# sh /home/wletsou/scripts/ukbb_pattern_count_summary.sh ukbb_bca_cases 16 52519188,52679188 2,j

# Haplotypes with p < 1e-XX
# test -f fisher_exact.ukbb_bca_cases.Results.txt && rm fisher_exact.ukbb_bca_cases.Results.txt; touch fisher_exact.ukbb_bca_cases.Results.txt; for file in fisher_exact.ukbb_bca_cases.Iteration*.txt; do awk '($6+0<1e-36){print $0}' $file >> fisher_exact.ukbb_bca_cases.Results.txt; done;

# translated significant haplotypes from column numbers to rsids
# bsub -P SJLIFE -J ukbb_haplotype_translate2_sub -oo ukbb_haplotype_translate2_sub.out -eo ukbb_haplotype_translate2_sub.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52519188-52679188.txt fisher_exact.ukbb_bca_cases.Results.txt 50"

# counts of translated haplotypes on each chromosome
# bsub -P SJLIFE -J ukbb_haplotype_model9_sub -oo ukbb_haplotype_model9_sub.out -eo ukbb_haplotype_model9_sub.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52519188-52679188.txt fisher_exact.ukbb_bca_cases.Results.translated.txt \"\" \"\" \"\" 16 52519188,52679188 50 \"\" \"\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1 /home/wletsou/scripts"

# add patterns to model
# bsub -P SJLIFE -J ukbb_haplotype_model9_iterate -oo ukbb_haplotype_model9_iterate.out -eo ukbb_haplotype_model9_iterate.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.sh allele_counts.txt /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt \"\" \"\" 1e-5 Significant_patterns.txt /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1 /home/wletsou/scripts"

# awk 'BEGIN{n=split("rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs45625634_C=0,rs45496893_T=0,rs34750829_A=0,rs9931232_A=1,rs45477396_C=0,rs75571494_A=0,rs4594251_T=0,rs4784227_T=1,rs11864809_C=0,rs118162666_A=0,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)_[A-Z]=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } ($3 in haplotype_snps){found+=1; printf "%s[%s]%s",$3,(haplotype_snps[$3]==1?$5:$4),(found<n?",":"\n")}' ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr16.52519188-52679188.haplotypes.hg38.vcf # h1 alleles
# rs78338509[A],rs12447430[T],rs45577538[G],rs8046985[A],rs45542333[G],rs78268044[C],rs45454402[A],rs45512493[A],rs45625634[T],rs45496893[C],rs34750829[G],rs9931232[A],rs45477396[T],rs75571494[G],rs4594251[C],rs4784227[T],rs11864809[T],rs118162666[G],rs45584434[G],rs45560737[C],rs45575339[C],rs45587544[G],rs78841172[T],rs76000465[T],rs111925335[A],rs75127968[G]

# reduction
# sh /home/wletsou/scripts/haplotype_rpart.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52519188-52679188.txt rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs45625634_C=0,rs45496893_T=0,rs34750829_A=0,rs9931232_A=1,rs45477396_C=0,rs75571494_A=0,rs4594251_T=0,rs4784227_T=1,rs11864809_C=0,rs118162666_A=0,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0 "" 1 # h1
#  reduced = rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0

# awk 'BEGIN{n=split("rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0",array,","); for (i=1;i<=n;i++) {a=gensub("(.*)_[A-Z]=(.*)","\\1","g",array[i]); b=gensub("(.*)=(.*)","\\2","g",array[i]); haplotype_snps[a]=b} } ($3 in haplotype_snps){found+=1; printf "%s[%s]%s",$3,(haplotype_snps[$3]==1?$5:$4),(found<n?",":"\n")}' ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr16.52519188-52679188.haplotypes.hg38.vcf # h1.reduced alleles
# rs78338509[A],rs12447430[T],rs45577538[G],rs8046985[A],rs45542333[G],rs78268044[C],rs45454402[A],rs45512493[A],rs34750829[G],rs9931232[A],rs4594251[C],rs4784227[T],rs45584434[G],rs45560737[C],rs45575339[C],rs45587544[G],rs78841172[T],rs76000465[T],rs111925335[A],rs75127968[G]

# check that h1 and h1.reduced have the same effect sizes
# sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52519188-52679188.txt "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs45625634_C=0,rs45496893_T=0,rs34750829_A=0,rs9931232_A=1,rs45477396_C=0,rs75571494_A=0,rs4594251_T=0,rs4784227_T=1,rs11864809_C=0,rs118162666_A=0,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "" "" "" "" h1 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1 /home/wletsou/scripts
# LRT p value = 1.234118e-41 on 1 degree(s) of freedom with HR = 1.265392e+00 (1.223759e+00 to 1.308441e+00), frequency = 2.139874e-01 (77478) = 2.549662e-01/2.118409e-01 (4595/72883)

# sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52519188-52679188.txt "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "" "" "" "" h1.reduced /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1 /home/wletsou/scripts
# LRT p value = 1.234118e-41 on 1 degree(s) of freedom with HR = 1.265392e+00 (1.223759e+00 to 1.308441e+00), frequency = 2.139874e-01 (77478) = 2.549662e-01/2.118409e-01 (4595/72883)

# effect of GWAS hit
# sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52519188-52679188.txt "rs4784227_T=1" "" "" "" "" rs4784227 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1 /home/wletsou/scripts
# LRT p value = 1.170830e-40 on 1 degree(s) of freedom with HR = 1.252320e+00 (1.212312e+00 to 1.293649e+00), frequency = 2.389303e-01 (86509) = 2.808234e-01/2.367358e-01 (5061/81448)

# bsub -P SJLIFE -J ukbb_haplotype_haplotype_LD2 -oo /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/ukbb_haplotype_haplotype_LD2.out -eo /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/ukbb_haplotype_haplotype_LD2.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ukbb_haplotype_haplotype_LD2.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52519188-52679188.txt \"rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0\" \"rs4784227_T=1\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1" # h1 (reduced) vs GWAS hit
# 0.213987	0.23893	0.213987
# rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0	rs4784227_T=1	0.867186	1

# plot SNPs

# echo "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" | sed 's/_[A-Z]=[0-9]//g' # snps in h1, no alleles
# convert .vcf to .bed at SNPs in h1 (reduced)
# sh /home/wletsou/scripts/vcf2bed.sh ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.hg19_chr16.52519188-52679188.haplotypes.hg38.vcf "rs78338509,rs12447430,rs45577538,rs8046985,rs45542333,rs78268044,rs45454402,rs45512493,rs34750829,rs9931232,rs4594251,rs4784227,rs45584434,rs45560737,rs45575339,rs45587544,rs78841172,rs76000465,rs111925335,rs75127968" "" h1.bed /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1 /home/wletsou/scripts # GRCh38 bed file

# mkdir phase2
# cd phase2

# get h1 and list of SNPs in remainder of TAD
# HAPLOTYPE_SNPS=$(echo "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" | sed 's/_[A-Z]=[0-9]//g') # snps in h1, no alleles

# awk 'BEGIN{OFS="\t"} ( ($2>52679188 && $2<=53075000) || "'$HAPLOTYPE_SNPS'" ~ $3){print $3}' /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/ukbb_genotyped_snps/ukbb_bca_cases+ukbb_bca_controls.chr16.52000000-53200000.haplotypes.vcf > TAD_snps.chr16.52679189-53075000.txt # http://3dgenome.fsm.northwestern.edu/view.php?method=Hi-C&species=human&assembly=hg19&source=inside&tissue=HMEC&type=Lieberman-raw&resolution=5&c_url=&transfer=&chr=chr16&start=52000000&end=53200000&sessionID=&browser=ucsc

# check that SNPs have good "imputation" quality
# awk 'NR==FNR{snp[$1]; next} ($8 in snp && $15==1){print $8}' TAD_snps.chr16.52679189-53075000.txt /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr16.51560000.53560000.qced.anno.info > TAD_snps.chr16.52679189-53075000.tmp
# test -f TAD_snps.chr16.52679189-53075000.tmp && mv TAD_snps.chr16.52679189-53075000.tmp TAD_snps.chr16.52679189-53075000.txt

# awk 'BEGIN{OFS="\t"} ($3==1){print $1,$1}' /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_cases.indiv
# awk 'BEGIN{OFS="\t"} ($3==0){print $1,$1}' /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt > ukbb_bca_controls.indiv

# vcf file for UKBB
# bsub -P SJLIFE -J ukbb_topmed_extract.chr16.52679189-53075000 -oo ukbb_topmed_extract.chr16.52679189-53075000.out -eo ukbb_topmed_extract.chr16.52679189-53075000.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_topmed_extract.sh TAD_snps.chr16.52679189-53075000.txt /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr16.51560000.53560000.qced.anno.info 16 52679189,53075000 /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/ukb.bca.topmedr2.chr16.51560000.53560000.qced.hg38.vcf.gz"

# vcf file for dbGaP
# bsub -P SJLIFE -J dbgap_haplotype_merge.chr16.52679189-53075000 -oo dbgap_haplotype_merge.chr16.52679189-53075000.out -eo dbgap_haplotype_merge.chr16.52679189-53075000.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/dbgap_haplotype_merge.sh dbgap28544_cases,dbgap28544_controls "/research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt" "/research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/chr16.qced.info" TAD_snps.chr16.52679189-53075000.txt 16 52679189,53075000 /research/projects/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/topmedr2/cleaned/vcfs/dbgap28544.topmedr2.cleaned.hg38.chr16.vcf.gz"

# check that there are no overlapping samples
# module load bcftools/1.10.2

# awk 'NR==FNR{seen[$1]; next} ($1 in seen){print $0}' <(bcftools query -l dbgap28544_cases+dbgap28544_controls.hg19_chr16.52679189-53075000.hg38.vcf.gz) <(bcftools query -l ukbb.topmed.hg19_chr16.52679189-53075000.hg38.vcf.gz) # dbGaP samples in UKBB, if any

# awk 'NR==FNR{seen[$1]; next} ($1 in seen){print $0}' <(bcftools query -l ukbb.topmed.hg19_chr16.52679189-53075000.hg38.vcf.gz) <(bcftools query -l dbgap28544_cases+dbgap28544_controls.hg19_chr16.52679189-53075000.hg38.vcf.gz) # UKBB samples in dbGaP, if any

# check for missing SNPs
# awk 'NR==FNR{seen[$3]; next} ($3 in seen==0){print $3}' <(bcftools view --no-header dbgap28544_cases+dbgap28544_controls.hg19_chr16.52679189-53075000.hg38.vcf.gz) <(bcftools view --no-header ukbb.topmed.hg19_chr16.52679189-53075000.hg38.vcf.gz) # SNPs in UKBB not in DRIVE, if any

# awk 'NR==FNR{seen[$3]; next} ($3 in seen==0){print $3}' <(bcftools view --no-header ukbb.topmed.hg19_chr16.52679189-53075000.hg38.vcf.gz) <(bcftools view --no-header dbgap28544_cases+dbgap28544_controls.hg19_chr16.52679189-53075000.hg38.vcf.gz) # SNPs in DRIVE not in UKBB, if any

# first round of overlaps (of h1-carrying chromosomes among cases)
# bsub -P SJLIFE -J ukbb_bca_overlap.chr16.52679189-53075000 -oo ukbb_bca_overlap.chr16.52679189-53075000.out -eo ukbb_bca_overlap.chr16.52679189-53075000.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_submission13.sh ukbb_bca_cases,ukbb_bca_controls,dbgap28544_cases,dbgap28544_controls 16 52679189,53075000 "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" TAD_snps.chr16.52679189-53075000.txt ukbb.topmed.hg19_chr16.52679189-53075000.hg38.vcf.gz,dbgap28544_cases+dbgap28544_controls.hg19_chr16.52679189-53075000.hg38.vcf.gz /research/rgs01/project_space/yasuigrp/EpiGenetics/common/drive4william/chr16.51560000.53560000.qced.anno.info 0.0001 \"\" \"\" \"\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts"

# set a threshold for remaining overlaps
# mv Pattern_combined_old.Iteration000.ukbb_bca_cases.chr16.52679189-53075000_2,j.txt Pattern_combined.Iteration000.ukbb_bca_cases.chr16.52679189-53075000_2,j.txt

# awk 'BEGIN{OFS="\t"} NR==FNR{count[$2]=$1; next} ($1 in count){print count[$1],$1}' Pattern_combined.Iteration000.ukbb_bca_cases.chr16.52679189-53075000_2,j.txt <(awk '($5+0<6.3e-6){print $0}' Pattern_differences.Iteration000.ukbb_bca_cases.chr16.52679189-53075000_2,j.p_vals.txt | awk 'BEGIN{array[0]} {delete array; n=split($1,array,","); if (n<snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]+0 || snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]+0==0) {snps[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]=n; seen[sprintf("%s\t%s\t%s\t%s",$2,$3,$4,$5)]=$1} } END{for (i in seen) {print seen[i]} }') > Pattern_combined.Iteration000.ukbb_bca_cases.chr16.52679189-53075000_2,j.tmp

# mv Pattern_combined.Iteration000.ukbb_bca_cases.chr16.52679189-53075000_2,j.txt Pattern_combined_old.Iteration000.ukbb_bca_cases.chr16.52679189-53075000_2,j.txt

# mv Pattern_combined.Iteration000.ukbb_bca_cases.chr16.52679189-53075000_2,j.tmp Pattern_combined.Iteration000.ukbb_bca_cases.chr16.52679189-53075000_2,j.txt

# remove old files from previous run (can skip if running iterations for the first time)
# for file in *Iteration00[1-9]*; do rm -r $file; done # keeps Iteration000 and removes Iterations 001-009; if Iterations 010 and above exist, need to rewrite code to remove them

# for file in fisher_exact*; do rm -r $file; done

# for file in Core*; do rm $file; done

# for file in Closed*; do rm $file; done

# Begin iterations:
# bsub -P SJLIFE -J pattern_overlap_iterate3.ukbb_bca_cases.chr16.52679189-53075000 -oo /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2/pattern_overlap_iterate3.ukbb_bca_cases.chr16.52679189-53075000.out -eo /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2/pattern_overlap_iterate3.ukbb_bca_cases.chr16.52679189-53075000.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/pattern_overlap_iterate3.sh haplotype_estimates.ukbb_bca_cases.chr16.52679189-53075000.subset.txt,haplotype_estimates.ukbb_bca_controls.chr16.52679189-53075000.subset.txt 16 52679189,53075000 100 50 2 2,j \"\" \"ukbb_bca_cases\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts"

# sh /home/wletsou/scripts/ukbb_pattern_count_summary.sh ukbb_bca_cases 16 52679189,53075000 2,j

# combined haplotype file for cases and controls
# awk 'BEGIN{OFS="\t"} {print $0}' ukbb_bca_cases.indiv ukbb_bca_controls.indiv > ukbb_bca_cases+ukbb_bca_controls.indiv
# awk 'NR==FNR{id[$1]; next} (NR!=FNR && FNR==1){print $0} (NR!=FNR && FNR>1 && $1 in id){print $0}' ukbb_bca_cases+ukbb_bca_controls.indiv haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.chr16.52679189-53075000.txt > haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52679189-53075000.txt

# forward selection
# translated in parallel, p <
# test -f fisher_exact.ukbb_bca_cases.Results.txt && rm fisher_exact.ukbb_bca_cases.Results.txt; touch fisher_exact.ukbb_bca_cases.Results.txt; for file in fisher_exact.ukbb_bca_cases.Iteration*.txt; do awk '($6+0<1e-2){print $0}' $file >> fisher_exact.ukbb_bca_cases.Results.txt; done;

# remove longer patterns with the same frequency, OR, and p-value as a shorter pattern
# awk 'BEGIN{array[0]} {delete array; n=split($1,array,","); if (n<snps[sprintf("%s\t%s\t%s\t%s",$3,$4,$5,$6)]+0 || snps[sprintf("%s\t%s\t%s\t%s",$3,$4,$5,$6)]+0==0) {snps[sprintf("%s\t%s\t%s\t%s",$3,$4,$5,$6)]=n; seen[sprintf("%s\t%s\t%s\t%s",$3,$4,$5,$6)]=$0} }  END{for (i in seen) {print seen[i]} }' fisher_exact.ukbb_bca_cases.Results.txt > fisher_exact.ukbb_bca_cases.Results.tmp; test -f fisher_exact.ukbb_bca_cases.Results.tmp && mv fisher_exact.ukbb_bca_cases.Results.tmp fisher_exact.ukbb_bca_cases.Results.txt

# bsub -P SJLIFE -J ukbb_haplotype_translate2_sub -oo ukbb_haplotype_translate2_sub.out -eo ukbb_haplotype_translate2_sub.err -R "rusage[mem=256]" -q standard "sh /home/wletsou/scripts/ukbb_haplotype_translate2_sub.sh haplotype_estimates.ukbb_bca_cases.chr16.52679189-53075000.subset.txt fisher_exact.ukbb_bca_cases.Results.txt 50"

# bsub -P SJLIFE -J ukbb_haplotype_model9_sub -oo ukbb_haplotype_model9_sub.out -eo ukbb_haplotype_model9_sub.err -R "rusage[mem=256]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_sub.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52679189-53075000.txt fisher_exact.ukbb_bca_cases.Results.translated.txt 1 \"rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0\" 1 16 52679189,53075000 50 \"\" \"\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts"

# bsub -P SJLIFE -J ukbb_haplotype_model9_iterate -oo ukbb_haplotype_model9_iterate.out -eo ukbb_haplotype_model9_iterate.err -R "rusage[mem=1000]" "sh /home/wletsou/scripts/ukbb_haplotype_model9_iterate.sh allele_counts.txt /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt 1 1 1e-5 Significant_patterns.txt /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts"

# HAPLOTYPE_LIST=$(awk '(NR>1 && $6+0<1e-5 && $3 > 1){printf "%s%s",(found>0?":":""),$2; found+=1} END{printf "\n"}' Significant_patterns.txt)

# sh /home/wletsou/scripts/ukbb_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/ukb.bca/hap.gds/ukb.bca.pheno.txt haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52679189-53075000.txt "rs12324961_C=1,rs3112586_T=1,rs4238751_G=1,rs80176221_G=0,rs7186592_G=1,rs118136222_T=0,rs118155214_C=0,rs77700051_C=0,rs117106695_C=0,rs118186521_C=0,rs16951525_C=0,rs117923020_A=0,rs79501461_T=0,rs77362202_T=0,rs117155117_G=0,rs62043164_G=0,rs17361619_C=0,rs77960456_G=0,rs117544975_C=0,rs1362413_A=0,rs59793894_G=0,rs74020833_C=0,rs17296337_G=0,rs79207645_A=0,rs117947188_C=0,rs61131526_A=0,rs117553291_G=0,rs17296940_C=0,rs76028093_C=0,rs17297072_T=0,rs79291013_G=0,rs117746278_A=0,rs17298178_C=0,rs79993873_C=0,rs78320298_G=0,rs62041381_C=0,rs117551222_C=0,rs7500472_A=0,rs72810045_A=0,rs12446016_G=0,rs150704739_T=0,rs142367768_T=0,rs4783791_G=0,rs117911537_T=0,rs1344484_C=0,rs16951771_A=0,rs73597288_A=0,rs6499147_C=1,rs12922473_G=0,rs78322982_A=0,rs6499148_A=1,rs62043282_A=0,rs113717347_A=0,rs62043323_A=0,rs12924810_A=0,rs80356178_A=0,rs34261756_C=0,rs117400881_C=0,rs1477029_A=1,rs7199322_G=0,rs4784276_C=0,rs77419902_T=0,rs78724580_A=0,rs72812119_G=0,rs72812122_A=0,rs76365437_A=0,rs13335861_G=0,rs4352043_C=1,rs36099013_A=0,rs72812136_T=0,rs79815180_T=0,rs74018575_C=0,rs34382563_T=0,rs12598778_C=0,rs79721160_A=0,rs117297253_A=0,rs3931698_T=1:rs12324961_C=0,rs3112586_T=1,rs12927162_G=0,rs118136222_T=0,rs118155214_C=0,rs77700051_C=0,rs117106695_C=0,rs118186521_C=0,rs16951525_C=0,rs117923020_A=0,rs79501461_T=0,rs117155117_G=0,rs62043164_G=0,rs17361619_C=0,rs77960456_G=0,rs113555472_A=0,rs59793894_G=0,rs4784253_G=0,rs74020833_C=0,rs17296337_G=0,rs4377151_G=0,rs79207645_A=0,rs117947188_C=0,rs61131526_A=0,rs117553291_G=0,rs17296940_C=0,rs76028093_C=0,rs17297072_T=0,rs79291013_G=0,rs8044632_C=0,rs117746278_A=0,rs9939896_T=0,rs79993873_C=0,rs78320298_G=0,rs62041381_C=0,rs117222080_T=0,rs117551222_C=0,rs7500472_A=0,rs72810045_A=0,rs12446016_G=0,rs73591253_T=0,rs113964117_A=0,rs142367768_T=0,rs4783791_G=0,rs117911537_T=0,rs1344484_C=0,rs73597288_A=0,rs12919591_G=0,rs9630629_C=0,rs78322982_A=0,rs6499148_A=1,rs62043282_A=0,rs113717347_A=0,rs62043323_A=0,rs12924810_A=0,rs80356178_A=0,rs34261756_C=0,rs117400881_C=0,rs1477029_A=1,rs7199322_G=0,rs4784276_C=0,rs77419902_T=0,rs78724580_A=0,rs72812119_G=0,rs72812122_A=0,rs62047145_T=0,rs76365437_A=0,rs13335861_G=0,rs36099013_A=0,rs72812136_T=0,rs74018575_C=0,rs34382563_T=0,rs12598778_C=0,rs56085964_C=1,rs7499206_C=1,rs12922187_G=0,rs8051957_G=1,rs8050297_A=1,rs79721160_A=0,rs79088862_A=0,rs4784287_T=1,rs117297253_A=0,rs3931698_T=1:rs12324961_C=0,rs3112586_T=1,rs12927162_G=0,rs80176221_G=0,rs118136222_T=0,rs118155214_C=0,rs77700051_C=0,rs11075587_C=0,rs117106695_C=0,rs118186521_C=0,rs16951525_C=0,rs117923020_A=0,rs79501461_T=0,rs77362202_T=0,rs117155117_G=0,rs62043164_G=0,rs17361619_C=0,rs77960456_G=0,rs117544975_C=0,rs1362413_A=0,rs113555472_A=0,rs59793894_G=0,rs4784253_G=0,rs74020833_C=0,rs117947188_C=0,rs61131526_A=0,rs117553291_G=0,rs17296940_C=0,rs76028093_C=0,rs17297072_T=0,rs79291013_G=0,rs117746278_A=0,rs17298178_C=0,rs79993873_C=0,rs78320298_G=0,rs62041381_C=0,rs117551222_C=0,rs7500472_A=0,rs72810045_A=0,rs12446016_G=0,rs73591253_T=0,rs150704739_T=0,rs113964117_A=0,rs1345327_C=0,rs142367768_T=0,rs4783791_G=0,rs117911537_T=0,rs1344484_C=0,rs16951771_A=0,rs73597288_A=0,rs12919591_G=0,rs9630629_C=1,rs6499147_C=1,rs12922473_G=0,rs78322982_A=0,rs6499148_A=1,rs62043282_A=0,rs113717347_A=0,rs62043323_A=0,rs12924810_A=0,rs80356178_A=0,rs34261756_C=0,rs117400881_C=0,rs7199322_G=0,rs78831406_T=0,rs4784276_C=0,rs78724580_A=0,rs72812119_G=0,rs72812122_A=0,rs76365437_A=0,rs13335861_G=0,rs4352043_C=1,rs36099013_A=0,rs72812136_T=0,rs79815180_T=0,rs72812140_T=0,rs34382563_T=0,rs12598778_C=0,rs56085964_C=0,rs7499206_C=0,rs8051957_G=1,rs8050297_A=1,rs79088862_A=0,rs117297253_A=0,rs3931698_T=1:rs3112586_T=1,rs4238751_G=1,rs80176221_G=0,rs7186592_G=1,rs118136222_T=0,rs118155214_C=0,rs77700051_C=0,rs11075587_C=1,rs117106695_C=0,rs118186521_C=0,rs16951525_C=0,rs117923020_A=0,rs79501461_T=0,rs77362202_T=0,rs117155117_G=0,rs62043164_G=0,rs17361619_C=0,rs77960456_G=0,rs117544975_C=0,rs1362413_A=1,rs4784249_C=1,rs113555472_A=0,rs59793894_G=0,rs74020833_C=0,rs17296337_G=0,rs79207645_A=0,rs117947188_C=0,rs61131526_A=0,rs117553291_G=0,rs17296940_C=0,rs76028093_C=0,rs17297072_T=0,rs117746278_A=0,rs17298178_C=0,rs79993873_C=0,rs78320298_G=0,rs62041381_C=0,rs117222080_T=0,rs117551222_C=0,rs7500472_A=0,rs72810045_A=0,rs12446016_G=0,rs73591253_T=0,rs150704739_T=0,rs113964117_A=0,rs142367768_T=0,rs117911537_T=0,rs16951771_A=0,rs73597288_A=0,rs12919591_G=0,rs9630629_C=0,rs6499147_C=1,rs12922473_G=1,rs78322982_A=0,rs62043282_A=0,rs113717347_A=0,rs62043323_A=0,rs12924810_A=0,rs34261756_C=0,rs117400881_C=0,rs78831406_T=0,rs4784276_C=0,rs77419902_T=0,rs78724580_A=0,rs72812119_G=0,rs62047145_T=0,rs76365437_A=0,rs13335861_G=0,rs36099013_A=0,rs72812136_T=0,rs79815180_T=0,rs74018575_C=0,rs72812140_T=0,rs34382563_T=0,rs56085964_C=0,rs12922187_G=0,rs79721160_A=0,rs79088862_A=0,rs4784287_T=1,rs117297253_A=0,rs56278234_T=0" "1" "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "1" "" h2-h5 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts
# Rscript /home/wletsou/scripts/ukbb_haplotype_model5.R file=/scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2/h2-h5.new_allele_counts.txt groups=1,2,3,4,5 verbose=1 n=4 n_covs=10
# [1] "X00000" "X10000" "X10001" "X10010" "X10100" "X11000"
# [1] "DT[,h0:=X00000]"
# [1] "DT[,h1.5:=X10000+X10001]"
# [1] "DT[,h1:=X10000]"
# [1] "DT[,h2:=X11000]"
# [1] "DT[,h3:=X10100]"
# [1] "DT[,h4:=X10010]"
# [1] "DT[,h5:=X10001]"
#
#    affected  h1.5    h1 h2 h3 h4  h5
# 1:        0 72741 72636 17 74 51 105
# 2:        1  4547  4525 13 20 15  22
#
#    affected      h1.5        h1           h2           h3           h4           h5
# 1:        0 0.2114281 0.2111229 0.0000494120 0.0002150875 0.0001482360 0.0003051917
# 2:        1 0.2523027 0.2510820 0.0007213406 0.0011097547 0.0008323161 0.0012207302
#
#     h1.5    h1 h2 h3 h4  h5
# 1: 77288 77161 30 94 66 127
#
#         h1.5        h1           h2           h3           h4           h5
# 1: 0.2134627 0.2131119 8.285736e-05 0.0002596197 0.0001822862 0.0003507628
#
# Surv(age,affected) ~ h1 + h2 + h3 + h4 + h5 + pc01 + pc02 + pc03 + pc04 + pc05 + pc06 + pc07 + pc08 + pc09 + pc10
# Call:
# coxph(formula = eval(parse(text = Y)), data = X)
#
#   n= 181034, number of events= 9011
#
#            coef  exp(coef)   se(coef)      z Pr(>|z|)
# h1    0.2243498  1.2515088  0.0171644 13.071  < 2e-16 ***
# h2    2.3267167 10.2442508  0.2777327  8.378  < 2e-16 ***
# h3    1.6527768  5.2214589  0.2239707  7.379 1.59e-13 ***
# h4    1.7469298  5.7369622  0.2585168  6.758 1.40e-11 ***
# h5    1.3498115  3.8566985  0.2135828  6.320 2.62e-10 ***
# pc01 -0.0007036  0.9992967  0.0069111 -0.102   0.9189
# pc02 -0.0143272  0.9857750  0.0071544 -2.003   0.0452 *
# pc03  0.0070480  1.0070728  0.0069042  1.021   0.3073
# pc04  0.0036264  1.0036330  0.0051771  0.700   0.4836
# pc05 -0.0007691  0.9992312  0.0022685 -0.339   0.7346
# pc06  0.0073484  1.0073754  0.0065749  1.118   0.2637
# pc07  0.0060906  1.0061092  0.0059062  1.031   0.3024
# pc08 -0.0018801  0.9981217  0.0058789 -0.320   0.7491
# pc09 -0.0022407  0.9977618  0.0022845 -0.981   0.3267
# pc10  0.0061152  1.0061339  0.0050968  1.200   0.2302
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#
#      exp(coef) exp(-coef) lower .95 upper .95
# h1      1.2515    0.79904    1.2101    1.2943
# h2     10.2443    0.09762    5.9439   17.6558
# h3      5.2215    0.19152    3.3663    8.0991
# h4      5.7370    0.17431    3.4565    9.5221
# h5      3.8567    0.25929    2.5375    5.8616
# pc01    0.9993    1.00070    0.9859    1.0129
# pc02    0.9858    1.01443    0.9720    0.9997
# pc03    1.0071    0.99298    0.9935    1.0208
# pc04    1.0036    0.99638    0.9935    1.0139
# pc05    0.9992    1.00077    0.9948    1.0037
# pc06    1.0074    0.99268    0.9945    1.0204
# pc07    1.0061    0.99393    0.9945    1.0178
# pc08    0.9981    1.00188    0.9867    1.0097
# pc09    0.9978    1.00224    0.9933    1.0022
# pc10    1.0061    0.99390    0.9961    1.0162
#
# Concordance= 0.541  (se = 0.003 )
# Likelihood ratio test= 295  on 15 df,   p=<2e-16
# Wald test            = 378.1  on 15 df,   p=<2e-16
# Score (logrank) test = 442.2  on 15 df,   p=<2e-16

# backward selection
# module load R/3.6.1
# Rscript /home/wletsou/scripts/ukbb_backward_selection2.R file=h2-h5.new_allele_counts.txt n=4 groups=1,2,3,4,5 threshold=0.05 n_covs=10
# remaining_vars         pval
# 1:             h2 1.716370e-08
# 2:             h3 2.431043e-07
# 3:             h4 2.479390e-06
# 4:             h5 8.870288e-06

# reduction
# sh /home/wletsou/scripts/haplotype_rpart.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52679189-53075000.txt "rs12324961_C=1,rs3112586_T=1,rs4238751_G=1,rs80176221_G=0,rs7186592_G=1,rs118136222_T=0,rs118155214_C=0,rs77700051_C=0,rs117106695_C=0,rs118186521_C=0,rs16951525_C=0,rs117923020_A=0,rs79501461_T=0,rs77362202_T=0,rs117155117_G=0,rs62043164_G=0,rs17361619_C=0,rs77960456_G=0,rs117544975_C=0,rs1362413_A=0,rs59793894_G=0,rs74020833_C=0,rs17296337_G=0,rs79207645_A=0,rs117947188_C=0,rs61131526_A=0,rs117553291_G=0,rs17296940_C=0,rs76028093_C=0,rs17297072_T=0,rs79291013_G=0,rs117746278_A=0,rs17298178_C=0,rs79993873_C=0,rs78320298_G=0,rs62041381_C=0,rs117551222_C=0,rs7500472_A=0,rs72810045_A=0,rs12446016_G=0,rs150704739_T=0,rs142367768_T=0,rs4783791_G=0,rs117911537_T=0,rs1344484_C=0,rs16951771_A=0,rs73597288_A=0,rs6499147_C=1,rs12922473_G=0,rs78322982_A=0,rs6499148_A=1,rs62043282_A=0,rs113717347_A=0,rs62043323_A=0,rs12924810_A=0,rs80356178_A=0,rs34261756_C=0,rs117400881_C=0,rs1477029_A=1,rs7199322_G=0,rs4784276_C=0,rs77419902_T=0,rs78724580_A=0,rs72812119_G=0,rs72812122_A=0,rs76365437_A=0,rs13335861_G=0,rs4352043_C=1,rs36099013_A=0,rs72812136_T=0,rs79815180_T=0,rs74018575_C=0,rs34382563_T=0,rs12598778_C=0,rs79721160_A=0,rs117297253_A=0,rs3931698_T=1" "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" 1 # h2
#  reduced = rs12324961_C=1,rs3112586_T=1,rs7186592_G=1,rs16951525_C=0,rs117923020_A=0,rs79501461_T=0,rs1362413_A=0,rs59793894_G=0,rs17296337_G=0,rs79207645_A=0,rs117553291_G=0,rs79291013_G=0,rs117746278_A=0,rs150704739_T=0,rs1344484_C=0,rs12922473_G=0,rs6499148_A=1,rs117400881_C=0,rs1477029_A=1,rs77419902_T=0,rs72812119_G=0,rs79815180_T=0

# sh /home/wletsou/scripts/haplotype_rpart.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52679189-53075000.txt "rs12324961_C=0,rs3112586_T=1,rs12927162_G=0,rs118136222_T=0,rs118155214_C=0,rs77700051_C=0,rs117106695_C=0,rs118186521_C=0,rs16951525_C=0,rs117923020_A=0,rs79501461_T=0,rs117155117_G=0,rs62043164_G=0,rs17361619_C=0,rs77960456_G=0,rs113555472_A=0,rs59793894_G=0,rs4784253_G=0,rs74020833_C=0,rs17296337_G=0,rs4377151_G=0,rs79207645_A=0,rs117947188_C=0,rs61131526_A=0,rs117553291_G=0,rs17296940_C=0,rs76028093_C=0,rs17297072_T=0,rs79291013_G=0,rs8044632_C=0,rs117746278_A=0,rs9939896_T=0,rs79993873_C=0,rs78320298_G=0,rs62041381_C=0,rs117222080_T=0,rs117551222_C=0,rs7500472_A=0,rs72810045_A=0,rs12446016_G=0,rs73591253_T=0,rs113964117_A=0,rs142367768_T=0,rs4783791_G=0,rs117911537_T=0,rs1344484_C=0,rs73597288_A=0,rs12919591_G=0,rs9630629_C=0,rs78322982_A=0,rs6499148_A=1,rs62043282_A=0,rs113717347_A=0,rs62043323_A=0,rs12924810_A=0,rs80356178_A=0,rs34261756_C=0,rs117400881_C=0,rs1477029_A=1,rs7199322_G=0,rs4784276_C=0,rs77419902_T=0,rs78724580_A=0,rs72812119_G=0,rs72812122_A=0,rs62047145_T=0,rs76365437_A=0,rs13335861_G=0,rs36099013_A=0,rs72812136_T=0,rs74018575_C=0,rs34382563_T=0,rs12598778_C=0,rs56085964_C=1,rs7499206_C=1,rs12922187_G=0,rs8051957_G=1,rs8050297_A=1,rs79721160_A=0,rs79088862_A=0,rs4784287_T=1,rs117297253_A=0,rs3931698_T=1" "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" 1 # h3
#  reduced = rs12324961_C=0,rs3112586_T=1,rs12927162_G=0,rs79501461_T=0,rs4784253_G=0,rs17296337_G=0,rs4377151_G=0,rs79207645_A=0,rs117553291_G=0,rs17297072_T=0,rs79291013_G=0,rs78320298_G=0,rs72810045_A=0,rs73591253_T=0,rs1344484_C=0,rs12919591_G=0,rs9630629_C=0,rs78322982_A=0,rs80356178_A=0,rs117400881_C=0,rs1477029_A=1,rs77419902_T=0,rs72812119_G=0,rs56085964_C=1,rs12922187_G=0,rs8051957_G=1,rs8050297_A=1,rs4784287_T=1

# sh /home/wletsou/scripts/haplotype_rpart.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52679189-53075000.txt "rs12324961_C=0,rs3112586_T=1,rs12927162_G=0,rs80176221_G=0,rs118136222_T=0,rs118155214_C=0,rs77700051_C=0,rs11075587_C=0,rs117106695_C=0,rs118186521_C=0,rs16951525_C=0,rs117923020_A=0,rs79501461_T=0,rs77362202_T=0,rs117155117_G=0,rs62043164_G=0,rs17361619_C=0,rs77960456_G=0,rs117544975_C=0,rs1362413_A=0,rs113555472_A=0,rs59793894_G=0,rs4784253_G=0,rs74020833_C=0,rs117947188_C=0,rs61131526_A=0,rs117553291_G=0,rs17296940_C=0,rs76028093_C=0,rs17297072_T=0,rs79291013_G=0,rs117746278_A=0,rs17298178_C=0,rs79993873_C=0,rs78320298_G=0,rs62041381_C=0,rs117551222_C=0,rs7500472_A=0,rs72810045_A=0,rs12446016_G=0,rs73591253_T=0,rs150704739_T=0,rs113964117_A=0,rs1345327_C=0,rs142367768_T=0,rs4783791_G=0,rs117911537_T=0,rs1344484_C=0,rs16951771_A=0,rs73597288_A=0,rs12919591_G=0,rs9630629_C=1,rs6499147_C=1,rs12922473_G=0,rs78322982_A=0,rs6499148_A=1,rs62043282_A=0,rs113717347_A=0,rs62043323_A=0,rs12924810_A=0,rs80356178_A=0,rs34261756_C=0,rs117400881_C=0,rs7199322_G=0,rs78831406_T=0,rs4784276_C=0,rs78724580_A=0,rs72812119_G=0,rs72812122_A=0,rs76365437_A=0,rs13335861_G=0,rs4352043_C=1,rs36099013_A=0,rs72812136_T=0,rs79815180_T=0,rs72812140_T=0,rs34382563_T=0,rs12598778_C=0,rs56085964_C=0,rs7499206_C=0,rs8051957_G=1,rs8050297_A=1,rs79088862_A=0,rs117297253_A=0,rs3931698_T=1" "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" 1 # h4
#  reduced = rs12324961_C=0,rs3112586_T=1,rs80176221_G=0,rs11075587_C=0,rs16951525_C=0,rs117923020_A=0,rs77362202_T=0,rs117544975_C=0,rs113555472_A=0,rs4784253_G=0,rs17297072_T=0,rs17298178_C=0,rs117551222_C=0,rs1345327_C=0,rs9630629_C=1,rs12922473_G=0,rs6499148_A=1,rs113717347_A=0,rs12924810_A=0,rs78831406_T=0,rs7499206_C=0,rs8051957_G=1,rs8050297_A=1

# sh /home/wletsou/scripts/haplotype_rpart.sh haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls.chr16.52679189-53075000.txt "rs3112586_T=1,rs4238751_G=1,rs80176221_G=0,rs7186592_G=1,rs118136222_T=0,rs118155214_C=0,rs77700051_C=0,rs11075587_C=1,rs117106695_C=0,rs118186521_C=0,rs16951525_C=0,rs117923020_A=0,rs79501461_T=0,rs77362202_T=0,rs117155117_G=0,rs62043164_G=0,rs17361619_C=0,rs77960456_G=0,rs117544975_C=0,rs1362413_A=1,rs4784249_C=1,rs113555472_A=0,rs59793894_G=0,rs74020833_C=0,rs17296337_G=0,rs79207645_A=0,rs117947188_C=0,rs61131526_A=0,rs117553291_G=0,rs17296940_C=0,rs76028093_C=0,rs17297072_T=0,rs117746278_A=0,rs17298178_C=0,rs79993873_C=0,rs78320298_G=0,rs62041381_C=0,rs117222080_T=0,rs117551222_C=0,rs7500472_A=0,rs72810045_A=0,rs12446016_G=0,rs73591253_T=0,rs150704739_T=0,rs113964117_A=0,rs142367768_T=0,rs117911537_T=0,rs16951771_A=0,rs73597288_A=0,rs12919591_G=0,rs9630629_C=0,rs6499147_C=1,rs12922473_G=1,rs78322982_A=0,rs62043282_A=0,rs113717347_A=0,rs62043323_A=0,rs12924810_A=0,rs34261756_C=0,rs117400881_C=0,rs78831406_T=0,rs4784276_C=0,rs77419902_T=0,rs78724580_A=0,rs72812119_G=0,rs62047145_T=0,rs76365437_A=0,rs13335861_G=0,rs36099013_A=0,rs72812136_T=0,rs79815180_T=0,rs74018575_C=0,rs72812140_T=0,rs34382563_T=0,rs56085964_C=0,rs12922187_G=0,rs79721160_A=0,rs79088862_A=0,rs4784287_T=1,rs117297253_A=0,rs56278234_T=0" "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" 1 # h5
#  reduced = rs3112586_T=1,rs4238751_G=1,rs77700051_C=0,rs11075587_C=1,rs77362202_T=0,rs17361619_C=0,rs1362413_A=1,rs117553291_G=0,rs17296940_C=0,rs17298178_C=0,rs117222080_T=0,rs12446016_G=0,rs142367768_T=0,rs12919591_G=0,rs9630629_C=0,rs12922473_G=1,rs62043282_A=0,rs113717347_A=0,rs12924810_A=0,rs117400881_C=0,rs62047145_T=0,rs76365437_A=0,rs13335861_G=0,rs74018575_C=0,rs34382563_T=0,rs56085964_C=0,rs12922187_G=0,rs79088862_A=0,rs56278234_T=0

# dbGaP replication

# combined haplotype file for cases and controls
# awk 'BEGIN{OFS="\t"} {print $0}' dbgap28544_cases.indiv dbgap28544_controls.indiv > dbgap28544_cases+dbgap28544_controls.indiv
# awk 'NR==FNR{id[$1]; next} (NR!=FNR && FNR==1){print $0} (NR!=FNR && FNR>1 && $1 in id){print $0}' dbgap28544_cases+dbgap28544_controls.indiv haplotype_estimates.ukbb_bca_cases+ukbb_bca_controls+dbgap28544_cases+dbgap28544_controls.chr16.52679189-53075000.txt > haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52679189-53075000.txt

# effect of GWAS hit
# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh  /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52679189-53075000.txt "rs4784227_T=1" "" "" "" "" rs4784227 /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts
# LRT p value = 1.808780e-50 on 1 degree(s) of freedom with OR = 1.229432e+00 (1.196436e+00 to 1.263338e+00), frequency = 2.644816e-01 (29276) = 2.822479e-01/2.433550e-01 (16971/12305)

# bsub -P SJLIFE -J dbgap28544_haplotype_haplotype_LD2 -oo dbgap28544_haplotype_haplotype_LD2.out -eo dbgap28544_haplotype_haplotype_LD2.err -R "rusage[mem=512]" "sh /home/wletsou/scripts/ukbb_haplotype_haplotype_LD2.sh haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52679189-53075000.txt \"rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0\" \"rs4784227_T=1\" /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2" # h1 (reduced) vs GWAS hit
# 0.227641	0.264482	0.227641
# rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0	rs4784227_T=1	0.81965	1


# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52679189-53075000.txt "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "" "" "" "" dbgap28544.h1.reduced /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts # h1 reduced
# LRT p value = 8.520011e-43 on 1 degree(s) of freedom with OR = 1.219175e+00 (1.184992e+00 to 1.254345e+00), frequency = 2.276407e-01 (25198) = 2.430814e-01/2.092793e-01 (14616/10582)

# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52679189-53075000.txt "rs12324961_C=1,rs3112586_T=1,rs4238751_G=1,rs80176221_G=0,rs7186592_G=1,rs118136222_T=0,rs118155214_C=0,rs77700051_C=0,rs117106695_C=0,rs118186521_C=0,rs16951525_C=0,rs117923020_A=0,rs79501461_T=0,rs77362202_T=0,rs117155117_G=0,rs62043164_G=0,rs17361619_C=0,rs77960456_G=0,rs117544975_C=0,rs1362413_A=0,rs59793894_G=0,rs74020833_C=0,rs17296337_G=0,rs79207645_A=0,rs117947188_C=0,rs61131526_A=0,rs117553291_G=0,rs17296940_C=0,rs76028093_C=0,rs17297072_T=0,rs79291013_G=0,rs117746278_A=0,rs17298178_C=0,rs79993873_C=0,rs78320298_G=0,rs62041381_C=0,rs117551222_C=0,rs7500472_A=0,rs72810045_A=0,rs12446016_G=0,rs150704739_T=0,rs142367768_T=0,rs4783791_G=0,rs117911537_T=0,rs1344484_C=0,rs16951771_A=0,rs73597288_A=0,rs6499147_C=1,rs12922473_G=0,rs78322982_A=0,rs6499148_A=1,rs62043282_A=0,rs113717347_A=0,rs62043323_A=0,rs12924810_A=0,rs80356178_A=0,rs34261756_C=0,rs117400881_C=0,rs1477029_A=1,rs7199322_G=0,rs4784276_C=0,rs77419902_T=0,rs78724580_A=0,rs72812119_G=0,rs72812122_A=0,rs76365437_A=0,rs13335861_G=0,rs4352043_C=1,rs36099013_A=0,rs72812136_T=0,rs79815180_T=0,rs74018575_C=0,rs34382563_T=0,rs12598778_C=0,rs79721160_A=0,rs117297253_A=0,rs3931698_T=1" "1" "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "1" "" dbgap28544.h2.full-length /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts # h2 full-length
# LRT p value = 8.714920e-01 on 1 degree(s) of freedom with OR = 1.484867e+00 (1.329547e-01 to 1.658331e+01), frequency = 2.710223e-05 (3) = 3.326237e-05/1.977692e-05 (2/1)

# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52679189-53075000.txt "rs12324961_C=1,rs3112586_T=1,rs7186592_G=1,rs16951525_C=0,rs117923020_A=0,rs79501461_T=0,rs1362413_A=0,rs59793894_G=0,rs17296337_G=0,rs79207645_A=0,rs117553291_G=0,rs79291013_G=0,rs117746278_A=0,rs150704739_T=0,rs1344484_C=0,rs12922473_G=0,rs6499148_A=1,rs117400881_C=0,rs1477029_A=1,rs77419902_T=0,rs72812119_G=0,rs79815180_T=0" "1" "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "1" "" dbgap28544.h2.reduced /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts # h2 reduced
# LRT p value = 4.022201e-01 on 1 degree(s) of freedom with OR = 5.699137e-01 (9.387118e-02 to 3.460079e+00), frequency = 4.517038e-05 (5) = 3.326237e-05/5.933075e-05 (2/3)

# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52679189-53075000.txt "rs12324961_C=0,rs3112586_T=1,rs12927162_G=0,rs118136222_T=0,rs118155214_C=0,rs77700051_C=0,rs117106695_C=0,rs118186521_C=0,rs16951525_C=0,rs117923020_A=0,rs79501461_T=0,rs117155117_G=0,rs62043164_G=0,rs17361619_C=0,rs77960456_G=0,rs113555472_A=0,rs59793894_G=0,rs4784253_G=0,rs74020833_C=0,rs17296337_G=0,rs4377151_G=0,rs79207645_A=0,rs117947188_C=0,rs61131526_A=0,rs117553291_G=0,rs17296940_C=0,rs76028093_C=0,rs17297072_T=0,rs79291013_G=0,rs8044632_C=0,rs117746278_A=0,rs9939896_T=0,rs79993873_C=0,rs78320298_G=0,rs62041381_C=0,rs117222080_T=0,rs117551222_C=0,rs7500472_A=0,rs72810045_A=0,rs12446016_G=0,rs73591253_T=0,rs113964117_A=0,rs142367768_T=0,rs4783791_G=0,rs117911537_T=0,rs1344484_C=0,rs73597288_A=0,rs12919591_G=0,rs9630629_C=0,rs78322982_A=0,rs6499148_A=1,rs62043282_A=0,rs113717347_A=0,rs62043323_A=0,rs12924810_A=0,rs80356178_A=0,rs34261756_C=0,rs117400881_C=0,rs1477029_A=1,rs7199322_G=0,rs4784276_C=0,rs77419902_T=0,rs78724580_A=0,rs72812119_G=0,rs72812122_A=0,rs62047145_T=0,rs76365437_A=0,rs13335861_G=0,rs36099013_A=0,rs72812136_T=0,rs74018575_C=0,rs34382563_T=0,rs12598778_C=0,rs56085964_C=1,rs7499206_C=1,rs12922187_G=0,rs8051957_G=1,rs8050297_A=1,rs79721160_A=0,rs79088862_A=0,rs4784287_T=1,rs117297253_A=0,rs3931698_T=1" "1" "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "1" "" dbgap28544.h3.full-length /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts # h3 full-length
# LRT p value = 2.894044e-01 on 1 degree(s) of freedom with OR = 7.015253e-01 (2.522935e-01 to 1.950655e+00), frequency = 1.355111e-04 (15) = 1.164183e-04/1.582153e-04 (7/8)

# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52679189-53075000.txt "rs12324961_C=0,rs3112586_T=1,rs12927162_G=0,rs79501461_T=0,rs4784253_G=0,rs17296337_G=0,rs4377151_G=0,rs79207645_A=0,rs117553291_G=0,rs17297072_T=0,rs79291013_G=0,rs78320298_G=0,rs72810045_A=0,rs73591253_T=0,rs1344484_C=0,rs12919591_G=0,rs9630629_C=0,rs78322982_A=0,rs80356178_A=0,rs117400881_C=0,rs1477029_A=1,rs77419902_T=0,rs72812119_G=0,rs56085964_C=1,rs12922187_G=0,rs8051957_G=1,rs8050297_A=1,rs4784287_T=1" "1" "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "1" "" dbgap28544.h3.reduced /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts # h3 reduced
# LRT p value = 7.799733e-02 on 1 degree(s) of freedom with OR = 5.887576e-01 (2.599748e-01 to 1.333343e+00), frequency = 2.168178e-04 (24) = 1.663119e-04/2.768768e-04 (10/14)

# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52679189-53075000.txt "rs12324961_C=0,rs3112586_T=1,rs12927162_G=0,rs80176221_G=0,rs118136222_T=0,rs118155214_C=0,rs77700051_C=0,rs11075587_C=0,rs117106695_C=0,rs118186521_C=0,rs16951525_C=0,rs117923020_A=0,rs79501461_T=0,rs77362202_T=0,rs117155117_G=0,rs62043164_G=0,rs17361619_C=0,rs77960456_G=0,rs117544975_C=0,rs1362413_A=0,rs113555472_A=0,rs59793894_G=0,rs4784253_G=0,rs74020833_C=0,rs117947188_C=0,rs61131526_A=0,rs117553291_G=0,rs17296940_C=0,rs76028093_C=0,rs17297072_T=0,rs79291013_G=0,rs117746278_A=0,rs17298178_C=0,rs79993873_C=0,rs78320298_G=0,rs62041381_C=0,rs117551222_C=0,rs7500472_A=0,rs72810045_A=0,rs12446016_G=0,rs73591253_T=0,rs150704739_T=0,rs113964117_A=0,rs1345327_C=0,rs142367768_T=0,rs4783791_G=0,rs117911537_T=0,rs1344484_C=0,rs16951771_A=0,rs73597288_A=0,rs12919591_G=0,rs9630629_C=1,rs6499147_C=1,rs12922473_G=0,rs78322982_A=0,rs6499148_A=1,rs62043282_A=0,rs113717347_A=0,rs62043323_A=0,rs12924810_A=0,rs80356178_A=0,rs34261756_C=0,rs117400881_C=0,rs7199322_G=0,rs78831406_T=0,rs4784276_C=0,rs78724580_A=0,rs72812119_G=0,rs72812122_A=0,rs76365437_A=0,rs13335861_G=0,rs4352043_C=1,rs36099013_A=0,rs72812136_T=0,rs79815180_T=0,rs72812140_T=0,rs34382563_T=0,rs12598778_C=0,rs56085964_C=0,rs7499206_C=0,rs8051957_G=1,rs8050297_A=1,rs79088862_A=0,rs117297253_A=0,rs3931698_T=1" "1" "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "1" "" dbgap28544.h4.full-length /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts # h4 full-length
# LRT p value = 3.002981e-01 on 1 degree(s) of freedom with OR = 7.730751e-01 (3.264932e-01 to 1.830498e+00), frequency = 1.897156e-04 (21) = 1.663119e-04/2.175461e-04 (10/11)

# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52679189-53075000.txt "rs12324961_C=0,rs3112586_T=1,rs80176221_G=0,rs11075587_C=0,rs16951525_C=0,rs117923020_A=0,rs77362202_T=0,rs117544975_C=0,rs113555472_A=0,rs4784253_G=0,rs17297072_T=0,rs17298178_C=0,rs117551222_C=0,rs1345327_C=0,rs9630629_C=1,rs12922473_G=0,rs6499148_A=1,rs113717347_A=0,rs12924810_A=0,rs78831406_T=0,rs7499206_C=0,rs8051957_G=1,rs8050297_A=1" "1" "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "1" "" dbgap28544.h4.reduced /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts # h4 reduced
# LRT p value = 5.125934e-01 on 1 degree(s) of freedom with OR = 9.249488e-01 (4.056053e-01 to 2.109268e+00), frequency = 2.077838e-04 (23) = 1.995742e-04/2.175461e-04 (12/11)

# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52679189-53075000.txt "rs3112586_T=1,rs4238751_G=1,rs80176221_G=0,rs7186592_G=1,rs118136222_T=0,rs118155214_C=0,rs77700051_C=0,rs11075587_C=1,rs117106695_C=0,rs118186521_C=0,rs16951525_C=0,rs117923020_A=0,rs79501461_T=0,rs77362202_T=0,rs117155117_G=0,rs62043164_G=0,rs17361619_C=0,rs77960456_G=0,rs117544975_C=0,rs1362413_A=1,rs4784249_C=1,rs113555472_A=0,rs59793894_G=0,rs74020833_C=0,rs17296337_G=0,rs79207645_A=0,rs117947188_C=0,rs61131526_A=0,rs117553291_G=0,rs17296940_C=0,rs76028093_C=0,rs17297072_T=0,rs117746278_A=0,rs17298178_C=0,rs79993873_C=0,rs78320298_G=0,rs62041381_C=0,rs117222080_T=0,rs117551222_C=0,rs7500472_A=0,rs72810045_A=0,rs12446016_G=0,rs73591253_T=0,rs150704739_T=0,rs113964117_A=0,rs142367768_T=0,rs117911537_T=0,rs16951771_A=0,rs73597288_A=0,rs12919591_G=0,rs9630629_C=0,rs6499147_C=1,rs12922473_G=1,rs78322982_A=0,rs62043282_A=0,rs113717347_A=0,rs62043323_A=0,rs12924810_A=0,rs34261756_C=0,rs117400881_C=0,rs78831406_T=0,rs4784276_C=0,rs77419902_T=0,rs78724580_A=0,rs72812119_G=0,rs62047145_T=0,rs76365437_A=0,rs13335861_G=0,rs36099013_A=0,rs72812136_T=0,rs79815180_T=0,rs74018575_C=0,rs72812140_T=0,rs34382563_T=0,rs56085964_C=0,rs12922187_G=0,rs79721160_A=0,rs79088862_A=0,rs4784287_T=1,rs117297253_A=0,rs56278234_T=0" "1" "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "1" "" dbgap28544.h5.full-length /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts # h5 full-length
# LRT p value = 7.284873e-02 on 1 degree(s) of freedom with OR = 2.163203e+00 (1.130407e+00 to 4.139616e+00), frequency = 4.065334e-04 (45) = 5.321980e-04/2.570999e-04 (32/13)

# sh /home/wletsou/scripts/dbgap28544_haplotype_model9_count_fit.sh /research/rgs01/project_space/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/for_william/phenotype/dbgap28544.pheno.txt haplotype_estimates.dbgap28544_cases+dbgap28544_controls.chr16.52679189-53075000.txt "rs3112586_T=1,rs4238751_G=1,rs77700051_C=0,rs11075587_C=1,rs77362202_T=0,rs17361619_C=0,rs1362413_A=1,rs117553291_G=0,rs17296940_C=0,rs17298178_C=0,rs117222080_T=0,rs12446016_G=0,rs142367768_T=0,rs12919591_G=0,rs9630629_C=0,rs12922473_G=1,rs62043282_A=0,rs113717347_A=0,rs12924810_A=0,rs117400881_C=0,rs62047145_T=0,rs76365437_A=0,rs13335861_G=0,rs74018575_C=0,rs34382563_T=0,rs56085964_C=0,rs12922187_G=0,rs79088862_A=0,rs56278234_T=0" "1" "rs78338509_G=0,rs12447430_C=0,rs45577538_T=0,rs8046985_A=1,rs45542333_T=0,rs78268044_T=0,rs45454402_G=0,rs45512493_G=0,rs34750829_A=0,rs9931232_A=1,rs4594251_T=0,rs4784227_T=1,rs45584434_A=0,rs45560737_A=0,rs45575339_T=0,rs45587544_A=0,rs78841172_C=0,rs76000465_C=0,rs111925335_C=0,rs75127968_A=0" "1" "" dbgap28544.h5.reduced /scratch_space/wletsou/sjlife/GWAS/UKBB_chr16.1/phase2 /home/wletsou/scripts # h5 reduced
# LRT p value = 1.223821e-01 on 1 degree(s) of freedom with OR = 1.980284e+00 (1.051638e+00 to 3.728969e+00), frequency = 4.155675e-04 (46) = 5.321980e-04/2.768768e-04 (32/14)
