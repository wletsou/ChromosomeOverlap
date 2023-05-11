# file="/Volumes/yasuigrp/EpiGenetics/common/yasuigrp_data/biobank/gds/ukb.bca.hap/ukb.bca.hap.chr11.gds"
# file="/Volumes/yasuigrp/EpiGenetics/common/yasuigrp_data/ProjHap/dbgap28544/topmed/hap2k/gds/hap2k.chr11.hg19.qced.gds"
# file="/Volumes/yasuigrp/EpiGenetics/common/yasuigrp_data/biobank/gds/ukb.bca.hap/ukb.bca.hap.chr16.gds"
# file="/Volumes/yasuigrp/EpiGenetics/common/yasuigrp_data/biobank/gds/ukb.bca.hap/ukb.bca.hap.chr10.gds"
# file="Z://ResearchHome/Groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/biobank/gds/ukb.bca.hap/ukb.bca.hap.chr21.gds"
# file="Z://ResearchHome/Groups/yasuigrp/projects/EpiGenetics/common/yasuigrp_data/biobank/gds/ukb.bca.hap/ukb.bca.hap.chr22.gds"

args = commandArgs(trailingOnly = TRUE) 
print(args)
cat("\n")
options("width"=300)
library(data.table)
library(stringr)
install.packages("BiocManager")
BiocManager::install("SeqArray")
BiocManager::install("Rsamtools")
BiocManager::install("gdsfmt")
library(BiocManager)
library(SeqArray)
library(Rsamtools)

for (i in 1:length(args)) { #check that inputs are of the form x = y
  if (regexpr('^[A-Za-z0-9_,.]{1,}(?=.*=)={1}?(?!.*=)',args[i],perl=TRUE)[1] == -1) { 
    stop(sprintf('Input %s must be of the form x = y',i)) 
  }
  eval(parse(text = stringr::str_trim(gsub('=+(.+)',"=\'\\1\'",args[i]))))
}

gds <- seqOpen(file, allow.duplicate=TRUE)
vcf.fn <- gsub(".gds$",".vcf.gz",file,perl=TRUE)
seqGDS2VCF(gds,vcf.fn,use_Rsamtools=TRUE)
seqClose(gds)

bgzip(vcf.fn,gsub(".gz",".bgz",vcf.fn)) # if vcf.fn is not bg-zipped

file.rename(gsub(".gz",".bgz",vcf.fn),vcf.fn)
