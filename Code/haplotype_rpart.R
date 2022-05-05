args = commandArgs(trailingOnly = TRUE) 
# print(args)
# cat("\n")
options("width"=300)
library(data.table)
library(stringr)
library(rpart)
for (i in 1:length(args)) { #check that inputs are of the form x = y
  if (regexpr('^[A-Za-z0-9_,.]+=(.*)',args[i],perl=TRUE)[1] == -1) {
    stop(sprintf('Input %s must be of the form x = y',i))
  }
  eval(parse(text = stringr::str_trim(gsub('=+(.+)',"=\'\\1\'",args[i]))))
}

if (exists("min_split")) { # positive integer, minimum number of items in split-off bin
  min_split <- as.numeric(min_split)
} else {
  min_split <- 10
}

DT <- as.data.table(read.table(file,header = TRUE))
vars <- colnames(DT)[3:ncol(DT)] # alleles, starting in column 3
alleles <- unlist(lapply(strsplit(haplotype,split = ","),function(X) gsub("=[0-9]","",X,perl=TRUE))) # alleles in haplotype (no values)

for (allele in alleles) {if (!(allele %in% vars)) {stop(sprintf("%s not found",allele))} } # stop if some allele is not a variable in DT

# https://cran.r-project.org/web/packages/rpart/vignettes/longintro.pdf
hapstat <- factor(DT$test_haplotype,levels = 0:1,labels = c("No","Yes"))
hap_formula <- sprintf("hapstat ~ %s",paste(alleles,collapse=" + "))
print(hap_formula)
hfit <- rpart(formula = eval(parse(text = hap_formula)),data = DT,control=rpart.control(cp=-1,minsplit=min_split,method='class',maxcompete=0))
print(hfit)
reduced.alleles <- as.data.table(hfit$splits,keep.rownames='allele')[count>0,]$allele
cat("\n",paste(unlist(lapply(unlist(strsplit(haplotype,","))[-1],function(X) if (gsub("=[0-9]","",X,perl=TRUE) %in% reduced.alleles) {X} )),collapse=","),"\n")


                         