library(optparse)
library(data.table)
library(survival)

option_list = list(make_option(c("-i","--included_haplotypes_file"),type="character",help="file containing chromosomes ids and counts of included haplotypes"),
                   make_option(c("-t","--test_haplotypes_files"),type="character",help="comma-separated list of files containing chromosome ids and counts of test haplotypes joined to each included haplotype",default=""),
                   make_option(c("-p","--phenotypes_file"),type="character",help="bca phenotypes with age and genotype pcs",default=""),
                   make_option(c("-k","--haplotype_number",type="numeric",help="column number (not including row names) of test haplotype in each file")),
                   make_option(c("-g","--groups"),type="character",help="grouping of  haplotype variables as a comma-separated list with hyphens for ranges",default = ""),
                   make_option(c("-n","--n_new_haplotypes"),type = "numeric",help = "number of new haplotypes being joined to model",default = 1),
                   make_option(c("-v","--verbose"),type = "numeric",help = "whether to print model output",default = 0),
                   make_option(c("-d","--directory"),type="character",help="output directory",default=getwd()))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

included_haplotypes_file <- opt$included_haplotypes_file # file containing chromosomes ids and counts of included haplotypes
test_haplotypes_files <- unlist(strsplit(opt$test_haplotypes_files,split = ",")) # comma-separated list of files containing chromosome ids and counts of test haplotypes joined to each included haplotype
phenotypes_file <- opt$phenotypes_file # bca phenotypes with age and genotype pcs
haplotype_number <- as.numeric(opt$haplotype_number) # column number (not including row names)of test haplotype in each file
groups <- unlist(strsplit(opt$groups,split = ",")) # grouping of  haplotype variables as a comma-separated list with hyphens for ranges
n <- as.numeric(opt$n_new_haplotypes) # number of new haplotypes being joined to model
verbose <- as.logical(as.numeric(opt$verbose)) # whether to print model output
directory <- opt$directory # output directory

eval.cox <- function(X,Y) {
  out <- tryCatch(coxph(formula = eval(parse(text = Y)),data = X),error = function(e) {NULL},warning = function(w) {NULL})
  return(out)
}

df <- read.table(included_haplotypes_file,header = TRUE)
for (i in seq_len(length(test_haplotypes_files))) {
  df_temp <- read.table(test_haplotypes_files[i],header = TRUE)
  print(as.data.table(df_temp[,c(1,haplotype_number + 1)])[,sum(.SD),.SDcols=2])
  df <- merge(df,df_temp[,c(1,haplotype_number + 1)], by = 1)
}
DT <- as.data.table(df)
DT <- DT[,.(haplotypes = paste(.SD,collapse = "")),by = sid,.SDcols = 2:ncol(DT)] # concatenate haplotype counts
DT <- as.data.table(data.frame(sid = DT$sid,as.data.frame.matrix(table(DT)))) # tabulate counts of each chromosome type for each chromosome (one 1 in each row)
DT$sid <- gsub("([^_]+)_.*","\\1",DT$sid)
DT <- DT[,lapply(.SD,sum),by = sid,.SDcols = 2:ncol(DT)] # sum chromosome types over each subject's two chromosomes

phenotypes <- read.table(phenotypes_file,header = TRUE)
phenotypes <- as.data.table(phenotypes[,grepl("sid",names(phenotypes)) | grepl("age.end",names(phenotypes)) | grepl("pc",names(phenotypes)) | grepl("^diag",names(phenotypes))])
phenotypes$sid <- as.character(phenotypes$sid)

DT <- merge(DT,phenotypes,by="sid")
DT <- setNames(DT,replace(names(DT),which(names(DT) == "diag"),"affected"))
DT <- setNames(DT,replace(names(DT),which(names(DT) == "age.end"),"age"))
if (verbose) {
  print(DT)
}

vars <- colnames(DT)[grepl("X[01]+",colnames(DT))] # haplotype variables
if (length(groups) == 0) {groups <- seq_len(length(vars))} # default is no grouping
covars <- colnames(DT)[grepl("pc",colnames(DT))]
n_covs <- length(covars)
covars_str <- paste(covars,collapse = " + ") # genotype principal components

pos <- lapply(vars,function(X) unlist(gregexpr("1",X))-1) # list all positions of "1" in the variable name, indicates which haplotypes are present

ll<-as.numeric(unlist(lapply(groups,function(X) gsub("([0-9]*)-.*","\\1",X,perl=TRUE)))) # lower limit of positions to be grouped
ul<-as.numeric(unlist(lapply(groups,function(X) gsub(".*-([0-9]*)","\\1",X,perl=TRUE)))) # upper limit of positions to be grouped

if (n<length(ul)) { # cases when there are INCLUDED_HAPLOTYPES in the model
  ranges <- lapply(1:(length(ll)-n),function(X) which(unlist(lapply(pos[1:(ul[length(ul)-n]+1)],function(Y) any(Y >= ll[X] & Y <= ul[X]) ) ) ) ) # indices in pos list where "1-positions" within the limits ll and ul exist; the 1-positions of variables in the last n groups are not considered; recall that pos has an index for the 0 haplotype, so start from position 2
  ranges <- c(ranges,lapply(length(ll)-(n-1):0,function(X) which(unlist(lapply(pos,function(Y) any(Y >= ll[X] & Y <= ul[X]) ) ) ) ) ) # add on the 1-positions of the variables in the last n groups
} else { # cases when there are only TEST_HAPLOTYPES in the model
  ranges <- lapply(1:length(ll),function(X) which(unlist(lapply(pos,function(Y) any(Y >= ll[X] & Y <= ul[X]) ) ) ) )
}
ranges <- ranges[unlist(lapply(ranges,length)) > 0] # remove groups with no representation among columns of DT

if (length(ll)-n>0) { # only if there are INCLUDED_HAPLOTYPES to which TEST_HAPLOTYPES can be joined
  ranges_combined <- unlist(lapply(1:(length(ll)-n),function(X) lapply(length(ll)-(n-1):0,function(Y) which(unlist(lapply(pos,function(Z) (any(Z >= ll[X] & Z <= ul[X]) & any(Z >= ll[Y] & Z <= ul[Y]) ) | (any(Z >= ll[X] & Z <= ul[X]) & !(any(Z > ul[length(ll)-n]) ) ) ) ) ) ) ),recursive=FALSE ) # find positions in pos list where each INCLUDED_HAPLOTYPE occurs with at most one other TEST_HAPLOTYPE
  ranges_combined <- ranges_combined[unlist(lapply(ranges_combined,length)) > 0]
} else {
  ranges_combined <- list()
}

# print(matrix(unlist(mapply(function(X,Y) lapply(pos,function(Z) (any(Z >= ll[X] & Z <= ul[X]) & any(Z >= ll[Y] & Z <= ul[Y]) ) | (any(Z >= ll[X] & Z <= ul[X]) & !(any(Z > ul[length(ll)-n]) ) ) ),1:(length(ll)-n),length(ll)-(n-1):0)),nrow=length(pos),byrow=FALSE) )
# print(apply(matrix(unlist(mapply(function(X,Y) lapply(pos,function(Z) (any(Z >= ll[X] & Z <= ul[X]) & any(Z >= ll[Y] & Z <= ul[Y]) ) | (any(Z >= ll[X] & Z <= ul[X]) & !(any(Z > ul[length(ll)-n]) ) ) ),1:(length(ll)-n),length(ll)-(n-1):0),2,which),nrow=length(pos),byrow=FALSE),2,which) )

if (sum(unlist(lapply(ranges,length))==0)==0) { # check that each variable category is represented
  if (verbose) {
    print(sprintf("DT[,h%d:=%s]",0,paste(vars[[which(unlist(lapply(pos,function(X) any(X < 0))))]],collapse="+")))
  }
  eval(parse(text = sprintf("DT[,h%d:=%s]",0,paste(vars[[which(unlist(lapply(pos,function(X) any(X < 0))))]],collapse="+")) ) ) # null category
  str0 <- "" # baseline string for model combining most-recent haplotype with the INCLUED_HAPLOTYPES
  if (length(ranges_combined) >= n) { # if each of n TEST_HAPLOTYPES is joined to each INCLUDED_HAPLOTYPE
    for (i in seq(n,n*(length(ll)-n),by = n)) { # get the most-recently joined haplotype on each INCLUDED_HAPLOTYPE
      if (verbose) {
        a <- i/n # which included_haplotype group the ith test_haplotype is being joined to
        b <- n+i/n # new variable formed by joining test_haplotype to included_haplotype a
        # print(sprintf("DT[,h%s:=%s]",paste(sort(unique(unlist(pos[ranges_combined[[i]]]))),collapse="."),paste(vars[ranges_combined[[i]]],collapse="+")))
        print(sprintf("DT[,h%s:=%s]",paste(c(a,b),collapse="."),paste(vars[ranges_combined[[i]]],collapse="+")))
      }
      a <- i/n # which included_haplotype group the ith test_haplotype is being joined to
      b <- length(ll)-1+i/n # new variable formed by joining test_haplotype to included_haplotype a
      eval(parse(text = sprintf("DT[,h%s:=%s]",paste(c(a,b),collapse="."),paste(vars[ranges_combined[[i]]],collapse="+"))))
      # str0 <- sprintf("%s%sh%s",str0,ifelse(nchar(str0)>0," + ",""),paste(sort(unique(unlist(pos[ranges_combined[[i]]]))),collapse="."))
      str0 <- sprintf("%s%sh%s",str0,ifelse(nchar(str0)>0," + ",""),paste(c(a,b),collapse="."))
    }
  }
  
  str <- "" # string for the model with the most-recent haplotype separate from the INCLUDED_HAPLOTYPES
  for (i in seq(1,length(ranges))) {
    if (verbose) {
      print(sprintf("DT[,h%d:=%s]",i,paste(vars[ranges[[i]]],collapse="+")))
    }
    eval(parse(text = sprintf("DT[,h%d:=%s]",i,paste(vars[ranges[[i]]],collapse="+"))))
    str <- sprintf("%s%sh%d",str,ifelse(nchar(str)>0," + ",""),i)
    if (exists("str0")==TRUE && i >= length(ranges)-(n-1) && i <= length(ranges) - 1) { # if there is a baseline model, add all TEST_HAPLOTYPES except the newest
      str0 <- sprintf("%s%sh%d",str0,ifelse(nchar(str0)>0," + ",""),i) # model combining the newest haplotype with the baseline haplotypes
    }
  }
  
  counts_separate <- DT[,lapply(.SD,sum),.SDcols=seq(ncol(DT)-length(ranges)+0,ncol(DT)),by=affected][order(affected),] # counts in each group
  freq_separate <- DT[,lapply(.SD,function(X) sum(X) / (2 * length(X))),.SDcols=seq(ncol(DT)-length(ranges)+0,ncol(DT)),by=affected][order(affected),]
  fisher.pvalue <- DT[,lapply(.SD,function(X) fisher.test(array(c(sum(X * (affected == 1)),2 * sum(affected == 1) - sum(X * (affected == 1)),sum(X * (affected == 0)),2 * sum(affected == 0) - sum(X * (affected == 0))),c(2,2)))$p.value),.SDcols=seq(ncol(DT)-length(ranges)+0,ncol(DT))]
  counts_combined <- DT[,lapply(.SD,sum),.SDcols=seq(ncol(DT)-length(ranges)+0,ncol(DT))] # counts combined
  freq_combined <- DT[,lapply(.SD,function(X) sum(X) / (2 * length(X))),.SDcols=seq(ncol(DT)-length(ranges)+0,ncol(DT))]
  if (verbose) { # Print haplotype counts and frequencies
    cat("\n")
    print(counts_separate) # counts in each group
    cat("\n")
    print(freq_separate)
    cat("\n")
    print(fisher.pvalue)
    cat("\n")
    print(counts_combined) # counts combined
    cat("\n")
    print(freq_combined)
    cat("\n")
  }
  
  if (exists("str0") == FALSE && n_covs == 0) { # null model if no included haplotypes and no covariates
    str0 <- "1"
  } else if (exists("str0") == FALSE && n_covs > 0) {
    str0 <- "" # no inlcuded haplotypes but covariates included
  }
  
  cox_formula0 <- sprintf("Surv(age,affected) ~ %s%s%s",str0,ifelse(nchar(str0)>0 && n_covs>0," + ",""),covars_str) # baseline model with newest haplotype combined with INCLUDED_HAPLOTYPES
  if (verbose) {
    cat(sprintf("%s\n",cox_formula0))
  }
  cox0 <- eval.cox(DT,cox_formula0)
  if (length(cox0) > 0) {
    if (verbose) {
      print(summary(cox0))
    }
    if (length(ranges) + n_covs > 1) { # more than one parameter in model 1 (below), i.e., more than 0 parameters in model 0
      deviance0 <- 2 * diff(cox0$loglik) # log likelihood of fitted model 1, 2 * (L_model - L_null)
      df0 <- nrow(summary(cox0)$coefficients)
    } else { # model with one fewer variables is the null model, only one value in cox0$loglik instead of 2
      deviance0 <- 0
      df0 <- 0
    }
    out_summary<-summary(cox0)
  }
  
  cox_formula1 <- sprintf("Surv(age,affected) ~ %s%s%s",str,ifelse(nchar(str)>0 && n_covs>0," + ",""),covars_str) # include interaction in model
  if (verbose) {
    cat(sprintf("%s\n",cox_formula1))
  }
  cox1 <- eval.cox(DT,cox_formula1)
  if (length(cox1) > 0) {
    if (verbose) {
      print(summary(cox1))
    }
    deviance1 <- 2 * diff(cox1$loglik) # log likelihood of fitted model 1, 2 * (L_model - L_null)
    df1 <- nrow(summary(cox1)$coefficients)
    out_summary<-summary(cox1)
  }
  idx <- which(rownames(out_summary$coefficients) == sprintf("h%d",length(ranges))) # row index of last variable
}
if (exists("deviance0") && exists("deviance1") && exists("df0") && exists("df1")) {
  p <- pchisq(abs(deviance1 - deviance0),abs(df1 - df0),lower.tail = FALSE)
  HR_ll <- exp(out_summary$coefficients[idx,1]+qnorm(0.025)*out_summary$coefficients[idx,3])
  HR_ul <- exp(out_summary$coefficients[idx,1]+qnorm(0.975)*out_summary$coefficients[idx,3])
  cat(sprintf("LRT p value = %0.6e on %d degree(s) of freedom with HR = %0.6e (%0.6e to %0.6e), frequency = %0.6e (%d) = %0.6e/%0.6e (%d/%d)\n",p,abs(df1-df0),out_summary$coefficients[idx,2],HR_ll,HR_ul,eval(parse(text = sprintf("freq_combined$h%s",length(ranges)))),eval(parse(text = sprintf("counts_combined$h%s",length(ranges)))),eval(parse(text = sprintf("freq_separate$h%s[2]",length(ranges)))),eval(parse(text = sprintf("freq_separate$h%s[1]",length(ranges)))),eval(parse(text = sprintf("counts_separate$h%s[2]",length(ranges)))),eval(parse(text = sprintf("counts_separate$h%s[1]",length(ranges))))))
} else {
  cat(sprintf("LRT p value = NA on NA degree(s) of freedom with HR = NA (NA to NA), frequency = NA (NA) = NA/NA (NA/NA)\n"))
}

  