args = commandArgs(trailingOnly = TRUE) 
# print(args)
# cat("\n")
options("width"=300)
library(data.table)
library(stringr)
library(survival)
for (i in 1:length(args)) { #check that inputs are of the form x = y
  if (regexpr('^[A-Za-z0-9_,.]+=(.*)',args[i],perl=TRUE)[1] == -1) {
    stop(sprintf('Input %s must be of the form x = y',i))
  }
  eval(parse(text = stringr::str_trim(gsub('=+(.+)',"=\'\\1\'",args[i]))))
}

eval.cox <- function(X,Y) {
  out <- tryCatch(coxph(formula = eval(parse(text = Y)),data = X),error = function(e) {NULL},warning = function(w) {NULL})
  return(out)
}

# file="/Volumes/wletsou/Chromosome_Overlap_results/UKBB_chr11.19/new_allele_counts_haplotypes.chr11.69231642-69431642.txt"

if (exists("n_covs")) { # number non-haplotype columns on the right-hand side of the table
  n_covs <- as.numeric(n_covs)
} else {
  n_covs <- 0 # none in addition to age and affected status
}

if (exists("verbose")) { # whether to print intermediate output
  verbose <- as.logical(as.numeric(verbose))
} else {
  verbose <- TRUE
}

DT <- as.data.table(read.table(file,header = TRUE))
vars <- colnames(DT)[2:(ncol(DT)-2-n_covs)] # haplotypes, eg. 00,01,10,11; excludes age and affected and covs (pcs)
if (verbose) {
  print(vars)
}
if (n_covs>0) { # string of covariates for model
  covars <- colnames(DT)[ncol(DT)-(2+n_covs-1):2]
  covars_str <- paste(covars,collapse=" + ")
} else {
  covars_str <- ""
}
pos <- lapply(vars,function(X) unlist(gregexpr("1",X))-1) # list all positions of "1" in the variable name, indicates which haplotypes are present

if (exists("n")) { # number of new haplotypes in the model (TEST_HAPLOTYPES)
  n <- as.numeric(n)
} else {
  n <- 1 # only the last variable is new
}

if (exists("groups")) {
  groups <- strsplit(groups,split = ",") # ranges of positions to be grouped, e.g., 1,2-4,5
} else {
  groups <- seq(1,length(vars))
}

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
    for (i in seq(n,length(ranges_combined),by = n)) { # get the most-recently joined haplotype on each INCLUDED_HAPLOTYPE
      if (verbose) {
        a <- i/n # which included_haplotype group the ith test_haplotype is being joined to
        b <- length(ll)-1+i/n # new variable formed by joining test_haplotype to included_haplotype a
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
  counts_combined <- DT[,lapply(.SD,sum),.SDcols=seq(ncol(DT)-length(ranges)+0,ncol(DT))] # counts combined
  freq_combined <- DT[,lapply(.SD,function(X) sum(X) / (2 * length(X))),.SDcols=seq(ncol(DT)-length(ranges)+0,ncol(DT))]
  if (verbose) { # Print haplotype counts and frequencies
    cat("\n")
    print(counts_separate) # counts in each group
    cat("\n")
    print(freq_separate)
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
      deviance0 <- 2 * diff(cox0$loglik) # log likelihood of fitted model 1, -2 * (L_model - L_null)
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
    deviance1 <- 2 * diff(cox1$loglik) # log likelihood of fitted model 1, -2 * (L_model - L_null)
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
