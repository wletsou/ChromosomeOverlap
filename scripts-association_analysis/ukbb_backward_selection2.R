# for getting the LRT p-values of each haplotype in the combined model, obtained by removing each haplotype from the model one-at-a-time

# Rscript ukbb_backward_selection2.R file=h2+h3.ukbb_discovery.new_allele_counts.txt n=2 groups=1,2,3 threshold=0.05 n_covs=10

args = commandArgs(trailingOnly = TRUE) 

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

if (exists("n_covs")) { # number non-haplotype columns on the right-hand side of the table
  n_covs <- as.numeric(n_covs)
} else {
  n_covs <- 0 # none in addition to age and affect status
}

if (exists("verbose")) { # whether to print intermediate output
  verbose <- as.logical(as.numeric(verbose))
} else {
  verbose <- TRUE
}

DT <- as.data.table(read.table(file,header = TRUE))
vars <- colnames(DT)[2:(ncol(DT)-2-n_covs)] # haplotypes, eg. 00,01,10,11; excludes age and affected
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

if (exists("n")) { # drop the last n variables from model 1 in comparisson to model 0
  n <- as.numeric(n)
} else {
  n <- 1 # drop the last variable if none specified
}

if (exists("groups")) {
  groups <- strsplit(groups,split = ",") # ranges of positions to be grouped, e.g., 1,2-4,5
} else {
  groups <- seq(1,length(vars))
}

if (exists("threshold")) { # p-value threshold for dropping variables from model
  threshold <- as.numeric(threshold)
} else {
  threshold <- 0.05 # drop variables with p >= 0.05
}

ll<-as.numeric(unlist(lapply(groups,function(X) gsub("([0-9]*)-.*","\\1",X,perl=TRUE)))) # lower limit of positions to be grouped
ul<-as.numeric(unlist(lapply(groups,function(X) gsub(".*-([0-9]*)","\\1",X,perl=TRUE)))) # upper limit of positions to be grouped
if (n<length(ul)) { # cases when there are INCLUDED_HAPLOTYPES in the model
  ranges <- lapply(1:(length(ll)-n),function(X) which(unlist(lapply(pos[1:(ul[length(ul)-n]+1)],function(Y) any(Y >= ll[X] & Y <= ul[X]) ) ) ) ) # indices in pos list where "1-positions" within the limits ll and ul exist; the 1-positions of variables in the last n groups are not considered
  ranges <- c(ranges,lapply(length(ll)-(n-1):0,function(X) which(unlist(lapply(pos,function(Y) any(Y >= ll[X] & Y <= ul[X]) ) ) ) ) ) # add on the 1-positions of the variables in the last n groups
} else { # cases when there are only TEST_HAPLOTYPES in the model
  ranges <- lapply(1:length(ll),function(X) which(unlist(lapply(pos,function(Y) any(Y >= ll[X] & Y <= ul[X]) ) ) ) )
}
ranges <- ranges[unlist(lapply(ranges,length)) > 0] # remove groups with no representation among columns of DT

if (sum(unlist(lapply(ranges,length))==0)==0) { # check that each variable category is represented
  if (verbose) {
    print(sprintf("DT[,h%d:=%s]",0,paste(vars[[which(unlist(lapply(pos,function(X) any(X < 0))))]],collapse="+")))
  }
  eval(parse(text = sprintf("DT[,h%d:=%s]",0,paste(vars[[which(unlist(lapply(pos,function(X) any(X < 0))))]],collapse="+")) ) ) # null haplotype
  str0 <- "" # baseline model including all haplotypes separately
  remaining_vars <- c() # c("h1")
  var_list <- c() # list of all haplotypes
  for (i in seq(1,length(ranges))) {
    if (verbose) {
      print(sprintf("DT[,h%d:=%s]",i,paste(vars[ranges[[i]]],collapse="+")))
    }
    eval(parse(text = sprintf("DT[,h%d:=%s]",i,paste(vars[ranges[[i]]],collapse="+")))) # add haplotype h_i to table/model
    str0 <- sprintf("%s%sh%d",str0,ifelse(nchar(str0)>0," + ",""),i)
    if (i > length(ranges) - n) { # only list variables corresponding to TEST_HAPLOTYPES
      remaining_vars <- c(remaining_vars,sprintf("h%d",i))
    }
    var_list <- c(var_list,sprintf("h%d",i)) # all variables, including INCLUDED_HAPLOTYPES
  }
  length_total <- 0 # number of merged haplotypes currently in ranges_combined 

  counts_separate <- DT[,lapply(.SD,sum),.SDcols=seq(ncol(DT)-length_total-length(ranges),ncol(DT)),by=affected][order(affected),] # counts in each group
  freq_separate <- DT[,lapply(.SD,function(X) sum(X) / (2 * length(X))),.SDcols=seq(ncol(DT)-length_total-length(ranges)+0,ncol(DT)),by=affected][order(affected),]
  counts_combined <- DT[,lapply(.SD,sum),.SDcols=seq(ncol(DT)-length_total-length(ranges)+0,ncol(DT))] # counts combined
  freq_combined <- DT[,lapply(.SD,function(X) sum(X) / (2 * length(X))),.SDcols=seq(ncol(DT)-length_total-length(ranges)+0,ncol(DT))]
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
  
  max_p_value <- 1
  LRT_pvals <- data.table(remaining_vars=remaining_vars,pval=1)
  while (length(remaining_vars)>0 & max_p_value >= threshold) {
    
    cox_formula0 <- sprintf("Surv(age,affected) ~ %s%s%s",str0,ifelse(nchar(str0)>0 && n_covs>0," + ",""),covars_str) # baseline model after excluding the most non-significant haplotype
    if (verbose) {
      cat(sprintf("%s\n",cox_formula0))
    }
    cox0 <- eval.cox(DT,cox_formula0) # evaluate baseline model for this round
    if (length(cox0) > 0) {
      if (verbose) {
        print(summary(cox0))
      }
      if (length(ranges) > 1) { # more than one parameter in model 1 (below), i.e., more than 0 parameter in model 0
        deviance0 <- 2 * diff(cox0$loglik) # log likelihood of fitted model 1, 2 * (L_model - L_null)
        df0 <- nrow(summary(cox0)$coefficients)
      } else {
        deviance0 <- 0
        df0 <- 0
      }
      out_summary0<-data.table(summary(cox0)$coefficients,keep.rownames="var")
    }

    str <- str0
    for (i in remaining_vars) {
      str <- gsub(sprintf("(?<=^)%s[' + ']*|(?<=[' + ']) %s([' + '])+|[' + ']*%s(?=$)",i,i,i),"\\1",str,perl = TRUE) # delete all remaining new variables i in model string
    }
    b <- unlist(strsplit(str,split="[ ]*[+][ ]*")) # split baseline string into variables
    for (i in remaining_vars) {
      str <- ""
      a <- gsub("h([0-9]*)","\\1",i,perl=TRUE) # number corresponding to remaining variable i
      combined <- lapply(b,function(X) c(unlist(strsplit(gsub("h([0-9.]*)","\\1",X,perl=TRUE),"[.]")),a)) # join number of remaining variable to numbers of variable in baseline model
      for (j in seq_along(b)) {
        if (verbose) {
          cat(sprintf("\nDT[,h%s := %s + %s]\n",paste(combined[[j]],collapse="."),b[[j]],i))
        }
        eval(parse(text = sprintf("DT[,h%s := %s + %s]",paste(combined[[j]],collapse="."),b[[j]],i))) # combine remaining variable i with each baseline variable
        length_total <- length_total+1 # added one more variable to table
        str <- sprintf("%s%s%s",str,ifelse(nchar(str)>0," + ",""),paste(c("h",paste(combined[[j]],collapse=".")),collapse=""))
      }
      
      for (j in remaining_vars) {
        str <- sprintf("%s%s%s",str,ifelse(nchar(str)>0," + ",""),j)
      }
      str <- gsub(sprintf("(?<=^)%s[' + ']*|(?<=[' + ']) %s([' + '])+|[' + ']*%s(?=$)",i,i,i),"\\1",str,perl = TRUE) # delete variable i in model string
      cox_formula1 <- sprintf("Surv(age,affected) ~ %s%s%s",str,ifelse(nchar(str)>0 && n_covs>0," + ",""),covars_str) # include interaction in model
      if (verbose) {
        cat(sprintf("\n%s\n",cox_formula1))
      }
      cox1 <- eval.cox(DT,cox_formula1)
      if (length(cox1) > 0) {
        if (verbose) {
          print(summary(cox1))
        }
        deviance1 <- 2 * diff(cox1$loglik) # log likelihood of fitted model 1, 2 * (L_model - L_null)
        df1 <- nrow(summary(cox1)$coefficients)
        out_summary<-data.table(summary(cox1)$coefficients,keep.rownames="var")
      }
      
      if (exists("deviance0") & exists("deviance1") & exists("df0") & exists("df1")) {
        p <- pchisq(abs(deviance1 - deviance0),abs(df1 - df0),lower.tail = FALSE)
        LRT_pvals[remaining_vars==i,]$pval<-p
        
        # Get stats for dropped variable using model0
        HR <- exp(out_summary0[var==i,coef])
        HR_ll <- exp(out_summary0[var==i,coef]+qnorm(0.025)*out_summary0[var==i,`se(coef)`])
        HR_ul <- exp(out_summary0[var==i,coef]+qnorm(0.975)*out_summary0[var==i,`se(coef)`])
        cat(sprintf("%s: LRT p value = %0.6e on %d degree(s) of freedom with HR = %0.6e (%0.6e to %0.6e), frequency = %0.6e (%d) = %0.6e/%0.6e (%d/%d)\n",i,p,abs(df1-df0),HR,HR_ll,HR_ul,eval(parse(text = sprintf("freq_combined$%s[1]",i))),eval(parse(text = sprintf("counts_combined$%s[1]",i))),eval(parse(text = sprintf("freq_separate$%s[2]",i))),eval(parse(text = sprintf("freq_separate$%s[1]",i))),eval(parse(text = sprintf("counts_separate$%s[2]",i))),eval(parse(text = sprintf("counts_separate$%s[1]",i)))))
      }
    }
    
    if (nrow(LRT_pvals)>0) { # get p-value for least-significant haplotype
      max_p_value <- LRT_pvals[pval==max(pval),]$pval
    }
    
    if (exists("dropped_hap")) {
      rm(dropped_hap) # clear last-removed haplotype
    }
    dropped_hap <- LRT_pvals[pval==max(pval) & pval >= threshold,]$remaining_vars # variable corresponding to maximum LRT p-value above threshold

    if (length(dropped_hap)>0) {
      for (i in remaining_vars) {
        str0 <- gsub(sprintf("(?<=^)%s[' + ']*|(?<=[' + ']) %s([' + '])+|[' + ']*%s(?=$)",i,i,i),"\\1",str0,perl = TRUE) # delete all remaining new variables i in model string
      }
      a <- gsub("h([0-9]*)","\\1",dropped_hap,perl=TRUE)
      b <- unlist(strsplit(str0,split="[ ]*[+][ ]*")) # split baseline string into variables
      remaining_vars <- setdiff(LRT_pvals$remaining_vars,dropped_hap) # remove dropped variable
      combined <- lapply(b,function(X) c(unlist(strsplit(gsub("h([0-9.]*)","\\1",X,perl=TRUE),"[.]")),a)) # join number of remaining variable to numbers of variable in baseline model
      str0 <- ""
      for (j in seq_along(b)) {
        str0 <- sprintf("%s%s%s",str0,ifelse(nchar(str0)>0," + ",""),paste(c("h",paste(combined[[j]],collapse=".")),collapse="")) # add combined variables to new baseline string
      }
      for (j in remaining_vars) {
        str0 <- sprintf("%s%s%s",str0,ifelse(nchar(str0)>0," + ",""),j)
      }
      
      counts_separate <- DT[,lapply(.SD,sum),.SDcols=seq(ncol(DT)-length_total-length(ranges),ncol(DT)),by=affected][order(affected),] # counts in each group
      freq_separate <- DT[,lapply(.SD,function(X) sum(X) / (2 * length(X))),.SDcols=seq(ncol(DT)-length_total-length(ranges)+0,ncol(DT)),by=affected][order(affected),]
      counts_combined <- DT[,lapply(.SD,sum),.SDcols=seq(ncol(DT)-length_total-length(ranges)+0,ncol(DT))] # counts combined
      freq_combined <- DT[,lapply(.SD,function(X) sum(X) / (2 * length(X))),.SDcols=seq(ncol(DT)-length_total-length(ranges)+0,ncol(DT))]
      if (verbose) { # Print haplotype counts and frequencies
        cat("\n")
        print(t(counts_separate)) # counts in each group
        cat("\n")
        print(t(freq_separate))
        cat("\n")
        print(t(counts_combined)) # counts combined
        cat("\n")
        print(t(freq_combined))
        cat("\n")
      }
      LRT_pvals <- LRT_pvals[!(remaining_vars==dropped_hap),]
    }
  }
  cox_formula0 <- sprintf("Surv(age,affected) ~ %s%s%s",str0,ifelse(nchar(str0)>0 && n_covs>0," + ",""),covars_str) # exclude interaction from model
  if (verbose) {
    cat(sprintf("\nNew baseline model:\n%s\n",cox_formula0))
  }
  cox0 <- eval.cox(DT,cox_formula0)
  if (length(cox0) > 0) {
    if (verbose) {
      print(summary(cox0))
    }
    if (length(ranges) > 1) { # more than one parameter in model 1 (below), i.e., more than 0 parameter in model 0
      deviance0 <- 2 * diff(cox0$loglik) # log likelihood of fitted model 1, -2 * (L_model - L_null)
      df0 <- nrow(summary(cox0)$coefficients)
    } else {
      deviance0 <- 0
      df0 <- 0
    }
    out_summary0<-data.table(summary(cox0)$coefficients,keep.rownames="var")
  }
  if (nrow(LRT_pvals)>0) {
    print(LRT_pvals)
  }
}

