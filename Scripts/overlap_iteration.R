library(optparse)
library(parallel)

#https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/
option_list = list(make_option(c("-l","--lower"),type="numeric",help="Index lower limit"), make_option(c("-u","--upper"),type="numeric",help="Index upper limit"), make_option(c("-z","--step_size"),type="numeric",help="Index step size"), make_option(c("-n","--n"),type="numeric",help="Total number of chromosomes"), make_option(c("-s","--string"),type="character",help="awk string to be substituted into"), make_option(c("-o","--sigma"),type="numeric",help="Number of chromosomes in tuple"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

ll <- as.numeric(opt$lower)
ul <- as.numeric(opt$upper)
n <- as.numeric(opt$n)
str <- opt$string
sigma <- as.numeric(opt$sigma)
step_size <- as.numeric(opt$step_size)
print(opt)
cat("\n")

index2combo2 <- function(I,n,sigma) {
  allow.repeats <- FALSE
  Ind_array <- array(I,dim=c(1,sigma+1)) #initate an array of the linear positions of the multi-index in sigma, sigma-1, ..., 1 dimensional tables
  multiindex <- array(0,dim=c(1,sigma+1))
  
  if (sigma == 1) {
    multiindex[2] <- I + 1 #multiindex corresponds to linear index when choosing a single element per draw
  } else {
    for (i in seq(1,sigma,by = 1)) { #loop through elements of the multi-index i1, i2, ... , i_sigma
      ind_temp <- 0
      k <- 0
      while (ind_temp <= 0) {
        #subtract the k layers of the sigma-i+1-dimensional table not containing the multiindex
        ind_temp <- (choose(n - multiindex[i],sigma-i+1) - Ind_array[i]) - choose(n-(multiindex[i] + k),sigma-i+1)
        k <- k + 1
      }
      #print(c(i,k,ind_temp))
      Ind_array[,(i+1):(sigma+1)] <- choose(n-(multiindex[i] + k-1),sigma-(i+1)+1) - ind_temp  #value of the linear index in a table of one fewer dimensions
      multiindex[i+1] <- multiindex[i] + k - 1 #record the layer k at which the linear index occurs, in sigma-i+1 dimensions
    }
  }
  if (allow.repeats) {
    multiindex <- multiindex - seq_along(multiindex) + 1 + 1
  }
  return(multiindex[,2:(sigma+1)]) #print the tuple, starting from 0, corresponding to the multi-index 1,2,...,sigma
}

# generate list of (chr_ll,chr_ul,index_start,index_end)
param_list <- lapply(seq(1,ceiling((ul-ll)/step_size)),function(X) c(index2combo2(ifelse(ll+(X-1)*step_size>ul,ul,ll+(X-1)*step_size),n,sigma),ifelse(ll+(X-1)*step_size>ul,ul,ll+(X-1)*step_size),ifelse(ll+X*step_size-1>ul,ul,ll+X*step_size-1)))
print(param_list)

print(lapply(param_list,function(X) as.list(c(rep(c(X[length(X)-1],rep(c(X[-c(length(X)-1,length(X))]),each=sigma)[-sigma*(length(X)-2)],X[length(X)]),2),ceiling(X[length(X)]/step_size))))[[1]])

# awk scripts
awk_list <- lapply(param_list,function(X) do.call(sprintf,c(str,as.list(c(X[length(X)-1],rep(c(X[-c(length(X)-1,length(X))]),each=sigma)[-sigma*(length(X)-2)],X[length(X)],ceiling(X[length(X)]/step_size)))))) # https://stackoverflow.com/questions/52121200/use-sprintf-with-a-vector-rather-than-a-variable-number-of-arguments-in-r
print(awk_list)
cat(sprintf("%d core%s\n",detectCores(),ifelse(detectCores()==1,"","s")))
# parallel evaluation
system.time(mclapply(awk_list,system,mc.cores=detectCores()))
traceback()