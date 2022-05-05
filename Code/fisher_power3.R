library(optparse)
library(statmod)
#https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/
option_list = list(make_option(c("-p","--p"),type="numeric",help="Frequency in cases"), make_option(c("-n","--n"),type="integer",help="Total chromosomes in cases"), make_option(c("-q","--q"),type="numeric",help="Frequency in controls"), make_option(c("-m","--m"),type="integer",help="Total chromosomes in controls"), make_option(c("-o","--OR"),type="numeric",help="Odds ratio to detect"), make_option(c("-a","--alpha"),type="numeric",help="False positive level",default=0.05), make_option(c("-i","--indiv"),action="store_true",help="Compute individual-level frequncies"), make_option(c("-c","--cases"),action="store_true",help="Use cases frequency as reference"),make_option(c("-s","--seed"),type="numeric",help="Random seed"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

p1 <- opt$p
n1 <- opt$n
p2 <- opt$q
n2 <- opt$m
OR <- opt$OR
alpha <- opt$alpha
indiv <- opt$indiv
cases <- opt$cases
seed <- opt$seed
if (length(indiv)==0) {
  indiv <- FALSE
}
if (length(cases)==0) {
  cases <- FALSE
}

if (indiv) { # number of individuals, not chromosomes
  n1 <- n1/2
  n2 <- n2/2
  p1 <- p1^2 + 2*p1*(1-p1^2) # fraction of homo- or heterozygotes
  p2 <- p2^2 + 2*p2*(1-p2^2) # fraction of homo- or heterozygotes
}

# Observed OR
if (indiv) {
  ft <- fisher.test(matrix(c(p1*n1,p2*n2,(1-p1)*n1,(1-p2)*n2),nrow=2)) # observed OR, subjects model
} else {
  ft <- fisher.test(matrix(c(x1,x2,n1-x1,n2-x2),nrow=2)) # observed OR, chromosome model
}

# Power to detect expected OR
if (cases) {
  p2 <- p1 / ( (1-p1) * OR + p1) # expected controls frequency
} else {
  p1 <- OR * p2 /(1 - p2 + OR * p2) # expected cases frequency (default)
}

z <- log(OR)/sqrt(1/n1/p1/(1-p1)+1/n2/p2/(1-p2))

if (length(seed)>0) {
  set.seed(seed)
}
pwr0 <- power.fisher.test(p1,p2,n1,n2,alpha = alpha,nsim=10000) # expected power
# pwr0 <- 1
if (z<0) {
  pwr <- pnorm(-z+qnorm(alpha/2))
} else {
  pwr <- 1-pnorm(-z+qnorm(1-alpha/2))
}

pwr <- pnorm(-z+qnorm(alpha/2)) + 1-pnorm(-z+qnorm(1-alpha/2))

cat(sprintf("OR = %0.6e (Observed %0.6e, p = %0.6e), power = %0.6e (%0.6e)",OR,ft$estimate,ft$p.value,pwr,pwr0),"\n")
# cat(sprintf("power = %0.6e (%0.6e)",pwr,pwr0))