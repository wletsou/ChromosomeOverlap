library(optparse)
library(statmod)
#https://www.r-bloggers.com/2015/09/passing-arguments-to-an-r-script-from-command-lines/
option_list = list(make_option(c("-x","--x"),type="integer",help="Counts in population 1"), make_option(c("-n","--n"),type="integer",help="Total in population 1"), make_option(c("-y","--y"),type="integer",help="Counts in population 2"), make_option(c("-m","--m"),type="integer",help="Total in population 2"),make_option(c("-a","--alpha"),type="numeric",help="False positive level",default=0.05))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

x1 <- opt$x
n1 <- opt$n
x2 <- opt$y
n2 <- opt$m
alpha <- opt$alpha

p1 <- x1/n1
p2 <- x2/n2

OR <- p1*(1-p2)/p2/(1-p1)
z <- log(OR)/sqrt(1/n1/p1/(1-p1)+1/n2/p2/(1-p2))

ft <- fisher.test(matrix(c(x1,x2,n1-x1,n2-x2),nrow=2))

pwr0 <- power.fisher.test(p1,p2,n1,n2,alpha = alpha,nsim=10000)

if (z<0) {
  pwr <- pnorm(-z+qnorm(alpha/2))
} else {
  pwr <- 1-pnorm(-z+qnorm(1-alpha/2))
}

pwr <- pnorm(-z+qnorm(alpha/2)) + 1-pnorm(-z+qnorm(1-alpha/2))

cat(sprintf("OR = %0.6e, p = %0.6e, power = %0.6e (%0.6e)",ft$estimate,ft$p.value,pwr,pwr0),"\n")