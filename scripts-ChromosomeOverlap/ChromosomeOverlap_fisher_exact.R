library(optparse)

option_list = list(make_option(c("-f","--file"),type="character",help="Haplotypes counts file"), make_option(c("-o","--output"),type="character",help="Output file string",default=""), make_option(c("-d","--directory"),type="character",help="Output directory",default=getwd()))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

f_in <- opt$file # file with haplotype name, a (cases count), b (cases non-count), c (controls count), d (controls non-count)
name <- opt$output
directory <- opt$directory

if (file.exists(f_in)) {
  df <- read.table(f_in,header=FALSE,sep="\t",row.names=1)
  res <- t(apply(df,1,function(X) {out=fisher.test(matrix(c(X[1],X[2],X[3],X[4]),nrow=2)); f1=X[1]/(X[1]+X[2]); f2=X[3]/(X[3]+X[4]); OR=as.numeric(out$estimate[1]); p=out$p.value[1]; return(rbind(f1,f2,OR,p))}))
  cat(sprintf("%s/fisher_exact%s.txt\n",directory,ifelse(nchar(name)>0,sprintf(".%s",name),name)))
  write.table(res,file=sprintf("%s/fisher_exact%s.txt",directory,ifelse(nchar(name)>0,sprintf(".%s",name),name)),quote=FALSE,sep="\t",row.names=TRUE,col.names=FALSE)
}