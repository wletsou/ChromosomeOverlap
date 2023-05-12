# for randomly permuting cases and controls such that the number of hetero- and homozygotes is preserved

library(optparse)

option_list = list(make_option(c("-f","--file"),action = "store",type = "character",help = "File of subject IDs"),
                   make_option(c("-n","--n_samples"),action = "store",type = "numeric",help = "Number of samples to draw as cases"),
                   make_option(c("-e","--no_header"),action = "store_false",type = "logical",help = "Indicates file has no header",default = TRUE),
                   make_option(c("-r","--replacement",action = "store_true",type = "logical",help = "Whether to sample with replacement"),default = FALSE),
                   make_option(c("-s","--seed"),action = "store",type = "numeric",help = "Random seed for sampling",default = 20200116),
                   make_option(c("-o","--output"),action = "store",type = "character",help = "Output file name",default = "row_sample.txt"),
                   make_option(c("-d","--directory"),action = "store",type = "character",help = "Output directory",default=getwd()))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

file.in <- opt$file
n <- opt$n_samples
header <- opt$no_header
repl <- opt$replacement
seed <- opt$seed
directory <- opt$directory
output <- opt$output

df <- read.table(file.in,header = header)
set.seed(seed)
df.sample <- df[sample(1:nrow(df),n,replace = repl),] # sample n subjects to be cases
write.table(df.sample,sprintf("%s/%s",directory,output),row.names = FALSE,col.names = FALSE,quote = FALSE)
