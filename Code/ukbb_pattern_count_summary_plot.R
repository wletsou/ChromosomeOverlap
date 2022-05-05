library(ggplot2)
library(ggrepel)
library(data.table)
library(stringr)

args = commandArgs(trailingOnly = TRUE)
for (i in 1:length(args)) { #check that inputs are of the form x = y
  if (regexpr('^[A-Za-z0-9_,.]+=(.*)',args[i],perl=TRUE)[1] == -1) {
    stop(sprintf('Input %s must be of the form x = y',i))
  }
  eval(parse(text = stringr::str_trim(gsub('=+(.+)',"=\'\\1\'",args[i]))))
}

# directory <- "Documents/Overlap paper"
if (!exists("directory")) {
  directory <- gsub("(.*)/[^/]*","\\1",fasta_file,perl = TRUE) # trim after the last /
  if (directory==fasta_file) { # gsub returns "file" if there is no / in "file"
    directory <- getwd()
  }
}

# file="/Volumes/wletsou/sjlife/GWAS/UKBB_chr11.7/REP_0/Pattern_counts.ukbb_bca_20200116_cases.chr11.69242608-69428155_2,j.txt"
# file="/Volumes/wletsou/Chromosome_Overlap_results/Pattern_counts.ukbb_bca_cases.chr11.69231642-69431642_2,j.txt"
filename1 <- gsub("(.*).txt","\\1.plot_1.png",gsub(".*/([^/]*)","\\1",file,perl = TRUE),perl = TRUE)
filename2 <- gsub("(.*).txt","\\1.plot_2.png",gsub(".*/([^/]*)","\\1",file,perl = TRUE),perl = TRUE)
filename3 <- gsub("(.*).txt","\\1.plot_3.png",gsub(".*/([^/]*)","\\1",file,perl = TRUE),perl = TRUE)

DT <-  as.data.table(read.table(file,header = FALSE,fill = TRUE))
setnames(DT,colnames(DT),c("Iteration","Total","Filtered","Core"))
step_size <- 10000 # secondary axis steps
multiple <- 50000 # round secondary axis maximum to nearest "multiple"
filtered_min <- 0
filtered_max <- floor((max(DT$Filtered) + multiple)/multiple)*multiple # maximum filtered patterns, rounded to nearest "multiple"

ggplot(data = DT[-1,],aes(x = Iteration)) + geom_line(mapping = aes(y = Total),color = "blue",size = 3) + geom_line(mapping = aes(y = Filtered * 1),color = "red",size = 3) + scale_y_continuous(name = "Total patterns",sec.axis = sec_axis(~./1,name = "Filtered patterns",breaks = seq(filtered_min,filtered_max,by = step_size))) + theme(text = element_text(size=28),axis.title.y = element_text(color = "blue"),axis.title.y.right = element_text(color = "red"),plot.margin=unit(c(0.5,0.5,0.5,0.5),"in"))

ggsave(filename=sprintf("%s/%s",directory,filename1),plot = last_plot(),width = 8.5,height = 11,dpi = 300)

ggplot(data = DT,aes(x = Iteration)) + geom_line(mapping = aes(y = Core),color = "blue",size = 3) + scale_y_continuous(name = "Closed patterns") + theme(text = element_text(size=28),axis.title.y = element_text(color = "blue"),plot.margin=unit(c(0.5,0.5,0.5,0.5),"in"))

ggsave(filename=sprintf("%s/%s",directory,filename2),plot = last_plot(),width = 8.5,height = 11,dpi = 300)

ggplot(data = DT,aes(x = Iteration)) + geom_line(mapping = aes(y = Total),color = "red",size = 3) + scale_y_continuous(name = "Total patterns") + theme(text = element_text(size=28),axis.title.y = element_text(color = "red"),plot.margin=unit(c(0.5,0.5,0.5,0.5),"in"))

ggsave(filename=sprintf("%s/%s",directory,filename3),plot = last_plot(),width = 8.5,height = 11,dpi = 300)
