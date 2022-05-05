args = commandArgs(trailingOnly = TRUE) 
print(args)
cat("\n")
options("width"=300)
library(data.table)
library(stringr)
library(ggplot2)
library(plotly)
library(scales)
for (i in 1:length(args)) { #check that inputs are of the form x = y
  if (regexpr('^[A-Za-z0-9_,.]+=(.*)',args[i],perl=TRUE)[1] == -1) {
    stop(sprintf('Input %s must be of the form x = y',i))
  }
  eval(parse(text = stringr::str_trim(gsub('=+(.+)',"=\'\\1\'",args[i]))))
}

save.widget <- function(plot,filename) {
  out <- tryCatch(htmlwidgets::saveWidget(as_widget(plot),filename),error = function(e) {NULL},warning = function(w) {NULL})
  return(out)
}

require(scales);
rev_log10_trans <- function() {
  scales::trans_new(
    name = "rev_log10", 
    transform = function(x) -log10(x), 
    inverse = function(x) 10^(-x));
}

rev_log10_breaks <- function(x) {
  -log10(x)
}

# directory <- "Documents/Overlap paper"
if (!exists("directory")) {
  directory <- gsub("(.*)/[^/]*","\\1",file,perl = TRUE) # trim after the last /
  if (directory==file) { # gsub returns "file" if there is no / in "file"
    directory <- getwd()
  }
}

# file="/Volumes/wletsou/sjlife/GWAS/UKBB_chr11.11/fisher_exact.ukbb_bca_cases.Iteration000-007.translated.txt"
# file="/Volumes/wletsou/sjlife/GWAS/UKBB_chr11.11/cox_hr.ukbb_bca_cases.Iteration000-007.translated.txt"

filename <- gsub("(.*).txt","\\1.html",gsub(".*/([^/]*)","\\1",file,perl = TRUE),perl = TRUE) # trim file extension and before last /, append .html extension

DT <- as.data.table(read.table(file,header = FALSE))
setnames(DT,colnames(DT),c("haplotype","iteration","cases","controls","OR","p_value"))

if (nrow(DT)==1) {
  DT <- rbind(DT,data.table(haplotype="null",iteration=DT$iteration,cases=0,controls=0,OR=1,p_value=1))
}

x_limits <- c(min(DT$p_value),max(DT$p_value))
y_limits <- c(min(DT$OR),max(DT$OR))
c_limits <- c(min(DT$iteration),max(DT$iteration))

p <- ggplot(DT,aes(x=p_value,y=OR,color=iteration,text = paste("cases: ",cases,"\ncontrols: ",controls))) + geom_point(size=3) + scale_x_continuous(trans='rev_log10',breaks = trans_breaks(function(x) -log10(x), function(x) 10 ^ (-x))(x_limits)) + labs(x=sprintf("Haplotype p value"),y=sprintf("Haplotype OR")) + scale_color_continuous(name="Iteration",breaks = seq(max(1,c_limits[2]),c_limits[1],by=-1)) + theme_gray(base_size=22) 

p <- ggplotly(p,tooltip = c("text","OR","p_value","iteration")) 

save.widget(p,sprintf("%s/%s",directory,filename))
cat(sprintf("Wrote file %s/%s\n",directory,filename))
