library(gpart)

pedfile <- "/Volumes/wletsou/Chromosome_Overlap_results/UKBB_chr10.1/phase2.1/LD_blocks/1000gp.10.120850000-121950000.hg38.ped"
mapfile <- "/Volumes/wletsou/Chromosome_Overlap_results/UKBB_chr10.1/phase2.1/LD_blocks/1000gp.10.120850000-121950000.hg38.map"
filename <- "/Volumes/wletsou/Chromosome_Overlap_results/UKBB_chr10.1/phase2.1/LD_blocks/1000gp.10.120850000-121950000.hg38.r2.heatmap"
output <- "/Volumes/wletsou/Chromosome_Overlap_results/UKBB_chr10.1/phase2.1/LD_blocks/1000gp.10.120850000-121950000.hg38.res.txt"

filename="/Volumes/wletsou/Chromosome_Overlap_results/UKBB_chr11.34/phase2/h2_Dprime_heatmap"
pedfile <- "/Volumes/wletsou/Chromosome_Overlap_results/UKBB_chr11.34/phase2/LD_blocks/1000gp.11.69050000-69700000.hg38.ped"
mapfile <- "/Volumes/wletsou/Chromosome_Overlap_results/UKBB_chr11.34/phase2/LD_blocks/1000gp.11.69050000-69700000.hg38.map"
filename <- "/Volumes/wletsou/Chromosome_Overlap_results/UKBB_chr11.34/phase2/LD_blocks/1000gp.11.69050000-69700000.hg38.r2.heatmap"
output <- "/Volumes/wletsou/Chromosome_Overlap_results/UKBB_chr11.34/phase2/LD_blocks/1000gp.11.69050000-69700000.hg38.res.txt"

pedfile <- "/Volumes/wletsou/Chromosome_Overlap_results/UKBB_chr16.2/phase2.1/LD_blocks/1000gp.16.52000000-53100000.hg38.ped"
mapfile <- "/Volumes/wletsou/Chromosome_Overlap_results/UKBB_chr16.2/phase2.1/LD_blocks/1000gp.16.52000000-53100000.hg38.map"
filename <- "/Volumes/wletsou/Chromosome_Overlap_results/UKBB_chr16.2/phase2.1/LD_blocks/1000gp.16.52000000-53100000.hg38.r2.heatmap"
output <- "/Volumes/wletsou/Chromosome_Overlap_results/UKBB_chr16.2/phase2.1/LD_blocks/1000gp.16.52000000-53100000.hg38.res.txt"

res2 <- BigLD(genofile=pedfile,SNPinfofile=mapfile,MAFcut=0.05,appendRare=FALSE)
LDblockHeatmap(genofile=pedfile,SNPinfofile=mapfile,blockresult=res2,assembly="GRCh37",geneDB="ucsc",LD="r2",filename=filename)
write.table(res2,file=output,quote=FALSE,col.names=TRUE,row.names=FALSE)



