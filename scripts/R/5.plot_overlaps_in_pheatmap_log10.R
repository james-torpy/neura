
#This is an R script that takes the combined overlap counts per million reads from RNA-seq sets and puts the log10 of
#these values into a heatmap

#load packages:
library(reshape2)
library(pheatmap)

#set up directory structure:
projectname = "neura/"
samplename = "stanley_pc"
inType = ".overlaps"

homeDir = "/home/jamtor/"
projectDir = paste0(homeDir, "projects/", projectname)
resultsDir = paste0(projectDir, "results/")
inDir = (paste0(resultsDir, samplename, inType))
print(inDir)

setwd(inDir)
load("sample_overlaps_boxplot.Rdata")

sample_overlaps_pheatmap = acast(sample_overlaps, Transcript_ID ~ Sample_name, value.var = "Overlap_count_per_million_reads")

pdf("/home/jamtor/projects/neura/results/stanley_pc.overlaps/sample_overlaps_pheatmap.pdf", width=16, height=12, onefile=FALSE)
sample_overlaps_pheatmap = sample_overlaps_pheatmap + 1
sample_overlaps_pheatmap_log10 = log10(sample_overlaps_pheatmap)
colnames(sample_overlaps_pheatmap_log10)<-gsub(".Rdata","",colnames(sample_overlaps_pheatmap_log10))
pheatmap(sample_overlaps_pheatmap_log10, main = "Expression of transcripts linked to schizophrenia in schizophrenia patients (log10 of reads per million reads)", fontsize = 14)
dev.off()

save.image("sample_overlaps_boxplot_and_pheatmap.Rdata")

#return to scripts directory:
setwd("~/projects/neura/scripts/R")