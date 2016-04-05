library(GenomicRanges)
library(ShortRead)
library(rtracklayer)

#set up directory/file structure:
projectname="neura"
gtfName = "neura_transcripts_hg38.gtf"

homeDir = "/home/jamtor/"
projectDir = paste0(homeDir, "projects/")
resultsDir = paste0(projectDir, projectname, "/results/"

gtf_homeDir = "/share/Temp/jamtor/"
gtf_projectDir = paste0(gtf_homeDir, "projects/"
transcriptsDir = paste0(gtf_projectDir, projectname, "raw_files/transcripts_files")
outFile = paste0(resultsDir,"whole_brain_overlap_fpkm.Rdata")

gtfFile = paste0(transcriptsDir, gtfName)

#load in the gene coordinates from the gtf file with included brain overlap counts:
genes = import(gtfFile)

#order genes by gene_id:
gene_list = split(genes, genes$gene_id)

#take gene_id and whole brain counts from the above and put it in a data frame:
gene_list_df = as.data.frame(c(gene_list[]))
brain_counts = unique(data.frame(gene_list_df$gene_id, gene_list_df$brain))
colnames(brain_exp) = c("gene_id", "brain_counts")

save(brain_counts, file=outFile)