
#load packages needed:
library(rtracklayer)

projectname = "neura"

#set up directory structure:
homeDir = "/Users/jamestorpy/clusterHome/"
projectDir = paste0(homeDir, "projects/", projectname, "/")
rawDir = paste0(projectDir, "raw_files")
inDir = (paste0(rawDir, "/annotations"))
inDir

setwd(inDir)

bedFile_hg38 <- "neura_transcripts_hg38.bed"
bed_hg38 <- import(bedFile)

gtfFile_hg19 <- "neura_transcripts_hg19.gtf"
gtf_hg19 <- import(gtfFile)

bed_hg38_df <- as.data.frame(bed_hg38)
head(bed_hg38_df)
gtf_hg19_df <- as.data.frame(gtf_hg19)
head(gtf_hg19_df)

gtf_hg38_df <- gtf_hg19_df
gtf_hg38_df$start <- bed_hg38_df$start
gtf_hg38_df$start <- bed_hg38_df$start
gtf_hg38_df$end <- bed_hg38_df$end
gtf_hg38_df$end <- bed_hg38_df$end
head(gtf_hg38_df)

gtf_hg38_gr <- makeGRangesFromDataFrame(gtf_hg38_df, keep.extra.columns = TRUE)
head(gtf_hg38_gr)

export.gff(gtf_hg38_gr, paste0(inDir, "/gtf_hg38.gff"))

#garvan_0820_coord <- subset(gtf_hg38_df, gtf_hg38_df$gene_id == "garvan_0820", select = c(gene_id, seqnames, transcript_id, start, end))

