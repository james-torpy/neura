#This script takes reads from an Rdata image containing reads from a number of different patients that
#overlapped transcripts from a GTF file. It creates a boxplot ordered by transcript with highest to lowest 
#median. The number of patients with no reads mapped to a transcript is shown by colour gradient. The plot
#is saved to a pdf file and the data is saved as an Rdata image in case debugging is needed.

#load packages:
library(reshape2)
library(ggplot2)

#set up directory structure:
projectname = "neura/"
samplename = "stanley_pc"
in_dataType = ".overlaps"
inFile = "sample_overlaps.Rdata"

homeDir = "/home/jamtor/"
projectDir = paste0(homeDir, "projects/", projectname)
resultsDir = paste0(projectDir, "results/")
inDir = (paste0(resultsDir, samplename, in_dataType))

writeLines("\n The inDir is:")
inDir

setwd(inDir)
load(inFile)

i=1
medians = data.frame(Transcript_ID = character(0), Median = numeric(0), stringsAsFactors = F)
rownames=rownames(sample_overlaps)
for (Transcript_ID in rownames) {
  numeric_values = as.numeric(sample_overlaps[Transcript_ID, ])
  median = median(numeric_values)
  print(paste0("The median of ", Transcript_ID, " is ", median))
  writeLines("\n")
  medians[i, "Transcript_ID"] = Transcript_ID
  medians[i, "Median"] = median
  i = i+1
}

medians = medians[order(medians[ ,2]),]

#make an additional column of dataframe for sample names and collapse sample name-specific columns:
sample_overlaps = melt(t(sample_overlaps))
#name columns:
names(sample_overlaps) = c("Sample_name", "Transcript_ID", "Overlap_count_per_million_reads")
#rearrange columns in better order:
sample_overlaps = sample_overlaps[c("Transcript_ID", "Sample_name", "Overlap_count_per_million_reads")]
#merge 'medians' and 'sample_overlaps' dataframes by Transcript_ID:
sample_overlaps = merge(sample_overlaps, medians, by="Transcript_ID")

transcript_ids = as.vector(unique(sample_overlaps$Transcript_ID))
i=1
zero_counts_df = data.frame(Transcript_ID = character(0), Zero_counts_count = numeric(0), stringsAsFactors = F)
for (t_id in transcript_ids) {
  cat(paste0("The transcript is: ", t_id))
  writeLines("\n")
  subset = subset(sample_overlaps, Transcript_ID == t_id)
  zero_counts_number = sum(subset$Overlap_count_per_million_reads == 0)
  cat(paste0("Number of patients with zero counts: ", zero_counts_number))
  writeLines("\n")
  writeLines("\n")
  zero_counts_df[i, "Transcript_ID"] = t_id
  zero_counts_df[i, "Zero_counts_count"] = zero_counts_number
  i = i+1
  }

#merge 'zero_counts_df' and 'sample_overlaps' dataframes by Transcript_ID:
sample_overlaps = merge(sample_overlaps, zero_counts_df, by="Transcript_ID")

#add an insignificant arbitary value to each overlap count so data can be put on log scale:
sample_overlaps_for_log = sample_overlaps
sample_overlaps_for_log$Overlap_count_per_million_reads = sample_overlaps_for_log$Overlap_count_per_million_reads + 1

#create boxplot of sample overlap counts with Transcript_IDs on x axis, order sorted by median for each transcript:
sample_overlaps_boxplot = ggplot(sample_overlaps_for_log, aes(x=reorder(Transcript_ID, -Median), y=Overlap_count_per_million_reads)) + geom_boxplot (aes(fill=Zero_counts_count)
) + scale_y_log10(
) + scale_fill_continuous(name="Number of zero counts \n per gene transcript"
) + ggtitle("Expression of genes linked to schizophrenia \n in schizophrenia patients"
) + xlab("Gene transcript ID"
) + ylab("Overlap of reads to gene transcripts (log10 of reads per million reads)"
) + theme(plot.title = element_text(face = "bold", size = 16
), axis.text.x = element_text(size = 10, angle=90, vjust=0.6))

#create pdf file of sample_overlaps_boxplot:
pdf("sample_overlaps_boxplot.pdf", width=16, height=12)
sample_overlaps_boxplot
dev.off()

#save data as image in case debugging is needed
save.image("sample_overlaps_boxplot.Rdata")

#return to scripts directory:
setwd("~/projects/neura/scripts/R")
