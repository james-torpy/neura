
#This is an R script that takes the combined overlap rpkms from RNA-seq sets in '2a.combine_overlap_counts_fpkm.R'
#and puts them into a boxplot

#load packages:
library(reshape2)
library(ggplot2)

#set up directory structure:
projectname = "neura"
samplename = ""
inType = ""

homeDir = "/home/jamtor/"
projectDir = paste0(homeDir, "projects/", projectname, "/")
resultsDir = paste0(projectDir, "results/")
inDir = resultsDir
cat(inDir)

setwd(inDir)
load("whole_brain_overlap_fpkm.Rdata")
load("neura_overlaps_fpkm.Rdata")

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
sample_overlaps_for_log$Overlap_count_per_million_reads = sample_overlaps$Overlap_count_per_million_reads + 0.000001

#create boxplot of sample overlap counts with Transcript_IDs on x axis, order sorted by median for each transcript:
sample_overlaps_boxplot = ggplot(sample_overlaps_for_log, aes(x=reorder(Transcript_ID, -Median), y=Overlap_count_per_million_reads)) + geom_boxplot (aes(fill=Zero_counts_count)
) + scale_y_log10(
) + scale_fill_continuous(name="Number of zero counts \n per gene transcript"
) + ggtitle("Expression of genes linked to schizophrenia in schizophrenia patients"
) + xlab("Overlap of reads to gene transcripts"
) + ylab("Gene transcript ID"
) + theme(plot.title = element_text(face = "bold", size = 16
), axis.text.x = element_text(angle=90, vjust=0.6
#) + geom_line(data = brain_counts
), legend.position(5))
sample_overlaps_boxplot
