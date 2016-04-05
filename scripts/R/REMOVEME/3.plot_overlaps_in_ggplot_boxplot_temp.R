#load packages:
library(reshape2)
library(ggplot2)

#set up directory structure:
projectname = "neura/"
samplename = "stanley_pc"
inType = ".overlaps"

homeDir = "/Users/jamestorpy/clusterHome/"
projectDir = paste0(homeDir, "projects/", projectname)
resultsDir = paste0(projectDir, "results/")
inDir = (paste0(resultsDir, samplename, inType))
inDir

setwd(inDir)
load("sample_overlaps.Rdata")

i=1
medians = data.frame(Transcript_ID = character(0), median = numeric(0), stringsAsFactors = F)
rownames=rownames(sample_overlaps)
for (Transcript_ID in rownames) {
  numeric_values = as.numeric(sample_overlaps[Transcript_ID, ])
  median = median(numeric_values)
  print(paste0("The median of ", Transcript_ID, " is ", median))
  writeLines("\n")
  medians[i, 1] = Transcript_ID
  medians[i, 2] = median
  i = i+1
}

medians = medians[order(medians[ ,2]),]

#make an additional column of dataframe for sample names and collapse sample name-specific columns:
sample_overlaps = melt(t(sample_overlaps))
#name columns:
names(sample_overlaps) = c("Sample_Name", "Transcript_ID", "Overlap_Count_Per_Million_Reads")
#rearrange columns in better order:
sample_overlaps = sample_overlaps[c("Transcript_ID", "Sample_Name", "Overlap_Count_Per_Million_Reads")]
#merge 'medians' and 'sample_overlaps' dataframes by Transcript_ID:
sample_overlaps = merge(sample_overlaps, medians, by="Transcript_ID")
#add an insignificant arbitary value to each overlap count so data can be put on log scale:
sample_overlaps_for_log = sample_overlaps
sample_overlaps_for_log$Overlap_Count_Per_Million_Reads = sample_overlaps$Overlap_Count_Per_Million_Reads + 0.000001
#create boxplot of sample overlap counts with Transcript_IDs on x axis, order sorted by median for each transcript:
sample_overlaps_boxplot = ggplot(sample_overlaps_for_log, aes(x=reorder(Transcript_ID, -median), y=Overlap_Count_Per_Million_Reads)) + geom_boxplot (
) + scale_y_log10(
) + theme(axis.text.x = element_text(angle=90, vjust=0.6))
sample_overlaps_boxplot
