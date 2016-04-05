#load packages:
library(reshape2)
library(ggplot2)

#set up directory structure:
projectname = "neura"
samplename = ""
inType = ""

#for shell, use following homeDir:
#homeDir = "/home/jamtor/"

#for RStudio, use following homeDir:
homeDir = "/Users/jamestorpy/clusterHome/"

projectDir = paste0(homeDir, "projects/", projectname, "/")
resultsDir = paste0(projectDir, "results/")
inDir = (paste0(resultsDir, samplename, inType))
inDir

setwd(inDir)
load("whole_brain_overlap_fpkm.Rdata")
load("neura_overlaps_fpkm.Rdata")

#fetch transcript names, patient names as vectors:
sample_overlaps_patients = names(sample_overlaps)
sample_overlaps_transcripts = sample_overlaps$C1rpkms.Rdata$genes

#make column 1 of sample_overlaps into the rownames, make column 2 into column 1:
for (transcript in sample_overlaps_patients) {
  rownames(sample_overlaps[[transcript]]) = sample_overlaps[[transcript]][ ,1]
  sample_overlaps[[transcript]][ ,1] = sample_overlaps[[transcript]][ ,2]
  colnames(sample_overlaps[[transcript]]) = c("gene_rpkmsPerTranscript", "")
  sample_overlaps[[transcript]] = subset(sample_overlaps[[transcript]], select="gene_rpkmsPerTranscript")
}

#change the current format of sample_overlaps (list of dataframes) to one dataframe, retaining column and rownames:
sample_overlaps = do.call("cbind", sample_overlaps)
colnames(sample_overlaps) = sample_overlaps_patients

i=1
medians = data.frame(Transcript_ID = character(0), Median = numeric(0), stringsAsFactors = F)
for (Transcript_ID in sample_overlaps_transcripts) {
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

###change original overlap counts per million to FPKM!
colnames(sample_overlaps) = c("Transcript_ID", "Sample_name", "Overlap_FPKM", "Median", "Zero_counts_count")
###

#add an insignificant arbitary value to each overlap count so data can be put on log scale if needed:
sample_overlaps_for_log = sample_overlaps
sample_overlaps_for_log$Overlap_FPKM = sample_overlaps$Overlap_FPKM + 1

#rename transcript ID column of brain_counts to 'Transcript_ID', counts to 'Brain_FPKM'
colnames(brain_counts) = c("Transcript_ID", "Overlap_FPKM")

brain_counts$Overlap_FPKM=as.numeric(levels(brain_counts$Overlap_FPKM))[brain_counts$Overlap_FPKM]

brain_counts_for_log = brain_counts
brain_counts_for_log$Overlap_FPKM = brain_counts$Overlap_FPKM + 1

#create boxplot of sample overlap counts with Transcript_IDs on x axis, order sorted by median for each transcript:
pdf(paste0(inDir, "pfc_boxplot_fpkm_with_capture-seq_wholebrain.pdf"), height = 17, width = 28)
sample_overlaps_boxplot = ggplot(data = sample_overlaps_for_log, aes(x=reorder(Transcript_ID, -Median), y=Overlap_FPKM)
) + scale_y_log10(
) + geom_boxplot (aes(fill=Zero_counts_count)
) + scale_fill_continuous(name="Number of zero counts \n per gene transcript"
) + geom_spoke(data = brain_counts_for_log, position = position_nudge(x = -0.38, y=0), aes(angle=0, radius=0.74, x=Transcript_ID, y=Overlap_FPKM
, colour = "Whole brain FPKMs from capture-seq")
) + scale_colour_manual(values = "red"
) + guides(colour=guide_legend(title=NULL)
) + ggtitle("Expression of novel transcripts in prefrontal cortex"
) + xlab("Gene transcript ID"
) + ylab("Gene expression / FPKM"
) + theme(plot.title = element_text(face = "bold", size = 30
), axis.title = element_text(size=20
), axis.text = element_text(size=20, angle=90, vjust=0.6
), legend.title = element_text(size=18
), legend.text = element_text(size=18
), legend.key.size = unit(1.5, "cm")
)
sample_overlaps_boxplot
dev.off()
