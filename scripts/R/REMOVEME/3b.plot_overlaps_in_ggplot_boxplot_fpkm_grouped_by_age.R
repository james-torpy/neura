#load packages:
library(reshape2)
library(ggplot2)

#set up directory structure:
projectname = "neura/"
samplename = "stanley_pc"
inType = ".overlaps"

homeDir = "/home/jamtor/"
projectDir = paste0(homeDir, "projects/", projectname)
resultsDir = paste0(projectDir, "results/")
rawDir = paste0(projectDir, "/raw_files/")
inDir = (paste0(resultsDir, samplename, inType))
inDir

setwd(inDir)
load("sample_overlaps_fpkm.Rdata")

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

#fetch patient information from excel file and create two dataframes sorted by age
patient_metadata = read.table(file = paste0(rawDir, "NeuropathologyConsortiumDemog.txt"), header = TRUE, sep="\t", colClasses = "character")
patient_metadata = data.frame(patient_metadata$Sample, patient_metadata$Age)
colnames(patient_metadata) = c("sample", "age")
patient_metadata_age_sorted = patient_metadata[order(patient_metadata$age),]

#subset data according to phenotype:
#sample_overlaps_control = subset(patient_metadata, patient_metadata$phenotype == "Unaffected control")
#sample_overlaps_schizophrenia = subset(patient_metadata, patient_metadata$phenotype == "Schiz.")
#sample_overlaps_bipolar = subset(patient_metadata, patient_metadata$phenotype == "BP")
#sample_overlaps_depression = subset(patient_metadata, patient_metadata$phenotype == "Dep.")

#add an insignificant arbitary value to each overlap count so data can be put on log scale if needed:
#sample_overlaps_for_log = sample_overlaps
#sample_overlaps_for_log$Overlap_count_per_million_reads = sample_overlaps$Overlap_count_per_million_reads + 1

#create boxplot of sample overlap counts with Transcript_IDs on x axis, order sorted by median for each transcript:
pdf("/home/jamtor/projects/neura/results/stanley_pc.overlaps/sample_overlaps_boxplot_fpkm.pdf", height = 14, width = 20)
sample_overlaps_boxplot = ggplot(sample_overlaps, aes(x=reorder(Transcript_ID, -Median), y=Overlap_count_per_million_reads)) + geom_boxplot (aes(fill=Zero_counts_count)
#) + scale_y_log10(
) + scale_fill_continuous(name="Number of zero counts \n per gene transcript"
) + ggtitle("Expression of genes linked to schizophrenia in schizophrenia patients"
) + xlab("Gene transcript ID"
) + ylab("Overlap of reads to transcripts (RPKM)"
) + theme(plot.title = element_text(face = "bold", size = 16
), axis.text.x = element_text(angle=90, vjust=0.6))
sample_overlaps_boxplot
dev.off()
