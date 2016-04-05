#This script takes a dataframe of overlap counts from all samples of project from script 2a or 2b and mapped to
#a GTF file of transcripts, groups them by age from a tab-delimited metadata file and creates a boxplot
#showing expression of each transcript by each age group.

#Has the option to split transcripts into two groups for clarity of plot - to fetch first half # lines 89-91,
#to fetch second half # lines 85-87, to fetch all # both

#load packages:
library(reshape2)
library(ggplot2)
library(grid)

#set up directory structure:
projectname = "neura"
#if fetching all files from project (not just one sample type), leave
#'samplename' blank (i.e. "")
samplename = ""
inType = ""
inFile_name = paste0(projectname, "_overlaps_fpkm.Rdata")
###to fetch first half of transcripts group make following "_a", for second half "_b", for all transcripts "" - complete this option
transcripts_group_suffix = "_b"

homeDir = "/Users/jamestorpy/clusterHome/"
projectDir = paste0(homeDir, "projects/", projectname)
resultsDir = paste0(projectDir, "/results/")
inDir = (paste0(resultsDir, samplename, inType))
inDir
inFile_name

setwd(inDir)
load(inFile_name)

#re-set up directory structure:
projectname = "neura"
samplename = ""
inType = ""

homeDir = "/Users/jamestorpy/clusterHome/"
projectDir = paste0(homeDir, "projects/", projectname)
resultsDir = paste0(projectDir, "/results/")
rawDir = paste0(projectDir, "/raw_files/")
inDir = (paste0(resultsDir, samplename, inType))
inDir

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

###find a better way to do the following:
#split sample_overlaps into two:
#to fetch first half:
#row_number = nrow(sample_overlaps)
#split_number = (row_number/2)
#sample_overlaps_a = sample_overlaps[1:split_number, ]

#to fetch second half:
row_number = nrow(sample_overlaps)
split_number = (row_number/2) + 1
sample_overlaps_b = sample_overlaps[split_number:row_number, ]
###

transcript_ids = as.vector(unique(sample_overlaps_b$Transcript_ID))
i=1
zero_counts_df = data.frame(Transcript_ID = character(0), Zero_counts_count = numeric(0), stringsAsFactors = F)
for (t_id in transcript_ids) {
  cat(paste0("The transcript is: ", t_id))
  writeLines("\n")
  subset = subset(sample_overlaps_b, Transcript_ID == t_id)
  zero_counts_number = sum(subset$Overlap_count_per_million_reads == 0)
  cat(paste0("Number of patients with zero counts: ", zero_counts_number))
  writeLines("\n")
  writeLines("\n")
  zero_counts_df[i, "Transcript_ID"] = t_id
  zero_counts_df[i, "Zero_counts_count"] = zero_counts_number
  i = i+1
}

#merge 'zero_counts_df' and 'sample_overlaps_b' dataframes by Transcript_ID:
sample_overlaps_b = merge(sample_overlaps_b, zero_counts_df, by="Transcript_ID")

#fetch patient information from excel file and create two dataframes sorted by age
patient_metadata = read.table(file = paste0(rawDir, "NeuropathologyConsortiumDemog.txt"), header = TRUE, sep="\t", colClasses = "character")
patient_metadata = data.frame(patient_metadata$Sample, patient_metadata$Age)
colnames(patient_metadata) = c("sample", "age")
patient_metadata_age_sorted = patient_metadata[order(patient_metadata$age),]

sample_names = patient_metadata_age_sorted$sample
Sample_name = ""
i=1
for (name in sample_names) {
  print(name)
  Sample_name[i] = paste0(name, "rpkms.Rdata")
  i=i+1
}

patient_metadata_age_sorted[ ,3] = Sample_name
colnames(patient_metadata_age_sorted) = c("sample", "age", "Sample_name")

sample_overlaps_b_age_sorted=merge(sample_overlaps_b, patient_metadata_age_sorted, by="Sample_name")
sample_overlaps_b_age_sorted=sample_overlaps_b_age_sorted[order(sample_overlaps_b_age_sorted$age),]
###find a better way to do the following:
#sample_overlaps_b_age_sorted[1:58,8] = "fetal"
#sample_overlaps_b_age_sorted[59:580,8] = "25-40"
#sample_overlaps_b_age_sorted[581:1334,8] = "41-55"
#sample_overlaps_b_age_sorted[1335:1595,8] = "56-70"
#colnames(sample_overlaps_b_age_sorted) = c("Sample_name", "Transcript_ID", "Overlap_count_per_million_reads", "Median", "Zero_counts_count", "sample", "age", "age_group")
###

#data = subset(sample_overlaps_b_age_sorted, Transcript_ID == "garvan_0531")
#data = data[ 1:7]
#data

i=1
for (age in sample_overlaps_b_age_sorted$age) {
  if (age < 1) {
    sample_overlaps_b_age_sorted[i,8] = "fetal"
  } else if ( (age >= 1) & (age <= 24) ) {
    sample_overlaps_b_age_sorted[i,8] = "1-24"
  } else if ( (age >= 25) & (age <= 40) ) {
    sample_overlaps_b_age_sorted[i,8] = "25-40"
  } else if ( (age >= 41) & (age <= 55) ) {
    sample_overlaps_b_age_sorted[i,8] = "41-55"
  } else if ( (age >= 56) & (age <= 70) ) {
    sample_overlaps_b_age_sorted[i,8] = "56-70"
  } else {
    sample_overlaps_b_age_sorted[i,8] = ">70"
  }
  i=i+1
}
  
###

###put data in wide format -> narrow format, separate df for median and age group, then merge

#create boxplot of sample overlap counts with Transcript_IDs on x axis, order sorted by median for each transcript:
pdf(paste0(inDir, "/pfc_boxplot_fpkm_byage", transcripts_group_suffix, ".pdf"), height = 18, width = 34)
sample_overlaps_b_boxplot = ggplot(sample_overlaps_b_age_sorted, aes(x=reorder(Transcript_ID, -Median), y=Overlap_count_per_million_reads)
) + geom_point(aes(colour=age_group, fill=age_group), position = position_jitterdodge(jitter.width = 1.0, jitter.height = 0, dodge.width=0.75)
) + scale_colour_manual(values = c("firebrick2", "dodgerblue3", "chartreuse3", "orange")
) + geom_boxplot (aes(fill=age_group)
#) + scale_y_log10(
) + scale_fill_manual(values = c("firebrick2", "dodgerblue3", "chartreuse3", "orange")
) + ggtitle("Expression of novel transcripts in prefrontal cortex"
) + xlab("Gene transcript ID"
) + ylab("Gene expression / FPKM"
) + theme(plot.title = element_text(face = "bold", size = 40), axis.title = element_text(size=30), axis.text.x = element_text(size = 30, angle=90, vjust=0.6
), axis.text.y = element_text(size = 30), legend.key.size = unit(2.5, "cm"), legend.key.height = unit(2, "cm"), legend.title = element_text(size = 30), legend.text = element_text(size = 30))
sample_overlaps_b_boxplot
dev.off()

###use below format for ggplot in future:
#data=sample_overlaps_b_age_sorted

#dataShort=data[,c("Transcript_ID","Overlap_count_per_million_reads","age_group")]
#p<-ggplot(dataShort,aes(Transcript_ID,Overlap_count_per_million_reads),fill=age_group)
#p<-p+geom_boxplot(aes(fill=age_group))
#p



