
#This R script loads output from script 1b_2.overlap_count_from_bam_counts_per_million_reads.R for each RNA-seq
#dataset and combines the counts into a single dataframe


#set up directory structure
projectname = "neura"
samplename = "stanley_pc"
in_dataType = ".overlaps"

homeDir = "/home/jamtor/"
projectDir = paste0(homeDir, "projects/", projectname)
resultsDir = paste0(projectDir, "/results")

#setup inDir/outFile:
inDir = paste0(resultsDir, "/", samplename, in_dataType)
outFile = paste0(inDir, "/sample_overlaps_fpkm.Rdata")

writeLines("\n The inDir is:")
inDir
writeLines("\n")
writeLines("\n The outFile is:")
outFile
writeLines("\n")


all_files = ( list.files(path = inDir, pattern="rpkms.Rdata", full.names = T) )

writeLines("The files used are:")
all_files
writeLines("\n")

sample_overlaps = lapply(all_files, function(x) {
  load(file = x)
  cat(".")
  get(ls()[ls()!= "filename"])
}
)

names(sample_overlaps) = basename(all_files)

writeLines("These are the sample overlaps:")
sample_overlaps

sample_overlaps=as.data.frame(sample_overlaps)

save.image(outFile)
