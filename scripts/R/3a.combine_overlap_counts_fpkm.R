
#This R script loads output from script 1b_1.overlap_count_from_bam_fpkm.R for each RNA-seq dataset and 
#combines the RPKM counts into a single dataframe

#set up directory structure
projectname = "neura"
#if fetching all files from project (not just one sample type), leave
#'samplename' and 'in_dataType' blank (i.e. "")
samplename = ""
in_dataType = ""

#for shell, use following homeDir:
#homeDir = "/home/jamtor/"

#for RStudio, use following homeDir:
homeDir = "/Users/jamestorpy/clusterHome/"

projectDir = paste0(homeDir, "projects/", projectname)
resultsDir = paste0(projectDir, "/results/")

#setup inPath/outFile:
inPath = resultsDir
outFile = paste0(inPath, projectname, "_overlaps_fpkm.Rdata")

writeLines("\n The inPath is:")
cat(inPath)
writeLines("\n")
writeLines("\n The outFile is:")
cat(outFile)
writeLines("\n")


all_files = ( list.files(path = inPath, pattern="rpkms.Rdata", full.names = TRUE, recursive = TRUE) )
sample_ids = basename(all_files)
#replace 'rpkms.Rdata' with '_rpkms' for sample IDs in 'sample_ids':
sample_ids = gsub("rpkms.Rdata", "_rpkms", sample_ids)

writeLines("The files used are:")
cat(sample_ids)
writeLines("\n")

#load all files in the 'all_files' list and put the dataframes of each file in the list 'sample_overlaps':
sample_overlaps = lapply(all_files, function(x) {
  load(file = x)
  cat(".")
  get(ls()[ls()!= "filename"])
}
)

names(sample_overlaps) = basename(all_files)

writeLines("These are the sample overlaps:")
print(sample_overlaps)

save.image(outFile)
