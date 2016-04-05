#!/home/jamtor/local/lib/r/R-3.2.2/bin/

#set up directory structure:
projectname = "neura/"
samplename = "stanley_pc"
inExt = ".star/"
genomename = "hg38_ercc/"

homeDir = "/share/Temp/jamtor/"
script_homeDir = "/home/jamtor/"

projectDir = paste0(homeDir, "projects/", projectname)
resultsDir = paste0(projectDir, "results/")
inPath = paste0(resultsDir, samplename, inExt)

script_projectDir = paste0(script_homeDir,"projects/", projectname)
scriptDir = paste0(script_projectDir,"/scripts/")
logDir = paste0(scriptDir,"R/logs")

system(paste0("mkdir -p ",logDir))

#fetch working directories of inFiles:
directories = ( list.files(inPath) )

for (directory_name in directories) {
        inDir = paste0(inPath, directory_name)
        writeLines("\n")
        writeLines("The directory used is:")
        print(inDir)

        sampleID = directory_name

#fetch filenames of bam files:
        files = list.files(inDir,pattern = "Aligned.sortedByCoord.out.bam$",full.names=T)
        #files = ( files[c(TRUE, FALSE)] )
        for (inFile in files) {
                writeLines("\n The sample used is:")
                print(sampleID)

                #Submit the job to the processing script
                #system((paste0("qsub -P GenomeInformatics -b y -j y -N sampleID -wd $logDir -pe smp 2 -V R --vanilla <",scriptDir,"/R/gene_count_from_bam_processing.r --bam=",inFile," --sampleID=",sampleID)))
                writeLines("\n")
                cat(paste0("qsub -P GenomeInformatics -b y -j y -N sampleID -wd $logDir -pe smp 2 -V R --vanilla <",scriptDir,"/R/gene_count_from_bam_processing.r --bam=",inFile," --sampleID=",sampleID))
                writeLines("\n")
        }
}