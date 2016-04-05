#!/home/jamtor/local/lib/r/R-3.2.2/bin/

#load packages needed:
library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome)
library(Rsamtools)

#set up directory structure:
projectname = "neura/"
samplename = "stanley_pc"
inExt = ".star/"
genomename = "hg38_ercc/"

homeDir = "/share/Temp/jamtor/"
projectDir = paste0(homeDir, "projects/", projectname)
resultsDir = paste0(projectDir, "results/")
inPath = paste0(resultsDir, samplename, inExt)

genomeDir = paste0("/home/jamtor/genomes/", genomename)

#annotation
gtfFile <- paste0(annotationDir,"transcript.gtf")
#load in the gene coordinates from the gtf file
genes<-import(gtfFile)

#gtfDir = paste0(resultsDir,"gtfs")
#system(paste0("mkdir -p ",gtfDir))
#set genes to be counted:

#geneList<-c("GAPDH","MYC","P53")
#for(gene in geneList){
#	#execute a grep command in bash through a system command
#	system(paste0("grep ",gene," ",genomeDir,"/gencode_v24.gtf > ",gtfDir,"/",gene,".gtf"))
#}


#assign column names for dataframe:
what = c("rname","strand","pos","qwidth")
#flag unmapped sequences to leave out:
flag = scanBamFlag(isUnmappedQuery=FALSE)
#define parameters of bam scan:
param = ScanBamParam(what=what,flag=flag)


#fetch working directories of inFiles:
directories = ( list.files(inPath) )

for (directory_name in directories) {
    inDir = paste0(inPath, directory_name)
    writeLines("\n The directory used is:")
    print(inDir)

	#set current working directory:
	setwd(inDir)

	#fetch filenames of bam files:
    files = list.files(pattern = "Aligned.sortedByCoord.out.bam")
    files = ( files[c(TRUE, FALSE)] )
    for (inFile in files) {
        writeLines("\n The file used is:")
        print(inFile)
		#index the inFile - file with your mapped library in bam format:
        indexBam(inFile)
		#define the inFile:
        bam=scanBam(inFile,param=param)

		#convert the inFile to GRanges object:
        results=bam[[1]]
        correct_chromosomes=!grepl("_",bam[[1]]$rname)
        seq_lengths=seqlengths(Hsapiens)[!grepl("_",names(seqlengths(Hsapiens)))]
        bam_gr=GRanges(
                seqnames = results$rname[correct_chromosomes],
                ranges = IRanges(start=bam[[1]]$pos[correct_chromosomes], width=bam[[1]]$qwidth[correct_chromosomes]),
                strand = bam[[1]]$strand[correct_chromosomes],
                seqlengths = seq_lengths
        )



		#split the gtf file by gene id
		geneL = split(inFile, inFile$gene_id)
			gene_countsPerTranscript = countOverlaps(geneL, bam_gr)
	        cat(paste0("This is the ", gene, " count"))
	        cat(gene_countsPerTranscript)
		}
    }
}

                                                                                        1,1           Top

