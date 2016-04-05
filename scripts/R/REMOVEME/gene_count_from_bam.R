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

#set genes to be counted:
gene1 = "GAPDH"
gene2 = "MYC"
gene3 = "P53"

#take .gtf files of co-ordinates of genes to count and make a list, in which 
#each element consists of all the exons for all the transcripts per gene:
        gene1_inFile = import(paste0(genomeDir, gene1, "_gencode.gtf"))
        gene2_inFile = import(paste0(genomeDir, gene2, "_gencode.gtf"))
        gene3_inFile = import(paste0(genomeDir, gene3, "_gencode.gtf"))

        writeLines("\n This is the gene1 .gtf file:")
        gene1_inFile
        writeLines("\n This is the gene2 .gtf file:")
        print(gene2_inFile)
        writeLines("\n This is the gene3 .gtf file:")
        gene3_inFile

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

        gene1_list = split(gene1_inFile, gene1_inFile$gene_id)
        gene2_list = split(gene2_inFile, gene2_inFile$gene_id)
        gene3_list = split(gene3_inFile, gene3_inFile$gene_id)

        gene1_countsPerTranscript = countOverlaps(gene1_list, bam_gr)
        gene2_countsPerTranscript = countOverlaps(gene2_list, bam_gr)
        gene3_countsPerTranscript = countOverlaps(gene3_list, bam_gr)

        writeLines("\n")
        print(paste0("This is the ", gene1, " count"))
        writeLines("\n")
        print(paste0("This is the ", gene2, " count"))
        writeLines("\n")
        print(paste0("This is the ", gene3, " count"))

        }
}

                                                                                        1,1           Top

