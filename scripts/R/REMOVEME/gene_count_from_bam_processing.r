#!/home/jamtor/local/lib/r/R-3.2.2/bin/

#load packages needed:
library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome)
library(Rsamtools)
library(R.utils)

#read arguments from the command line
args <- R.utils::commandArgs(asValues=TRUE)
if (!is.null(args[["bam"]])){inFile = args$bam}
if (!is.null(args[["sampleID"]])){sampleID = args$sampleID}

gtfFile = "/share/Temp/jamtor/projects/neura/raw_files/neura_transcripts_hg38.gtf"

resultsDir = paste0("/home/jamtor/projects/neura/results/overlaps/")
system(paste0("mkdir -p ",resultsDir))


writeLines("\n This is the gtfFile:")
gtfFile

#load in the gene coordinates from the gtf file
genes = import(gtfFile)

#assign column names for dataframe:
what = c("rname","strand","pos","qwidth")
#flag unmapped sequences to leave out:
flag = scanBamFlag(isUnmappedQuery=FALSE)
#define parameters of bam scan:
param = ScanBamParam(what=what,flag=flag)

###need to grab inFiles from other script! index the inFile - file with your mapped library in bam format:
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

gene_list = split(genes, genes$gene_id)


gene_countsPerTranscript = countOverlaps(gene_list, bam_gr)

writeLines("\n")
print(paste0("This is the gene count"))
save(gene_countsPerTranscript,file=paste0(resultsDir,sampleID,".Rdata"))
