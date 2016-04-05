#load packages needed:
library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome)
library(Rsamtools)
library(reader)

#set up directory structure:
projectname = "neura/"
samplename = "stanley_pc/"

homeDir = "/Users/jamestorpy/clusterHome/"
projectDir = paste0(homeDir, "projects/", projectname)
resultsDir = paste0(projectDir, "results/")
inDir = paste0(resultsDir, samplename, "bam_files")

#set current working directory:
setwd(inDir)
getwd()

#assign column names for dataframe:
what = c("rname","strand","pos","qwidth")
#flag unmapped sequences to leave out:
flag = scanBamFlag(isUnmappedQuery=FALSE)
#define parameters of bam scan:
param = ScanBamParam(what=what,flag=flag)

#fetch filenames of bam files:
i=0
files=list.files()
for (file in files) {
  writeLines("\n The file used is:")
  print(file)

#index the inFile - file with your mapped library in bam format:
  indexBam(file)
  }




#index the inFile - file with your mapped library in bam format:
indexBam(cat.path(dir=inDir, )
  paste0(inDir, "case_proliferative_1/case_proliferative_1.Sorted_Aligned.out.bam")
#define the inFile:
bam=scanBam("/Users/jamestorpy/clusterHome/projects/Grant/results/4.star-cuffdiff_ercc_protocol/Grant.star/case_proliferative_1/case_proliferative_1.Sorted_Aligned.out.bam",param=param)

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

#take .gtf files of co-ordinates of genes GAPDH, MYC, p53 and make a list, in which each element consists of all the exons for all the transcripts per gene
GAPDH = import("/Users/jamestorpy/clusterHome/genomes/hg38_ercc/GAPDH_gencode.gtf")
MYC = import("/Users/jamestorpy/clusterHome/genomes/hg38_ercc/MYC_gencode.gtf")
P53 = import("/Users/jamestorpy/clusterHome/genomes/hg38_ercc/p53_gencode.gtf")

GAPDH_list = split(GAPDH, GAPDH$gene_id)
MYC_list = split(MYC, MYC$gene_id)
P53_list = split(P53, P53$gene_id)

countsPerTranscript = countOverlaps(GAPDH_list, bam_gr)

