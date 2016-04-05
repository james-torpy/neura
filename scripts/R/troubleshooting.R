load("/home/jamtor/projects/neura/results/gingeras/ENCFF001RN_troubleshooting.Rdata")

library(GenomicRanges)
library(ShortRead)
library(rtracklayer)
library("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome)
library(Rsamtools)
library(R.utils)
library(edgeR)

bam_gr=GRanges(
seqnames = results$rname[correct_chromosomes],
ranges = IRanges(start=bam[[1]]$pos[correct_chromosomes], width=bam[[1]]$qwidth[correct_chromosomes]),
strand = bam[[1]]$strand[correct_chromosomes],
seqlengths = seq_lengths
        )