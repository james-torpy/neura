#!/bin/bash

#make modules to load here

module load gi/samtools/1.2

#number of cores
numcores=15

#genome directories
genomeName="hg38_ercc"
annotationName="gencode.v24.annotation"

homeDir="/home/jamtor"
genomeDir="$homeDir/genomes/$genomeName"
genomeFile="$genomeDir/$genomeName.fa"
annotationFile="$homeDir/genomes/$genomeName/$annotationName.gtf"
outDir="$genomeDir/star_ref_76bpreads"

mkdir -p $outDir

#log directory

logDir="/home/jamtor/projects/neura/scripts/QC/logs"
mkdir -p $logDir

echo This is the genomeDir:
echo $genomeDir
echo -e
echo This is the genomeFile:
echo $genomeFile
echo -e
echo This is the annotationFile:
echo $annotationFile
echo -e
echo This is the outDir
echo $outDir
echo -e

#generate the star reference files:
star_ref_line="$homeDir/local/bin/STAR \
	--runMode genomeGenerate \
	--sjdbGTFfile $annotationFile --sjdbOverhang 74 --genomeDir $genomeDir \
	--genomeFastaFiles $genomeFile --runThreadN $numcores \
	--outFileNamePrefix $outDir"

echo This is the star_ref_line:
echo $star_ref_line

#submit job with name 'RSEM_count_$sample' to 15 cluster cores:
qsub -N STAR_ref_$genomeName -wd $logDir -b y -j y -R y -pe smp $numcores -V $star_ref_line
