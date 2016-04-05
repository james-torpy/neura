#!/bin/bash

#Link star output bam files from one directory to another directory with a unique sample name
#for each link

#define number of cluster cores used:
numcores=1

#make directory hierachy:
inPath="/share/Temp/jamtor/projects/neura/results/stanley_pc.star/"
outPath="/home/jamtor/projects/neura/results/bam_files"

logDir="/share/Temp/jamtor/projects/neura/scripts/logs"

mkdir -p $logDir

#fetch directory names of bam files, uniqueIDs for each file, and define where file links 
#will be outputted to:
for inDir in $inPath/*; do
	inFile="$inDir/Aligned.sortedByCoord.out.bam";
	uniqueID=`basename $inDir`;
	outFile="$outPath/$uniqueID.Aligned.sortedByCoord.out.bam"

	echo -e
	echo This is the inFile:
	echo $inFile
	echo -e
	echo This is the uniqueID:
	echo $uniqueID
	echo -e
	echo This is the outFile:
	echo $outFile
	echo -e

	link_line="ln -s $inFile $outFile"

	echo -e
	echo This is the link_line:
	echo $link_line
	echo -e

#submit to cluster:
	qsub -N LINK_$uniqueID -b y -wd $logDir -j y -R y -pe smp $numcores -V $link_line
done;