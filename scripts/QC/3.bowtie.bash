 #!/bin/bash

#This script builds a bowtie index from a human genome gtf file, for use with tophat

module load pethum/bowtie2/prebuilt/2.2.6
module load pethum/tophat/prebuilt/2.1.0
module load gi/samtools/1.2

#number of cores
numcores=24

#genome directories
genomeName="hg38_ercc"

homeDir="/home/jamtor/"
genomeDir="$homeDir/genomes/$genomeName"
genomeFile="$genomeDir/$genomeName.fa"

#log directory
projectname="neura"

projectDir="$homeDir/projects/$projectname"
scriptsPath="$projectDir/scripts/QC"
logDir="$scriptsPath/logs"

mkdir -p $logDir

echo This is the genomeDir:
echo $genomeDir
echo -e
echo This is the genomeFile:
echo $genomeFile
echo -e
echo This is the logDir:
echo $logDir
echo -e

#generate the bowtie 2 index files:
bowtie_line="bowtie2-build $genomeFile hg38_ercc"

echo This is the bowtie_line:
echo $bowtie_line

#submit job with name 'BOWTIE_build_$genomeName' to 10 cluster cores:
qsub -N BOWTIE_build_$genomeName -wd $logDir -b y -cwd -j y -R y -P GenomeInformatics -pe smp $numcores -V $bowtie_line

#move log files to logDir:
#move_line=mv *$genomeDir/BOWTIE_build_$genomeName.o* $logDir

#remove uneeded log files:
#remove_line=rm *$genomeDir/BOWTIE_build_$genomeName.po*

#submit above lines to cluster, holding until bowtie build job is done:
#qsub -N MOVE_logs -wd $genomeDir -b y -cwd -j y -R y -pe smp 1 -V $move_line

#qsub -N REMOVE_logs -wd $genomeDir-b y -cwd -j y -R y -pe smp 1 -V $remove_line


