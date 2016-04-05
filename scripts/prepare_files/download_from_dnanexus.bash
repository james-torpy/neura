#!/bin/bash

for i in {22..31}; do
	download_line_1="dx download PFC_$i\_*_1.fastq.tar.gz";
	download_line_2="dx download PFC_$i\_*_2.fastq.tar.gz";

	echo -e
	echo This is the download_line_1:
	echo $download_line_1
	echo -e
	echo This is the download_line_2
	echo $download_line_2

	j=$(($i-1))
	
	echo -e
	echo This is j:
	echo $j

	qsub -N dl1_$i -hold_jid dl1_$j -wd /share/Temp/jamtor/projects/neura/raw_files/stanley_pc -b y -j y -R y -pe smp 1 -V $download_line_1;

	qsub -N dl2_$i -hold_jid dl2_$j -wd /share/Temp/jamtor/projects/neura/raw_files/stanley_pc -b y -j y -R y -pe smp 1 -V $download_line_2;
done;
