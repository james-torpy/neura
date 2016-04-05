#!/bin/bash

workDir="/share/Temp/jamtor/projects/neura//share/Temp/jamtor/projects/neura/raw_files/stanley_pc"

files=( $(ls *.tar.gz) )

echo -e
echo These are the file:
echo ${files[@]}

for file in ${files[@]}; do
	untar_line="tar -xvf $file";
	
	echo -e
	echo This is the untar_line:
	echo $untar_line

	 qsub -N untar_$file -wd $workDir -b y -j y -R y -pe smp 1 -V $untar_line;
done;
