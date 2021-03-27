#!/usr/bin/env bash

#rRNA/index database

if [ $# -eq 0 ]
  then
    echo "Supply number of cores!"
fi

CORES=$1

aDir=$(pwd)
cd references/BOWTIE_rRNA_UniVec/
bowtie-build --threads $CORES allrRNA_univec.fasta rRNA_UniVec
cd $aDir

FILES=$(ls temp/out/)
mkdir temp/mapped
mkdir temp/alignment
mkdir temp/without_rRNA
mkdir temp/log
touch temp/log/bowtie.log

for file in $FILES
do
	echo $file >> temp/log/bowtie.log
	bowtie -S -n 1 --threads $CORES --phred33-quals references/BOWTIE_rRNA_UniVec/rRNA_UniVec \
	temp/out/$file --un temp/without_rRNA/$file --al temp/mapped/$file > temp/alignment/$file 2>> temp/log/bowtie.log

done

python3 scripts/summary_from_bowtie.py temp/log/bowtie.log

RESULT=$?
if [ $RESULT -eq 0 ]; then
  cp temp/log/summary_bowtie.tab log/C_rRNA.tab
  rm -r -f temp/log 
  rm -r -f temp/alignment 
  rm -r -f temp/mapped
  rm -r -f temp/out 
else
  echo "Job failed."
fi

