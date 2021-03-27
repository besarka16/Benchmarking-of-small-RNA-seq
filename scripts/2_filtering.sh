#!/usr/bin/env bash

if [ $# -eq 0 ]
  then
    echo "Supply number of cores!"
fi

CORES=$1
FILES=$(ls temp)
mkdir temp/long
mkdir temp/short
mkdir temp/log
mkdir temp/out

for file in $FILES
do
 cutadapt -m 16 -M 28 -j $CORES --too-short-output temp/short/$file --too-long-output temp/long/$file -o temp/out/$file temp/$file > temp/log/$file.log
done

gunzip temp/long/*
gunzip temp/short/*
python3 scripts/summary_from_filter.py temp/short/ temp/long/
RESULT=$?
if [ $RESULT -eq 0 ]; then 
     rm -r -f temp/log
     rm -r -f temp/long
     rm -r -f temp/short
     rm temp/*fastq*
else
  	echo "Generation of log failed."
fi

FILES=$(ls temp/out | grep "UMINEXTflex")
for file in $FILES
do
	gzip temp/out/$file
done

