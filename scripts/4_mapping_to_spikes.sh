#!/usr/bin/env bash

if [ $# -eq 0 ]
  then
    echo "Supply number of cores!"
fi

CORES=$1

aDir=$(pwd)
cd references/BOWTIE_Spikes
#rRNA/index database
bowtie-build ExiSEQ-NGS-QC-Spike-ins_plus_custom.fa spikes
cd $aDir

FILES=$(ls temp/without_rRNA/)
mkdir temp/mapped
mkdir temp/alignment
mkdir temp/without_spikes
mkdir temp/log
touch temp/log/bowtie.log

for file in $FILES
do
        echo $file >> temp/log/bowtie.log
	bowtie -S -n 0 --threads $CORES --phred33-quals references/BOWTIE_Spikes/spikes \
	temp/without_rRNA/$file --un temp/without_spikes/$file --al temp/mapped/$file > \
	temp/alignment/$file 2>> temp/log/bowtie.log
done
cat temp/log/bowtie.log | sed '/^[[:space:]]*$/d' >> temp/log/bowtie2.log
python3 scripts/summary_from_bowtie.py temp/log/bowtie2.log

RESULT=$?
if [ $RESULT -eq 0 ]; then
  cp temp/log/summary_bowtie.tab log/D_spikes.tab
  rm -r -f temp/log 
  rm -r -f temp/alignment 
  rm -r -f temp/mapped
  rm -r -f temp/out 
  rm -r -f temp/without_rRNA 
else
  echo "Job failed."
fi

