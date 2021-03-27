#!/usr/bin/env bash

if [ $# -eq 0 ]
  then
    echo "Supply number of cores!"
fi

CORES=$1

#MirXplore


FILES=$(ls temp/mapped_miRXplore/ | grep 'UMI' | grep 'miRXplore')
for file in $FILES
do	

	samtools index temp/mapped_miRXplore/$file
done

mkdir temp/mapped_miRXplore/deduplicated_bams/
mkdir temp/mapped_miRXplore/deduplication_log/

for file in $FILES
do
	umi_tools dedup -I temp/mapped_miRXplore/$file \
	-S temp/mapped_miRXplore/deduplicated_bams/$file -L temp/mapped_miRXplore/deduplication_log/$file.log 
done

#miRNA - mature

FILES=$(ls temp/mapped_mature/ | grep 'UMI'| grep 'plasma')
for file in $FILES
do	

	samtools index temp/mapped_mature/$file
done

mkdir temp/mapped_mature/deduplicated_bams/
mkdir temp/mapped_mature/deduplication_log/

for file in $FILES
do
	umi_tools dedup -I temp/mapped_mature/$file \
	-S temp/mapped_mature/deduplicated_bams/$file -L temp/mapped_mature/deduplication_log/$file.log 
done
