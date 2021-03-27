#!/usr/bin/env bash

if [ $# -eq 0 ]
  then
    echo "Supply number of cores!"
fi

CORES=$1

#indexing genomes (genomeSAindexNbases as min(14,log2(genomeLength)/2-1))
#genomeChrBinNbits as min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)])

mkdir references/STAR_miRXplore/

STAR --runThreadN $CORES --runMode genomeGenerate --genomeSAindexNbases 6 --genomeChrBinNbits 6 --genomeDir references/STAR_miRXplore/ \
--genomeFastaFiles references/miRXplore.fa

mkdir temp/mapped_miRXplore

gunzip temp/without_spikes/*

FILES=$(ls temp/without_spikes)

for file in $FILES
do
	STAR --runThreadN $CORES --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0 \
        --outFilterMatchNminOverLread 0 --alignIntronMax 1 --alignEndsType EndToEnd \
	--genomeDir references/STAR_miRXplore/ --readFilesIn temp/without_spikes/$file --outFileNamePrefix temp/mapped_miRXplore/$file \
	--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 32000000000 --outReadsUnmapped Fastx
done

aDir=$(pwd)

cd temp/mapped_miRXplore
mkdir 1_log_final
mv *Log.final.out 1_log_final/
mkdir 2_log_out
mv *Log.out 2_log_out/
mkdir 3_progress_log
mv *progress.out 3_progress_log/
mkdir 4_SJ_out
mv *SJ.out.tab 4_SJ_out/
mkdir 5_STARtmp
mv *gz_STARtmp 5_STARtmp/
mkdir 6_unmapped
mv *mate1 6_unmapped/

cd $aDir

python3 scripts/summary_from_STAR.py temp/mapped_miRXplore/1_log_final/
cp temp/mapped_miRXplore/1_log_final/1_Summary_log.tab log/E_miRXplore.tab

