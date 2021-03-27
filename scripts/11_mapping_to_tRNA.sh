#!/usr/bin/env bash

if [ $# -eq 0 ]
  then
    echo "Supply number of cores!"
fi

CORES=$1

#indexing genomes (genomeSAindexNbases as min(14,log2(genomeLength)/2-1))
#genomeChrBinNbits as min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)])


#tRNA
mkdir references/STAR_tRNA/

STAR --runThreadN $CORES --runMode genomeGenerate --genomeSAindexNbases 5 --genomeChrBinNbits 4 --limitGenomeGenerateRAM 60000000000 --genomeDir references/STAR_tRNA/ \
--genomeFastaFiles references/hg38-tRNAs.fa

mkdir temp/mapped_tRNA/

FILES=$(ls isomiRROR-plasma/data/ | grep "unmapped.fastq")

for file in $FILES
do
	STAR --runThreadN $CORES --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0 \
	--outFilterMatchNminOverLread 0 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMultimapNmax 40 \
	--genomeDir references/STAR_tRNA/ \
	--readFilesIn isomiRROR-plasma/data/$file \
	--outFileNamePrefix temp/mapped_tRNA/$file. \
	--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 60000000000 --outReadsUnmapped Fastx
done

working_dir=$(pwd)

cd temp/mapped_tRNA/
mkdir 1_log_final
mv *Log.final.out 1_log_final/
mkdir 2_log_out
mv *Log.out 2_log_out/
mkdir 3_progress_log
mv *progress.out 3_progress_log/
mkdir 4_SJ_out
mv *SJ.out.tab 4_SJ_out/
mkdir 5_STARtmp
mv *tmp 5_STARtmp/
mkdir 6_unmapped
mv *mate1 6_unmapped/

cd $working_dir

python3 scripts/summary_from_STAR.py temp/mapped_tRNA/1_log_final/
cp temp/mapped_tRNA/1_log_final/1_Summary_log.tab log/Ch_tRNA.tab
