#!/usr/bin/env bash

if [ $# -eq 0 ]
  then
    echo "Supply number of cores!"
fi

CORES=$1

#indexing genomes (genomeSAindexNbases as min(14,log2(genomeLength)/2-1))
#genomeChrBinNbits as min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)])


mkdir references/STAR_ncRNA/
gunzip reference/non_coding_without_miRNA.fa.gz

STAR --runThreadN $CORES --runMode genomeGenerate --genomeSAindexNbases 14 --genomeChrBinNbits 10 --limitGenomeGenerateRAM 60000000000 --genomeDir references/STAR_ncRNA/ \
--genomeFastaFiles references/non_coding_without_miRNA.fa


mkdir temp/mapped_ncRNA/

FILES=$(ls temp/mapped_piRNA/6_unmapped/)

for file in $FILES
do
	STAR --runThreadN $CORES --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0 \
	--outFilterMatchNminOverLread 0 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMultimapNmax 40 \
	--genomeDir references/STAR_ncRNA/ \
	--readFilesIn temp/mapped_piRNA/6_unmapped/$file \
	--outFileNamePrefix temp/mapped_ncRNA/$file. \
	--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 60000000000 --outReadsUnmapped Fastx
done

working_dir=$(pwd)

cd temp/mapped_ncRNA/
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

python3 scripts/summary_from_STAR.py temp/mapped_ncRNA/1_log_final/
cp temp/mapped_ncRNA/1_log_final/1_Summary_log.tab log/J_ncRNA.tab

