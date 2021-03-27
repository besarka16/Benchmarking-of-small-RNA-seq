#!/usr/bin/env bash

if [ $# -eq 0 ]
  then
    echo "Supply number of cores!"
fi

CORES=$1

#indexing genomes (genomeSAindexNbases as min(14,log2(genomeLength)/2-1))
#genomeChrBinNbits as min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)])

#mature miRNA
mkdir references/STAR_mature/

STAR --runThreadN $CORES --runMode genomeGenerate --genomeSAindexNbases 7 --genomeChrBinNbits 6 --limitGenomeGenerateRAM 60000000000 --genomeDir references/STAR_mature/ \
--genomeFastaFiles references/hsa_mature.fa

mkdir temp/mapped_mature

FILES=$(ls temp/mapped_genome/fastq)

for file in $FILES; do 
    mv -- "temp/mapped_genome/fastq/$file" "temp/mapped_genome/fastq/${file%.fastq.gz.Aligned.sortedByCoord.out.bam.fastq}.fastq"
done

FILES=$(ls temp/mapped_genome/fastq)

for file in $FILES
do
	STAR --runThreadN $CORES --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 15 --outFilterScoreMinOverLread 0 \
	--outFilterMatchNminOverLread 0 --alignIntronMax 1 --alignEndsType EndToEnd \
	--genomeDir references/STAR_mature/ \
	--readFilesIn temp/mapped_genome/fastq/$file \
	--outFileNamePrefix temp/mapped_mature/$file \
	--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 32000000000 --outReadsUnmapped Fastx
done

working_dir=$(pwd)

cd temp/mapped_mature/
mkdir 1_log_final
mv *Log.final.out 1_log_final/
mkdir 2_log_out
mv *Log.out 2_log_out/
mkdir 3_progress_log
mv *progress.out 3_progress_log/
mkdir 4_SJ_out
mv *SJ.out.tab 4_SJ_out/
mkdir 5_STARtmp
mv *fastq_STARtmp 5_STARtmp/
mkdir 6_unmapped
mv *mate1 6_unmapped/

cd $working_dir

python3 scripts/summary_from_STAR.py temp/mapped_mature/1_log_final/
cp temp/mapped_mature/1_log_final/1_Summary_log.tab log/G_mature.tab

