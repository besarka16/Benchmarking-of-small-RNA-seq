#!/usr/bin/env bash

if [ $# -eq 0 ]
  then
    echo "Supply number of cores!"
fi

CORES=$1

#indexing genomes (genomeSAindexNbases as min(14,log2(genomeLength)/2-1))
#genomeChrBinNbits as min(18,log2[max(GenomeLength/NumberOfReferences,ReadLength)])


mkdir references/STAR_genome/


STAR --runThreadN $CORES --runMode genomeGenerate --genomeDir references/STAR_genome/ \
 --genomeFastaFiles \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.10.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.11.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.12.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.13.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.14.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.15.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.16.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.17.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.18.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.19.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.1.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.20.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.21.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.22.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.2.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.3.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.4.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.5.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.6.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.7.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.8.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.9.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.MT.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.nonchromosomal.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.X.fa \
references/HSA_Genome_GRCh38/Homo_sapiens.GRCh38.dna_rm.chromosome.Y.fa 
 
mkdir temp/mapped_genome

FILES=$(ls temp/without_spikes)

for file in $FILES
do
	STAR --runThreadN $CORES --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 15 --outFilterScoreMinOverLread 0 \
	--outFilterMatchNminOverLread 0 --alignIntronMax 1 --alignEndsType EndToEnd --outFilterMultimapNmax 40 \
	--genomeDir references/STAR_genome/ \
	--readFilesIn temp/without_spikes/$file \
	--outFileNamePrefix temp/mapped_genome/$file. \
	--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 32000000000 --outReadsUnmapped Fastx
done

working_dir=$(pwd)

cd temp/mapped_genome
mkdir 1_log_final
mv *Log.final.out 1_log_final/
mkdir 2_log_out
mv *Log.out 2_log_out/
mkdir 3_progress_log
mv *progress.out 3_progress_log/
mkdir 4_SJ_out
mv *SJ.out.tab 4_SJ_out/
mkdir 5_STARtmp
mv *gz._STARtmp 5_STARtmp/
mkdir 6_unmapped
mv *mate1 6_unmapped/

cd $working_dir

mkdir  temp/mapped_genome/fastq/
FILES=$(ls temp/mapped_genome/ | grep 'Aligned.sortedByCoord.out.bam')
for file in $FILES
do	
	samtools fastq temp/mapped_genome/$file > temp/mapped_genome/fastq/$file.fastq
done

python3 scripts/summary_from_STAR.py temp/mapped_genome/1_log_final/
cp temp/mapped_genome/1_log_final/1_Summary_log.tab log/F_genome.tab
