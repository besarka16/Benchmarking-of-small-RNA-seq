#!/usr/bin/env bash

#Trimming

if [ $# -eq 0 ]
  then
    echo "Supply directory containing fastq samples!"
fi

if [ $# -eq 1 ]
  then
    echo "Supply number of cores!"
fi

DIR=$1
CORES=$2

mkdir ./temp
mkdir ./temp/log

#lexogen
FILES=$(ls $DIR | grep "Lexogen")

for file in $FILES
do
	cutadapt -a "TGGAATTCTCGGGTGCCAAGG" -j $CORES -m 1 -q 20 -o ./temp/$file $DIR/$file > ./temp/log/$file.log
done

#norgen
FILES=$(ls $DIR | grep "Norgen")

for file in $FILES
do
	cutadapt -a "TGGAATTCTCGGGTGCCAAGG" -j $CORES -m 1 -q 20 -o ./temp/$file $DIR/$file > ./temp/log/$file.log
done

#realseq
FILES=$(ls $DIR | grep "RealSeq")

for file in $FILES
do
	cutadapt -a "TGGAATTCTCGGGTGCCAAGG" -j $CORES -m 1 -q 20 -o ./temp/$file $DIR/$file > ./temp/log/$file.log
done

#qiaseq
FILES=$(ls $DIR | grep "QIAseq")

for file in $FILES
do
	cutadapt -a "AACTGTAGGCACCATCAAT" -j $CORES -m 1 -q 20 -o ./temp/$file $DIR/$file > ./temp/log/$file.log
	umi_tools extract --extract-method=regex --bc-pattern=".*(?P<discard_1>AACTGTAGGCACCATCAAT){s<=1}(?P<umi_1>.{12})(?P<discard_2>.*)" -L ./temp/log/UMI$file.log --stdin=$DIR/$file --stdout ./temp/UMI$file
done

#nextflex
FILES=$(ls $DIR | grep "NEXTflex")
mkdir temp/nextflex
mkdir temp/nextflex_UMI
touch list_for_UMI.txt

for file in $FILES
do
  cutadapt -O 7 -a "N{4}TGGAATTCTCGGGTGCCAAGG" -u 4 -j $CORES -m 1 -q 20 -o ./temp/$file $DIR/$file > ./temp/log/$file.log
  cutadapt -O 7 -a "TGGAATTCTCGGGTGCCAAGG" -j $CORES -m 1 -q 20 -o ./temp/nextflex/UMI$file $DIR/$file > ./temp/log/UMI$file.log
done

gunzip temp/nextflex/*
ls temp/nextflex >> list_for_UMI.txt

inDir=temp/nextflex
outDir=temp/nextflex_UMI

  #reading every 4th line starting with line 2, get first 4 characters of sequence
  awk2='NR%4==2'
  < list_for_UMI.txt parallel -P4 "cat $inDir/{} | awk '$awk2' | cut -d' ' -f2 | cut -c1-4 > $outDir/first4_{}.txt"

  #reading every 4th line starting with line 2, get last 4 characters of sequence
  < list_for_UMI.txt  parallel -P4 "cat $inDir/{} | awk '$awk2' | sed 's/^.*\(.\{4\}\)/\1/' > $outDir/last4_{}.txt"

  #pasting first UMI 4 nuc. with last UMI 4 nuc.
  < list_for_UMI.txt parallel -P4 "paste -d'\0' $outDir/first4_{}.txt $outDir/last4_{}.txt > $outDir/UMI_{}.txt"

  #quadruple UMIs
  < list_for_UMI.txt parallel -P4 "awk '{for(i=0;i<4;i++)print}' $outDir/UMI_{}.txt >$outDir/quad_UMI_{}.txt"

  # add an "_" to the front of every UMI line
  awk3='$0="_"$0'
  < list_for_UMI.txt parallel -P4 "awk '$awk3'  $outDir/quad_UMI_{}.txt > $outDir/final_UMI_{}.txt"

  # add the UMI to the fastq file identifier line
  awk4='{getline p<f} (NR%4==1){$1=$1" "$2;$2=p}1'
  < list_for_UMI.txt parallel -P4 "awk '$awk4' OFS= f=$outDir/final_UMI_{}.txt $inDir/{} > $outDir/NEXT_{}_UMItools_R1.fq"

  #remove reads from fastq with Ns in the UMI:
  < list_for_UMI.txt parallel -P4 "sed -e '/_N\|_.*N/,+3d' $outDir/NEXT_{}_UMItools_R1.fq > $outDir/NEXT_Ns_rem_{}_UMItools_R1.fq"

  #remove random 4 base pair seqs that make up the UMI from the fastq read sequence line:
  < list_for_UMI.txt parallel -P4 "cutadapt -u 4 -o $outDir/trim2_{}_forUMI_tools.fq $outDir/NEXT_Ns_rem_{}_UMItools_R1.fq"
  < list_for_UMI.txt parallel -P4 "cutadapt -m 1 -u -4 -o $outDir/trimmed_{}_forUMI_tools.fq $outDir/trim2_{}_forUMI_tools.fq"

  #remove space form the identifier of the fastq
  < list_for_UMI.txt parallel  -P4 "sed 's/ /-/' $outDir/trimmed_{}_forUMI_tools.fq > temp/{}"

rm -r -f temp/nextflex
rm -r -f temp/nextflex_UMI
rm list_for_UMI.txt

#smarter
FILES=$(ls $DIR | grep "SMARTer")

for file in $FILES
do
	cutadapt -a "AAAAAAAAAAA" -u 3 -j $CORES -m 1 -q 20 -o ./temp/$file $DIR/$file > ./temp/log/$file.log
done

#edgeseq
FILES=$(ls $DIR | grep "EdgeSeq")

for file in $FILES
do
	cutadapt -a "GATCGGAAGAGCACACGTCTGAACTC" -u 3 -j $CORES -m 1 -q 20 -o ./temp/$file $DIR/$file > ./temp/log/$file.log
done

python3 ./scripts/summary_from_trimming.py ./temp/log/

RESULT=$?
if [ $RESULT -eq 0 ]; then 
        mkdir ./log/
        cp temp/log/1_summary_log.tab ./log/A_trimming.tab
        rm -r -f temp/log
else
  	echo "Generation of log failed."
fi

