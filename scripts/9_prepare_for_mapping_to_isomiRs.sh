#!/usr/bin/env bash

working_dir=$(pwd)

#preparation of inputs
cp -r isomiRROR-master isomiRROR-plasma
cp -r isomiRROR-master isomiRROR-miRXplore

cp isomiRROR-plasma/templates/samples.csv isomiRROR-plasma/
cp isomiRROR-miRXplore/templates/samples.csv isomiRROR-miRXplore/

cp references/hairpin_miRXplore.fa isomiRROR-miRXplore/refs/
cp references/miRXplore.fa isomiRROR-miRXplore/refs/

cp references/hsa_hairpin.fa isomiRROR-plasma/refs/
cp references/hsa_mature.fa isomiRROR-plasma/refs/

echo "#References
reference_folder: refs
mature_db_file: miRXplore
hairpin_ref_file: hairpin_miRXplore
#Input variables
min_length: 16
max_length: 26" > isomiRROR-miRXplore/config.yaml

echo "#References
reference_folder: refs
mature_db_file: hsa_mature
hairpin_ref_file: hsa_hairpin
#Input variables
min_length: 16
max_length: 28" > isomiRROR-plasma/config.yaml

cp temp/mapped_miRXplore/*miRXplore*.bam isomiRROR-miRXplore/data
cp temp/mapped_mature/*plasma*.bam isomiRROR-plasma/data

FILES=$(ls isomiRROR-miRXplore/data/ )
for file in $FILES
do
samtools fastq isomiRROR-miRXplore/data/$file > isomiRROR-miRXplore/data/$file.fastq
done

rm isomiRROR-miRXplore/data/*.bam

FILES=$(ls isomiRROR-plasma/data/ )
for file in $FILES
do
samtools fastq isomiRROR-plasma/data/$file > isomiRROR-plasma/data/$file.fastq
done

rm isomiRROR-plasma/data/*.bam

cd isomiRROR-miRXplore/data/

for file in *.fastq
do
  mv "$file" "${file/.fastq.gzAligned.sortedByCoord.out.bam.fastq/.fastq}"
done

cd $working_dir

cd isomiRROR-plasma/data/

for file in *.fastq
do
  mv "$file" "${file/.fastqAligned.sortedByCoord.out.bam.fastq/.fastq}"
done

cd $working_dir

echo "samples" > isomiRROR-miRXplore/samples.csv
ls isomiRROR-miRXplore/data/ | sed 's/.fastq//g' >> isomiRROR-miRXplore/samples.csv

echo "samples" > isomiRROR-plasma/samples.csv
ls isomiRROR-plasma/data/ | sed 's/.fastq//g' >> isomiRROR-plasma/samples.csv

