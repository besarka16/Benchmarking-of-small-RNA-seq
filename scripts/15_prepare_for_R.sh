#!/usr/bin/env bash

mkdir R/src/mapped_mature/
mv temp/mapped_mature/*.bam R/src/mapped_mature

mkdir R/src/mapped_miRXplore/
mv temp/mapped_miRXplore/*.bam R/src/mapped_miRXplore

mkdir R/src/mapped_ncRNA/
mv temp/mapped_ncRNA/*.bam R/src/mapped_ncRNA

rm src/mapped_mature/UMI*
rm src/mapped_miRXplore/UMI*

mv mapped_mature/deduplicated_bams/* mapped_mature/
mv mapped_miRXplore/deduplicated_bams/* mapped_miRXplore/

mkdir R/src/log_mapping
mv log/* R/src/log_mapping
mv isomiRROR-miRXplore/isomir_readcount.txt R/src/log_mapping/miRXplore_isomiRs.txt
mv isomiRROR-plasma/isomir_readcount.txt R/src/log_mapping/plasma_isomiRs.txt
