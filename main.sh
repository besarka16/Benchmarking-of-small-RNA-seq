#!/usr/bin/env bash

#first argument specifies number of cores available for usage (8)

#Trimming
./scripts/1_trimming.sh /home/sarka/BTU_39/\!\!Users/SARKA/miRNA_kits/GEO_submission_200428/Samples/ 8

#Filtering
./scripts/2_filtering.sh 8

#remove rRNA
./scripts/3_remove_rRNA.sh 8

#remove spikes
./scripts/4_mapping_to_spikes.sh 8

#map to miRXplore
./scripts/5_mapping_to_miRXplore.sh 8 

#map to genome
./scripts/6_mapping_to_genome.sh 8 

#map to mature
./scripts/7_mapping_to_mature.sh 8

#map to hairpin
./scripts/8_mapping_to_hairpin.sh 8 

#map to isomiRs
./scripts/9_prepare_for_mapping_to_isomiRs.sh

./scripts/10_mapping_to_isomiRs.sh 8

#map to tRNA
./scripts/11_mapping_to_tRNA.sh 8 

#map to piRNA
./scripts/12_mapping_to_piRNA.sh 8

#map to ncRNA
./scripts/13_mapping_to_ncRNA_ensembl.sh 8

#deduplication of UMI
./scripts/14_deduplicate_UMI.sh

#prepare for R
./scripts/15_prepare_for_R.sh