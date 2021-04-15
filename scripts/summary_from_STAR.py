#!/usr/bin/env python

#takes one command line argumnet - directory, where log files (*Final.log.out) from STAR are stored and makes summary

import os
import sys

directory= sys.argv[1]

output = "Sample\tInput reads\tUnique\tUnique%\tMultimapping\tMultimapping%\tToo many mismatches%\tToo many loci%\tToo short%\tunmapped other%\tMismatch rate\tMapped_Length\n"

for file in os.listdir(directory):
	filepath = directory + file
	if "UMI" in filepath:
		output += "UMI" + file.split("_")[0].replace("UMI","") + "_" + file.split("_")[1] + "_" + file.split("_")[2].split(".")[0] + "\t"
		print(output)
	else:
		output += file.split("_")[0] + "_" + file.split("_")[1] + "_" + file.split("_")[2].split(".")[0] + "\t" 
		print(output)
	with open (filepath, mode = "r", encoding = "utf-8") as f:
		for line in f:
			if "Number of input reads" in line:
				input_reads = line.strip().split("|")[1].split("\t")[1]
			if "Uniquely mapped reads number" in line:
				unique = line.strip().split("|")[1].split("\t")[1]
			if "Number of reads mapped to multiple loci" in line:
				multi = line.strip().split("|")[1].split("\t")[1]
			if "too many mismatches" in line:
				mismatches = line.strip().split("|")[1].split("\t")[1]
			if "% of reads mapped to too many loci" in line:
				too_many_loci = line.strip().split("|")[1].split("\t")[1]
			if "too short" in line:
				too_short = line.strip().split("|")[1].split("\t")[1]
			if "of reads unmapped: other" in line:
				other_unmapped = line.strip().split("|")[1].split("\t")[1]
			if "Mismatch rate per base" in line:
				mismatch_rate = line.strip().split("|")[1].split("\t")[1]
			if "Average mapped length" in line:
				mapped_length = line.strip().split("|")[1].split("\t")[1]

	output += input_reads + "\t"
	output += unique + "\t"
	output += str(round((int(unique) / int(input_reads))*100,2))  + "\t"
	output += multi + "\t"
	output += str(round((int(multi) / int(input_reads))*100,2)) + "\t"
	output += mismatches.replace("%","") + "\t"
	output += too_many_loci.replace("%","") + "\t"
	output += too_short.replace("%","") + "\t"
	output += other_unmapped.replace("%","") + "\t"
	output += mismatch_rate.replace("%","") + "\t"
	output += mapped_length.replace("%","") + "\n"

filepath = directory + "1_Summary_log.tab"
with open (filepath, mode = "w", encoding = "utf-8") as out:
	out.write(output)
