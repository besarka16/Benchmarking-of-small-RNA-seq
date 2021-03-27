#!/usr/bin/env python

#takes one command line argumnet - log file from bowtie

import os
import sys

file=sys.argv[1]

output = "Sample\tInput reads\tAligned reads\tAligned reads %\tFailed to align\tFailed to align %\n"

with open(file, mode = "r", encoding = "utf-8") as f:
	for line in f:
		if ".fastq.gz" in line:
			if "UMI" in line:
				sample = line.strip().split("_")[0].replace("UMI","") + "_" + line.strip().split("_")[1] + "_UMI"
			else:
				sample = line.strip().split("_")[0] + "_" + line.strip().split("_")[1]
		elif "processed" in line:
			input_reads = line.strip().split(" ")[3]
		elif "reported alignment:" in line:
			aligned = line.strip().split(" ")[8]
			alignedp = round((int(aligned)/int(input_reads))*100,2)
			failed = int(input_reads) - int(aligned)
			failedp = round((int(failed)/int(input_reads))*100,2)
			
			output += sample + "\t"
			output += input_reads + "\t"
			output += aligned + "\t"
			output += str(alignedp) + "\t"
			output += str(failed) + "\t"
			output += str(failedp) + "\n"

with open("temp/log/summary_bowtie.tab", mode="w", encoding = "utf-8") as out:
	out.write(output)