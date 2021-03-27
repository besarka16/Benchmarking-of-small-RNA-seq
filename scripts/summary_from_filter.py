#!/usr/bin/env python

#takes one command line argumnet - 2 directories where too short and too long reads are stored


import os
import sys

def file_len(fname):
	count = 0
	with open(fname) as f:
		for count,line in enumerate(f,1):
			continue
	return count 

shortDir=sys.argv[1]
longDir=sys.argv[2]

output = "Sample\tFiltered_out\tShort_reads\tShort_reads_%\tLong_reads\tLong_reads_%\n"
for file in os.listdir(shortDir):
	filepath_short = shortDir + file
	filepath_long = longDir + file
	short = file_len(filepath_short)/4
	long = file_len(filepath_long)/4
	filteredout = short + long
	shortp = (short/filteredout)*100
	longp = (long/filteredout)*100
	if "UMI" in filepath_short:
		output += file.split("_")[0].replace("UMI","") + "_" + file.split("_")[1] + "_UMI" + "\t"
	else:
		output += file.split("_")[0] + "_" + file.split("_")[1] + "\t" 
	output += str(filteredout) + "\t" + str(short) + "\t" + str(shortp) + "\t" + str(long) + "\t" + str(longp) + "\n"

with open("log/B_filtering.tab", mode="w", encoding = "utf-8") as out:
	out.write(output)
