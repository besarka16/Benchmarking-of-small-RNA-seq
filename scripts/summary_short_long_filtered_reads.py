#!/usr/bin/env python

import os
import sys

directory=sys.argv[1]
output = "Sample\tShort reads\tShort reads%\tLong reads\tLong reads%\n"
output_distribution = "Sample\tLength_of_read\n"

for file in os.listdir(directory):
	filepath = directory + file
	seq_bool = False
	short_counter = 0
	long_counter = 0
	good_length_counter= 0
	sample = file[0:25]
	with open(filepath, mode = "r", encoding = "utf-8") as f:
		print(filepath)
		for line in f:
			if line.startswith("@"):
				seq_bool = True
			else:
				if seq_bool == True:
					output_distribution += sample + "\t" + str(len(line.strip())) + "\n"
					if len(line.strip()) > 28:
						long_counter += 1
						seq_bool = False
					elif len(line.strip()) < 16:
						short_counter += 1
						seq_bool = False
					else:
						seq_bool = False
						good_length_counter += 1
	f.close()
	print(sample)
	print(short_counter)
	print(long_counter)
	print(good_length_counter)

	srp = round(((short_counter/(long_counter + short_counter + good_length_counter))*100),2)
	lrp = round(((long_counter/(long_counter + short_counter + good_length_counter))*100),2)
	new_line = sample + "\t" + str(short_counter) + "\t" + str(srp)+ "\t" + str(long_counter) + "\t" + str(lrp)+ "\n"
	output += new_line

with open("short_long_filtered_reads.tab", mode = "w", encoding = "utf-8") as out:
	out.write(output)

with open("read_length_distribution.tab", mode = "w", encoding = "utf-8") as out2:
	out2.write(output_distribution)
