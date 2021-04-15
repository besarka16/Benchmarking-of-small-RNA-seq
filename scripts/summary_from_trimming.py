#!/usr/bin/env python

#takes one command line argumnet - directory where all log files from cutadapt are stored

import os
import sys

directory=sys.argv[1]

output = "Sample\tInput reads\tReads with adaptor\tReads with adaptor %\tReads too short after trimming\tReads too short after trimming %\n"
for file in os.listdir(directory):
    filepath = directory + file
    if "UMIQ" in file:
        output += "UMI" + file.split("_")[0].replace("UMI","") + "_" + file.split("_")[1] + "_" + file.split("_")[2].split(".")[0] + "\t"
        with open (filepath, mode = "r", encoding = "utf-8") as f:
            for line in f:
                if "INFO Input Reads" in line:
                    input_reads = line.strip().split()[5]
                elif "INFO Reads output" in line:
                    adaptors = line.strip().split()[5]

                    input_reads = input_reads.replace(",","",3)
                    adaptors = adaptors.replace(",","",3)
                    adaptorsp = round((int(adaptors)/int(input_reads))*100,2)

                    output += input_reads + "\t"
                    output += str(adaptors) + "\t"
                    output += str(adaptorsp) + "\n"
    else:
        with open (filepath, mode = "r", encoding = "utf-8") as f:
            if "UMI" in file:
                output += "UMI" + file.split("_")[0].replace("UMI","") + "_" + file.split("_")[1] + "_" + file.split("_")[2].split(".")[0] + "\t"
            else:
                output += file.split("_")[0] + "_" + file.split("_")[1] + "_" + file.split("_")[2].split(".")[0] + \t"
            for line in f:
                if "Total reads processed" in line:
                    input_reads = line.strip().split()[3]
                elif "Reads with adapters" in line:
                    adaptors = line.strip().split()[3]
                elif "Reads written" in line:
                    written = line.strip().split()[4]
                    input_reads = input_reads.replace(",","",2)
                    adaptors = adaptors.replace(",","",3)
                    written = written.replace(",","",3)
                    short = int(input_reads) - int(written)
                    adaptorsp = round((int(adaptors)/int(input_reads))*100,2)
                    shortp = round((int(short)/int(input_reads))*100,2)
                    output += input_reads + "\t"
                    output += str(adaptors) + "\t"
                    output += str(adaptorsp) + "\t"
                    output += str(short) + "\t"
                    output += str(shortp) + "\n"


outfile = directory + "1_summary_log.tab"

with open(outfile, mode="w", encoding = "utf-8") as out:
    out.write(output)
