#!/usr/bin/env python3

#
# Description: Script to remove all contigs in an assembly file that are smaller than a given size. Will assume sample is located in default location
#
# Usage: python3 ./removeShortContigs.py -i input_assembly_file -t threshold_to_trim_all_contigs_smaller_than -s
#
# Output location: -s SPAdes default_config.sh_output_location/run_ID/sample_name/sample_name/SPAdes/sample_name_scaffolds_trimmed.fasta or
#	-s SPAdes default_config.sh_output_location/run_ID/sample_name/sample_name/Unicycler_assemblies/sample_name_uni_assembly/sample_name_plasmid_assembly_trimmed.fasta
#
# Modules required: None
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

import sys
import glob
import fileinput
import argparse

# Parse all arguments from command line
def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to trim contigs')
	parser.add_argument('-i', '--input', required=True, help='input assembly filename')
	parser.add_argument('-t', '--threshold', required=True, help='threshold size to trim below')
	parser.add_argument('-s', '--source', required=True, help='Source filetype (normal_SPAdes or plasFlow)')
	return parser.parse_args()

# Script that will trim fasta files of any sequences that are smaller than the threshold

def trim_assembly(input_assembly, trim_threshold, input_type):
	assembly=open(input_assembly,'r')
	trimmed_assembly=input_assembly+".TRIMMED.fasta"
	trimmed_output=open(trimmed_assembly, 'w')
	line=assembly.readline().strip()
	total_size=0
	total_no_size=0
	total_cuts=0
	while line != '':
		if line [0] == ">":
			line_sections=line.split("_")
			sequence=""
			#print (line_sections[3], "vs", trim_threshold)
			if input_type == "normal_SPAdes":
				contig_size = line_sections[-3].split('|')[0]
			elif input_type == "plasFlow":
				contig_size = line_sections[-3]
			else:
				print("Unknown input type:", input_type)
			if int(contig_size) > int(trim_threshold):
				header=line
				line=assembly.readline().strip()
				while line != '' and line[0] != ">":
					sequence=sequence+'\n'+line
					line=assembly.readline().strip()
				#print(len(sequence),"+",total_size,"=", total_size + len(sequence))
				total_size = total_size + len(sequence)
				#print ("Yes",header)
				trimmed_output.write(header)
				trimmed_output.write(sequence+'\n')
			else:
				#print ("Nope",line)
				total_cuts+=1
				total_no_size=total_no_size+int(contig_size)
				line=assembly.readline().strip()
		else:
			line=assembly.readline().strip()
	trimmed_output.close
	assembly.close
	print("size:", total_size, "cut:", total_cuts,"contigs", total_no_size, "bps")

args = parseArgs()
trim_assembly(args.input, args.threshold, args.source)
