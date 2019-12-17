#!/usr/bin/env python

#
# Description: Script to convert our seqlog xlsx file to a command line manageable tsv file
#   *** Must be run in python2, currently. Working on converting it to python3
#
# Usage: python2 ./xlsx_converter.py input_xlsx_file sheet_name_to_convert
#
# Output location: standard out
#
# Modules required: None
#
# v1.0 (10/3/2019)
#


#from __future__ import print_function
import os,sys,csv,pandas as pd,argparse

# Parse all arguments from command line
def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to trim contigs')
	parser.add_argument('-i', '--input', required=True, help='input excel filename')
	parser.add_argument('-r', '--run', required=True, help='Run ID to match t')
	parser.add_argument('-s', '--sheet', required=True, help='sheetname')
	parser.add_argument('-o', '--output', required=True, help='Output file to export to')
	return parser.parse_args()

def do_conversion(excel_filename, sheetname_in, output_name, run_name):
	seqlog = pd.read_excel(excel_filename, sheet_name=sheetname_in)  #usecols['CDC Aliquot ID (Miseq_ID)','Output Folder Name']
	matching_isolates=[]
	for index, row in seqlog.iterrows():
		print(index,row)
		if row['Output Folder Name'] == run_name:
			matching_isolates.append(run_name+"/"+row['CDC Aliquot ID (Miseq_ID)'])
	print("Matching rows: {0}".format(len(matching_isolates)))
	for match in matching_isolates:
		print(match+"\n")

args = parseArgs()
do_conversion(args.input, args.sheet, args.output, args.run)
