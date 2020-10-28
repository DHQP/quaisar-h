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
# v1.0.2 (05/12/2020)
#


#from __future__ import print_function
import os,sys,csv,pandas as pd,argparse

# Parse all arguments from command line
def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to trim contigs')
	parser.add_argument('-i', '--input', required=True, help='input excel filename')
	parser.add_argument('-l', '--list', required=True, help='original list to match against')
	parser.add_argument('-r', '--run', required=True, help='Run ID to match to')
	parser.add_argument('-s', '--sheet', required=True, help='sheetname')
	parser.add_argument('-o', '--output', required=True, help='Output file to export to')
	return parser.parse_args()

def create_sample_dict(list_in):
	samples=[]
	#f = open(list_in, "r")
	#line=f.readline()
	#sample=samples.append(str(line.split("/")[1]))
	#return samples

	with open(list_in) as fp:
		line = fp.readline().strip()
		cnt = 1
		while line:
			samples.append(str(line.split("/")[1]))
			line = fp.readline().strip()
			cnt += 1
	return samples


def do_conversion(excel_filename, sheetname_in, output_name, run_name, sample_list):
	sample_count=len(sample_list)
	for thing in sample_list:
		print(":"+thing+":")
	if sample_count == 0:
		print("No samples added from original list, cant compare nothing to Seqlog, exiting")
		exit()
	#for id in range(0, sample_count):
	#	print(id, sample_list[id])
	seqlog = pd.read_excel(excel_filename, sheet_name=sheetname_in)  #usecols['CDC Aliquot ID (Miseq_ID)','Output Folder Name']
	matching_isolates=[]
	missing_samples=0
	if 'OSII WGS ID (HQ)' in seqlog.columns:
		column_title='OSII WGS ID (HQ)'
	elif 'CDC Local Aliquot ID or Outbreak ID' in seqlog.columns:
		column_title='CDC Local Aliquot ID or Outbreak ID'
	elif 'CDC Aliquot ID (Miseq ID)' in seqlog.columns:
		column_title='CDC Aliquot ID (Miseq ID)'

	for index, row in seqlog.iterrows():
		#print(index,row)
		if row['Output Folder Name'] == run_name:
			print(column_title, str(row[column_title]))
			if column_title == "OSII WGS ID (HQ)":
				if str(row[column_title]) in sample_list:
					matching_isolates.append(str(run_name)+"/"+str(row[column_title]))
				# Catch if name is changed due to FAILED/META/REFERENCE
				elif str(row[column_title])+"_FAILED" in sample_list:
					matching_isolates.append(str(run_name)+"/"+str(row[column_title])+"_FAILED")
				elif str(row[column_title]) == "nan":
					print(column_title, str(row['CDC Local Aliquot ID or Outbreak ID']))
					if str(row['CDC Local Aliquot ID or Outbreak ID']) != "nan":
						if str(row['CDC Local Aliquot ID or Outbreak ID']) in sample_list:
							matching_isolates.append(str(run_name)+"/"+str(row['CDC Local Aliquot ID or Outbreak ID']))
						elif str(row['CDC Local Aliquot ID or Outbreak ID'])+"_FAILED" in sample_list:
							matching_isolates.append(str(run_name)+"/"+str(row['CDC Local Aliquot ID or Outbreak ID'])+"_FAILED")
						else:
							print("No match for non-empty", str(row['CDC Local Aliquot ID or Outbreak ID']), "in OSII WGS IS (HQ)-less CDC Local Aliquot ID or Outbreak ID")
							missing_samples+=1
				else:
					print("No match for empty", str(row['CDC Local Aliquot ID or Outbreak ID']), "in OSII WGS IS (HQ)-less CDC Local Aliquot ID or Outbreak ID")
					missing_samples+=1
			elif column_title == 'CDC Local Aliquot ID or Outbreak ID' or column_title == 'CDC Aliquot ID (Miseq ID)':
				if str(row[column_title]) in sample_list:
					matching_isolates.append(str(run_name)+"/"+str(row[column_title]))
				# Catch if name is changed due to FAILED/META/REFERENCE
				elif str(row[column_title])+"_FAILED" in sample_list:
					matching_isolates.append(str(run_name)+"/"+str(row[column_title])+"_FAILED")
				else:
					print("No match for", str(row[column_title]),"in CDC (Local) Aliquot ID [or Outbreak ID|Miseq ID]")
					missing_samples+=1
	print("Matching rows: {0}".format(len(matching_isolates)))
	print("Missing samples:", missing_samples)
	summary_out=open(output_name, 'w')
	for match in matching_isolates:
		summary_out.write(match+"\n")
	for i in range(0, missing_samples):
		summary_out.write("MISSING_SAMPLE\n")
	summary_out.close()

args = parseArgs()
do_conversion(args.input, args.sheet, args.output, args.run, create_sample_dict(args.list))
