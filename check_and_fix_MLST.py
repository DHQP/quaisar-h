#!/usr/bin/env python3

#
# Description: Script to go through an isolates MLST output and check and reformat any instances of partial or multiple profiles found. SUB for profile or allele means that there is
# 	a novel allele or the profile has not been assigned yet, hence something needs to be submitted to pubmlst. AU (Allele Unknown) implies that an allele can not be determined
# 	no profile can be determined from the current assembly/qc_reads
#
# Usage: is python3 ./check_and_fix_MLST.py -i input_MLST_file -t filetype_of_MLST_input_file (standard or srst2)
#
# Output location: default_config.sh_output_location/run_ID/sample_name/MLST/
#
# Modules required: None
#
# v1.0.1 (1/29/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

import sys
import os
import glob
import math
import argparse
import itertools as it
from pathlib import Path

def parseArgs(args=None):
	parser = argparse.ArgumentParser(description='Script to check MLST types for duplicate alleles and implications on final typing')
	parser.add_argument('-i', '--input', required=True, help='input mlst filename')
	parser.add_argument('-t', '--filetype', required=True, help='filetype of mlst file (standard or srst2)')
	return parser.parse_args()

# main function that looks if all MLST types are defined for an outptu mlst file
def do_MLST_check(input_MLST_file, MLST_filetype):
	# Must check if input_MLST_file has more than 1 line, different versions of MLST make different outputs
	MLST_changed_file="/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/updated_MLSTs.txt"
	blanks_file="/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/blank_MLSTs.txt"
	filepath=input_MLST_file[::-1].split("/")[2:4]
	#print(filepath)
	for i in range(0, len(filepath)):
		#print(filepath[i])
		filepath[i]=filepath[i][::-1]
	filepath=filepath[::-1]
	filepath="/".join(filepath)

	bad_types=['-', 'AU', 'SUB', 'NAF']
	# Outdated now because all old versions should be fixed, but keeping for a little while longer
	change_to_SUB=["NAM","PAM","NID","NAM&PAM", "PAM&NAM", "NF"]
	change_to_AU=["AMI", "NAM&AMI", "PAM&AMI", "NAM&PAM&AMI", "PAM&NAM&AMI"]
	types=""
	schemes=[]
	MLST_file=open(input_MLST_file,'r')
	MLST_line=MLST_file.readline().strip()
	MLST_items=MLST_line.split("	")
	#print("\n".join(MLST_items))
	if MLST_filetype == "standard":
		sample=MLST_items[0]
		db_name=MLST_items[1]
		MLST_temp_type=MLST_items[2]
		allele_list=[]
		allele_names=[]
		allele_count=len(MLST_items)
		# Default list size in case it is Empty
		list_size=0
		for allele in range(3, allele_count):
			#print(MLST_items[allele])
			allele_Identifier=MLST_items[allele].split("(")[0]
			alleles=MLST_items[allele].split("(")[1].split(")")[0].split(",")
			allele_names.append(allele_Identifier)
			allele_list.append(alleles)
		MLST_file.close()
		#allele_list=[['1'], ['3'], ['189','3'], ['2'], ['2'], ['96','107'], ['3']]
	elif MLST_filetype == "srst2":
		genus=input_MLST_file.split("_")[-2]
		species=input_MLST_file.split("_")[-1].split(".")[0]
		db_name=find_DB_taxonomy(genus, species)
		allele_list=[]
		allele_names=[]
		for i in range(2, len(MLST_items)):
			if MLST_items[i] == "mismatches":
				break
			else:
				allele_names.append(MLST_items[i])
		MLST_line_two=MLST_file.readline().strip()
		MLST_items_second=MLST_line_two.split("	")
		MLST_temp_type=MLST_items_second[1]
		sample=MLST_items_second[0]
		for i in range(2, 2+len(allele_names)):
			if '*' in MLST_items_second[i] or '?' in MLST_items_second[i] or '-' in MLST_items_second[i]:
				allele_list.append(MLST_items_second[i].split(","))
				print ("Appending non-int:", MLST_items_second[i].split(","))
			else:
				allele_list.append(MLST_items_second[i].split(","))
				print ("Appending int:", MLST_items_second[i].split(","))
		MLST_file.close()
	else:
		print("Unknown MLST filetype, can not continue")
		exit()
	print("RAW OLD types:", MLST_temp_type)
	MLST_temp_type=MLST_temp_type.replace("/", ",").replace("|",",")
	if "," not in MLST_temp_type:
		mlstype=[MLST_temp_type]
		mlstype_str = str(MLST_temp_type)
	else:
		mlstype=MLST_temp_type.split(",")
		mlstype_str = "|".join(mlstype)
	for i in range(0, len(mlstype)):
		if mlstype[i] in change_to_AU:
			mlstype[i]="AU"
		if mlstype[i] in change_to_SUB or "*" in mlstype[i] or "?" in mlstype[i]:
			mlstype[i]="SUB"
		if mlstype[i] not in bad_types:
			mlstype[i] = int(mlstype[i])
	#mlstype.sort()
	mlstype=sorted(mlstype, key=lambda x: str(x))

	print("Current MLST type:", mlstype)
	list_size=len(allele_list)
	print("Allele_names:", allele_names)
	print("Alleles_found:", allele_list)
	#for allele_index in range(0,len(allele_list)):
	#	allele_list[allele_index]=allele_list[allele_index].sort()
	if list_size == 7:
		schemes = it.product(allele_list[0], allele_list[1], allele_list[2], allele_list[3], allele_list[4], allele_list[5], allele_list[6])
	elif list_size == 8:
		schemes = it.product(allele_list[0], allele_list[1], allele_list[2], allele_list[3], allele_list[4], allele_list[5], allele_list[6], allele_list[7])
	else:
		print("Unknown size "+str(list_size)+" of allele_list")
	schemes=(list(schemes))
	print("All possible schemes:")
	print(schemes)
	#print(*schemes, sep = "\n")
	#print()
	checking=False
	for profile_index in range(0, len(schemes)):
		#print(profile_index, schemes[profile_index])
		temp_scheme=[]
		for temp_allele in schemes[profile_index]:
			temp_scheme.append(temp_allele)
		schemes[profile_index]=temp_scheme
		#print(profile_index, schemes[profile_index])
	if len(schemes) == 0:
		print("No schemes found???, probably needs to be deleted, just adding to blanks list currently")
		#os.remove(args.input)
		blanks=open(blanks_file,'a+')
		blanks.write(filepath+"	has no scheme database determined...checked against wrong or unknown DB\n")
		blanks.close()
	elif len(schemes) == 1:
		if mlstype[0] not in bad_types:
	 		print("This sample is singular and defined\n")
		else:
			print("This sample is singular and UNdefined\n")
			new_types=get_type(schemes, allele_names, db_name, MLST_filetype)
			checking=True
	elif len(schemes) > 1:
		if "NAF" in mlstype:
			print("This sample had no alleles found last time, must be very poor quality or compared to the wrong database\n")
			new_types=get_type(schemes, allele_names, db_name, MLST_filetype)
			checking=True
		elif "-" not in mlstype and "AU" not in mlstype and "SUB" not in mlstype:
			if len(schemes) == len(mlstype):
				print("This sample is a multiple and defined\n")
			elif len(schemes) > len(mlstype):
				print("Not enough types to match schemes, checking")
				new_types=get_type(schemes, allele_names, db_name, MLST_filetype)
				checking=True
			elif len(schemes) < len(mlstype):
				print("Not enough schemes to match types, checking")
				new_types=get_type(schemes, allele_names, db_name, MLST_filetype)
				checking=True
		else:
			print("This sample is a multiple and something is UNdefined")
			new_types=get_type(schemes, allele_names, db_name, MLST_filetype)
			checking=True
	print("Old types:", mlstype)

	if checking:
		new_types=sorted(new_types, key=lambda x: str(x))
		print("New types:", new_types)
		if mlstype != new_types or "-" in mlstype or "AU" in mlstype or "SUB" in mlstype or "NAF" in mlstype:
			for i in range(0, len(new_types)):
				#print(new_types[i])
				#if new_types[i] == -1:
				#	print("Found a -1")
				#	new_types[i] = "-"
				new_types[i] = str(new_types[i])
			new_types.sort()
			new_types='|'.join(new_types)
			print("Updating MLST types in", input_MLST_file, "from", mlstype_str, "to", new_types)
			MLST_temp_types=new_types
			# Log any incomplete/strange types found
			if '-' in MLST_temp_types or 'SUB' in MLST_temp_types:
				if MLST_temp_types.count("-", 0, len(MLST_temp_types)) + MLST_temp_types.count("SUB", 0, len(MLST_temp_types)) + MLST_temp_types.count("AU", 0, len(MLST_temp_types)) == 1:
					problem=["Profile_undefined"]
				elif "NAF" in MLST_temp_type:
					problem="NO_Alleles_for_profile"
				else:
					problem=["Profiles_undefined"]
				for i in range(0, len(schemes)):
					for j in range(0, len(schemes[i])):
						if "-" in schemes[i][j] or "AU" in schemes[i][j] or "SUB" in schemes[i][j] or  "~" in schemes[i][j] or "?" in schemes[i][j] or "*" in schemes[i][j]:
							if problem[0] == "Profile_undefined" or problem[0] == "Profiles_undefined":
								problem[0]="Allele(s)-"+str(allele_names[j])
							else:
								if allele_names[j] not in problem and "Allele(s)-"+str(allele_names[j]) not in problem:
									problem.append(allele_names[j])
				if problem[0] == "Profile_undefined" or problem[0] == "Profiles_undefined":
					print("Must submit profile(s) for :", filepath)
				elif problem[0] == "NO_Alleles_for_profile":
					print("No alleles found for sample, maybe try checking that samples was run against proper DB?")
				else:
					if MLST_filetype == "standard":
						print("Investigate/Submit allele or maybe try srst2 to fix allele issue on:", filepath)
					else:
						print("Investigate/Submit allele to fix allele issue on:", filepath)
				blanks=open(blanks_file,'a+')
				if MLST_filetype == "standard":
					blanks.write(filepath+"	standard:"+",".join(problem)+"	"+"	".join(MLST_items[1:])+"\n")
				elif MLST_filetype == "srst2":
					blanks.write(filepath+"	srst2:"+",".join(problem)+"	"+"	".join(MLST_items_second[1:])+"\n")
				blanks.close()
			# Change original type to new type(s) depending on source filetype
			if MLST_filetype == "standard":
				MLST_items[2]=MLST_temp_types
				MLST_file=open(input_MLST_file,'w')
				MLST_file.write('	'.join(MLST_items))
				MLST_file.close()
			elif MLST_filetype == "srst2":
				MLST_items_second[1]=MLST_temp_types
				MLST_file=open(input_MLST_file,'w')
				MLST_file.write('	'.join(MLST_items)+"\n")
				MLST_file.write('	'.join(MLST_items_second))
				MLST_file.close()
			MLST_changed_file_handler=open(MLST_changed_file,'a+')
			MLST_changed_file_handler.write(filepath+"	"+db_name+"	"+mlstype_str+" to "+new_types+"\n")
			MLST_changed_file_handler.close()
		else:
			print(input_MLST_file, "is as good as it gets with type", mlstype)
	else:
		print("Sticking with already found mlstype", mlstype,"\n")

# Uses the local copy of DB file to look up actual ST type
def get_type(list_of_profiles, list_of_allele_names, DB_file, source_filetype):
	types=["Not_initialized"]
	#print(list_of_profiles,":", list_of_allele_names)
	if source_filetype == "srst2":
		if DB_file == "abaumannii":
			for i in range(0,len(list_of_allele_names)):
				list_of_allele_names[i] = "Oxf_"+list_of_allele_names[i]
		elif DB_file == "abaumannii_2":
			for i in range(0,len(list_of_allele_names)):
				list_of_allele_names[i] = "Pas_"+list_of_allele_names[i]
		else:
			print("No adjustments needed to names in list_of_allele_names")
	full_db_path="/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases/pubmlsts/"+DB_file+"/"+DB_file+".txt"
	with open(full_db_path,'r') as scheme:
		profile_size=0
		types = [-1] * len(list_of_profiles)
		#print("Size:", len(types), " &  contents:", types)
		for line in scheme:
			db_line=line.strip()
			db_items=db_line.split("	")
			if db_items[0] == "ST":
				for item in db_items:
					if item != "clonal_complex" and item != "species":
						profile_size+=1
					else:
						break
				#print(db_items[1:profile_size])
				#print(list_of_allele_names)
				if db_items[1:profile_size] != list_of_allele_names:
					print("Allele names DO NOT match...We'll have to fix this if it ever comes up")
					print("db:"+"	".join(db_items))
					print("list:"+"	".join(list_of_allele_names))
				#else:
				#	print("Allele names match, yay!")
			else:
				for index in range(0,len(types)):
					current_profile=db_items[1:profile_size]
					type(current_profile)
					type(list_of_profiles)
					#print(current_profile)
					#print(list_of_profiles[index])
					if current_profile == list_of_profiles[index]:
						print("Match-"+str(db_items[0]), current_profile)
						types[index] = int(db_items[0])
						break
	types.sort()
	for i in range(0, len(types)):
		#print("type_check#", len(types), ":", types[i])
		if types[i] == -1:
			passed="true"
			# AU - Allele(s) missing and not close to anything in current database, therefore can not do anything further
			# SUB - One or more alleles and/or profiles ned to be submitted for classification to proper DB scheme
			if len(set(list_of_profiles[i])) == 1 and list_of_profiles[i] == '-':
				types[i]="NAF"
				continue
			for locus in list_of_profiles[i]:
				#print("Test:",locus)
				if '?' in locus or '~' in locus or '*' in locus:
					passed="false"
					if types[i] != "AU":
						types[i]="SUB"
				#print(types[i])
				elif '-' in locus:
					passed="false"
					types[i]="AU"
			if passed == "true":
				types[i] = "SUB"
		#else:
		#	types[i] = (types[i])
	return types

def find_DB_taxonomy(genus, species):
	if genus == "Acinetobacter":
		if species == "baumannii#1-Oxford":
			#print("Waiting for confirmation of success for abaumannii#1")
			species="baumannii"
		elif species == "baumannii#2-Pasteur":
			#print("Waiting for confirmation of success for abaumannii#2")
			species="baumannii_2"
		else:
			print("Extra baumannii. Not a known DB variant #1/#2")
	elif genus == "Escherichia":
		if species == "coli#1-Achtman":
			#print("Waiting for confirmation of success for abaumannii#1")
			species="coli"
		elif species == "coli#2-Pasteur":
			#print("Waiting for confirmation of success for abaumannii#2")
			species="coli_2"
		else:
			print("Extra coli. Not a known DB variant #1/#2")
	elif genus == "Burkholderia" and species == "cepacia":
		return "bcc"
	db_test_species=str(genus[0:1]).lower()+species
	#print("Test_species_DB=", db_test_species)
	species_exists = os.path.exists('/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases/pubmlsts/'+db_test_species)
	genus_exists = os.path.exists('/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases/pubmlsts/'+genus.lower())
	species_path = Path('/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases/pubmlsts/'+db_test_species)
	genus_path = Path('/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases/pubmlsts/'+db_test_species)
	if species_path.exists():
		print("Found species DB:", db_test_species)
		return db_test_species
	elif genus_path.exists():
		print("Found genus DB:", genus.lower())
		return genus
	else:
		print("No database found for", genus, species)
		exit()



print("Parsing MLST file ...")
args = parseArgs()
if os.stat(args.input).st_size > 0:
	do_MLST_check(args.input, args.filetype) #, "/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases/mlst/abaumannii_Pasteur.txt") #sys.argv[3])
else:
	print(args.input,"has an empty mlst file, so it will be deleted")
	os.remove(args.input)
