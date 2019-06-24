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
	MLST_temp_type=MLST_temp_type.replace("/", ",")
	if "," not in MLST_temp_type:
		mlstype_str = [MLST_temp_type]
		mlstype=[MLST_temp_type]
		for i in range(0, len(mlstype)):
			if mlstype[i] != '-':
				mlstype[i] = int(mlstype[i])
		mlstype.sort()
	else:
		if "," in MLST_temp_type:
			mlstype=MLST_temp_type.split(",")
		else:
			mlstype=MLST_temp_type.split("|")
	print("Current MLST type:", mlstype, "\n")
	list_size=len(allele_list)
	print("Allele_names:", allele_names)
	print("Alleles_found:", allele_list, "\n")
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
	print()
	checking=False
	for profile_index in range(0, len(schemes)):
		#print(profile_index, schemes[profile_index])
		temp_scheme=[]
		for temp_allele in schemes[profile_index]:
			temp_scheme.append(temp_allele)
		schemes[profile_index]=temp_scheme
		#print(profile_index, schemes[profile_index])
	if len(schemes) == 0:
		print("No schemes found???")
	elif len(schemes) == 1:
		if mlstype[0] != "-":
	 		print("This sample is singular and defined\n")
		else:
			print("This sample is singular and UNdefined\n")
			new_types=get_type(schemes, allele_names, db_name)
			checking=True
	elif len(schemes) > 1:
		if "-" not in mlstype:
			if len(schemes) == len(mlstype):
				print("This sample is a multiple and defined\n")
			elif len(schemes) > len(mlstype):
				print("Not enough types to match schemes")
				new_types=get_type(schemes, allele_names, db_name)
				checking=True
			elif len(schemes) < len(mlstype):
				print("Not enough schemes to match types")
				new_types=get_type(schemes, allele_names, db_name)
				checking=True
		else:
			print("This sample is a multiple and something is UNdefined")
			new_types=get_type(schemes, allele_names, db_name)
			checking=True
	print("Old types:", mlstype, "\n")
	filepath=input_MLST_file[::-1].split("/")[2:4]
	#print(filepath)
	for i in range(0, len(filepath)):
		#print(filepath[i])
		filepath[i]=filepath[i][::-1]
	filepath=filepath[::-1]
	filepath="/".join(filepath)
	if checking:
		print("New types:", new_types, "\n")
		if mlstype != new_types:
			for i in range(0, len(new_types)):
				#print(new_types[i])
				#if new_types[i] == -1:
				#	print("Found a -1")
				#	new_types[i] = "-"
				new_types[i] = str(new_types[i])
			#new_types.sort()
			new_types='/'.join(new_types)
			print("Updating MLST types in", input_MLST_file, "from", ",".join(mlstype_str), "to", new_types)
			MLST_temp_types=new_types
			# Log any incomplete/strange types found
			if '-' in MLST_temp_types or 'ND' in MLST_temp_types or 'NF' in MLST_temp_types:
				if MLST_temp_types.count("-", 0, len(MLST_temp_types)) + MLST_temp_types.count("ND", 0, len(MLST_temp_types)) + MLST_temp_types.count("NF", 0, len(MLST_temp_types)) == 1:
					problem=["Profile_undefined"]
				else:
					problem=["Profiles_undefined"]
				for i in range(0, len(schemes)):
					for j in range(0, len(schemes[i])):
						if "-" in schemes[i][j] or "~" in schemes[i][j] or "?" in schemes[i][j] or "*" in schemes[i][j]:
							if problem[0] == "Profile_undefined" or problem[0] == "Profiles_undefined":
								problem[0]="Allele(s)-"+str(allele_names[j])
							else:
								if allele_names[j] not in problem and "Allele(s)-"+str(allele_names[j]) not in problem:
									problem.append(allele_names[j])
				if problem[0] == "Profile_undefined" or problem[0] == "Profiles_undefined":
					print("Must submit profile(s) for :", filepath)
				else:
					if MLST_filetype == "standard":
						print("Investigate/Submit allele or maybe try srst2 to fix allele issue on:", filepath)
					else:
						print("Investigate/Submit allele to fix allele issue on:", filepath)
				blanks_file="/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/blank_MLSTs.txt"
				blanks=open(blanks_file,'a+')
				blanks.write(filepath+"	"+",".join(problem)+"	"+"	".join(MLST_items[1:])+"\n")
				blanks.close()
			# Change original type to new type(s) depending on source filetype
			if MLST_filetype == "standard":
				MLST_items[2]=MLST_temp_types
				MLST_file=open(input_MLST_file,'w')
				MLST_file.write('	'.join(MLST_items))
				MLST_file.close()
			elif MLST_filetype == "srst2":
				MLST_items[1]=MLST_temp_types
				MLST_file=open(input_MLST_file,'w')
				MLST_file.write('	'.join(MLST_items))
				MLST_file.write('	'.join(MLST_items_second))
				MLST_file.close()
			MLST_changed_file="/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/updated_MLSTs.txt"
			MLST_changed_file_handler=open(MLST_changed_file,'a+')
			MLST_changed_file_handler.write(filepath+"	"+db_name+"	"+",".join(mlstype_str)+" to "+new_types+"\n")
			MLST_changed_file_handler.close()
		else:
			print(input_MLST_file, "is as good as it gets with type", mlstype)

# Uses the local copy of DB file to look up actual ST type
def get_type(list_of_profiles, list_of_allele_names, DB_file):
	types=["Not_initialized"]
	full_db_path="/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases/pubmlsts/"+DB_file+"/"+DB_file+".txt"
	with open(full_db_path,'r') as f:
		profile_size=0
		types = [-1] * len(list_of_profiles)
		#print("Size:", len(types), " &  contents:", types)
		for line in f:
			db_line=line.strip()
			db_items=db_line.split("	")
			if db_items[0] == "ST":
				for item in db_items:
					if item != "clonal_complex" and item != "species":
						profile_size+=1
					else:
						break
				print(db_items[1:profile_size])
				print(list_of_allele_names)
				if db_items[1:profile_size] == list_of_allele_names:
					print("Allele names match, yay!")
				else:
					print("We'll have to fix this if it ever comes up")
					print("db: "+db_items)
					print("list:"+allele_names)
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
		if types[i] == -1:
			types[i] = "-"
		else:
			types[i] = str(types[i])
	return types

def find_DB_taxonomy(genus, species):
	if genus == "Acinetobacter":
		if species == "baumannii#1":
			print("Waiting for confirmation of success for abaumannii#1")
		elif species == "baumannii#2":
			print("Waiting for confirmation of success for abaumannii#2")
		else:
			print("Waiting for confirmation of filenames for abuamanniis")
	elif genus == "Escherichia":
		if species == "coli#1":
			print("Waiting for confirmation of success for abaumannii#1")
		elif species == "coli#2":
			print("Waiting for confirmation of success for abaumannii#2")
		else:
			print("Waiting for confirmation of filenames for ecolis")
	elif genus == "Burkholderia" and species == "cepacia":
		return "bcc"
	else:
		db_test_species=str(genus[0:1]).lower()+species
		print("Test_species_DB=", db_test_species)
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


print("Parsing MLST file ...\n")
args = parseArgs()
do_MLST_check(args.input, args.filetype) #, "/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases/mlst/abaumannii_Pasteur.txt") #sys.argv[3])