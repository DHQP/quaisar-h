#!/bin/sh -l

#$ -o get_Assemblies_from_Instruments.out
#$ -e get_Assemblies_from_Instruments.err
#$ -N get_Assemblies_from_Instruments
#$ -cwd
#$ -q short.q

# Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Will find all assembly files (.fna or .fasta) within the given folder
#
# Usage: ./get_Assemblies_from_folder.sh full_path_Output_directory folder_with_Assemblies
#
# Output location: default_config.sh_output_location/run_ID
#
# Modules required: None
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty project name supplied to $0, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./get_Assemblies_from_folder.sh  full_path_Output_directory location_of_Assemblies"
	echo "Output by default is downloaded to ${processed}/run_ID and copied to ${processed}/run_ID/sample_name/Assembly"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "Empty source supplied to $0, exiting"
	exit 1
fi


# Sets folder to where files will be downloaded to
OUTDATADIR="${1}"
if [ ! -d "${OUTDATADIR}" ]; then
	echo "Creating $OUTDATADIR"
	mkdir -p "${OUTDATADIR}"
fi

if [ -f "${OUTDATADIR}/${1}_list.txt" ]; then
	rm "${OUTDATADIR}/${1}_list.txt"
fi

# Goes through given folder
echo "${2}"
for file in ${2}/*
do
	# Check if file is a recognized assembly format externsion
	if [[ "${file}" = *.fasta ]] || [[ "${file}" = *.fna ]]; then
		filename=$(basename -- "$file")
		extension="${filename##*.}"
		sample="${filename%.*}"

		mkdir -p ${OUTDATADIR}/${sample}/Assembly
		cp ${file} ${OUTDATADIR}/${sample}/Assembly/${sample}.fasta
		cp ${OUTDATADIR}/${sample}/Assembly/${sample}.fasta ${OUTDATADIR}/${sample}/Assembly/scaffolds.fasta
		echo -e "${1}/${sample}" >> "${OUTDATADIR}/${1}_list.txt"
		#python3 ${shareScript}/removeShortContigs.py -i ${OUTDATADIR}/${sample}/Assembly/${sample}.fasta -t 500 -s normal_SPAdes
	else
		echo "${file} is not an fna or fasta file, not acting on it"
	fi
	# Invert list so that the important isolates (for us at least) get run first
	if [[ -f "${OUTDATADIR}/${1}_list.txt" ]]; then
		sort -k2,2 -t'/' -r "${OUTDATADIR}/${1}_list.txt" -o "${OUTDATADIR}/${1}_list.txt"
	fi
done

#Script exited gracefully (unless something else inside failed)
exit 0
