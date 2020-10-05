#!/bin/sh -l

#$ -o get_Assemblies_from_Instruments.out
#$ -e get_Assemblies_from_Instruments.err
#$ -N get_Assemblies_from_Instruments
#$ -cwd
#$ -q short.q

#
# Description: Will find all assembly files (.fna or .fasta) within the given folder
#
# Usage: ./get_Assemblies_from_folder.sh -p run_ID -i folder_with_Assemblies -o output_directory [-c config_file_path]
#
# Output location: output_directory/run_ID
#
# Modules required: None
#
# v1.0.1 (08/18/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./get_Assemblies_from_folder.sh -i input_folder -o output folder -p project_ID [-c path_to_config_file]"
	echo "Output is saved to ${processed}/run_ID/ where processed is retrieved from config file, either default or imported"
}

# Parse command line options
options_found=0
while getopts ":h?i:p:o:c:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		i)
			echo "Option -i triggered, argument = ${OPTARG}"
			input_folder=${OPTARG};;
		p)
			echo "Option -p triggered, argument = ${OPTARG}"
			project=${OPTARG};;
		o)
			echo "Option -o triggered, argument = ${OPTARG}"
			output_folder=${OPTARG};;
		c)
			echo "Option -c triggered, argument = ${OPTARG}"
			config=${OPTARG};;
		:)
			echo "Option -${OPTARG} requires as argument";;
		h)
			show_help
			exit 0
			;;
	esac
done

# Show help info for when no options are given
if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
fi

# Checks for proper argumentation
if [[ -z "${input_folder}" ]] || [[ ! -d "${input_folder}" ]]; then
	echo "Empty/non-existent input folder (${input_folder}) supplied, exiting"
	exit 1
elif [[ -z "${output_folder}" ]]; then
	echo "Empty/non-existent output folder (${output_folder}) supplied, exiting"
	exit 1
elif [[ -z "${project}" ]]; then
	echo "Empty project ID supplied, can not place files correctly, must exit"
	exit 1
fi

if [[ -f "${config}" ]]; then
	echo "Loading special config file - ${config}"
	. "${config}"
else
	echo "Loading default config file"
	if [[ ! -f "./config.sh" ]]; then
		cp ./config_template.sh ./config.sh
	fi
	. ./config.sh
	cwd=$(pwd)
	config="${cwd}/config.sh"
fi

# Sets folder to where files will be downloaded to
OUTDATADIR="${output}/${project}"
if [ ! -d "${OUTDATADIR}" ]; then
	echo "Creating $OUTDATADIR"
	mkdir -p "${OUTDATADIR}"
fi

if [ -f "${OUTDATADIR}/${project}_list.txt" ]; then
	rm "${OUTDATADIR}/${project}_list.txt"
fi

# Goes through given folder
echo "${input_folder}"
for file in ${input_folder}/*
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
