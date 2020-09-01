#!/bin/sh -l

#$ -o mashlist.out
#$ -e mashlist.err
#$ -N mashlist
#$ -cwd
#$ -q short.q

#
# Description: Script to create mashtree of specified isolates that were processed by Quaisar pipeline (proper folder structures)
#
# Usage: ./mashtree_of_list.sh -i absolute_path_to_list -d absolute_output_directory -o tree_output_name [-c path_to_config_file]
#
# Output location: parameter
#
# Modules required: perl/5.16.1-MT, mashtree/0.29
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml perl/5.16.1-MT mashtree/0.29

#  Function to print out help blurb
show_help () {
	echo "Usage is ./mashtree_of_list.sh -i path_to_list -d output_directory -o tree_output_name [-c path_to_config_file]"
	echo "Output is saved to output_directory/tree_output_namepath_to_folder"
}

options_found=0
while getopts ":h?i:d:o:c:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		i)
			echo "Option -i triggered, argument = ${OPTARG}"
			input=${OPTARG};;
		d)
			echo "Option -d triggered, argument = ${OPTARG}"
			outdir=${OPTARG};;
		o)
			echo "Option -o triggered, argument = ${OPTARG}"
			output_file=${OPTARG};;
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

if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
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

# create output directory if it does not exist
if [[ ! -d ${outdir} ]]; then
	mkdir -p ${outdir}
fi

# Copy over all fasta files from original locations to the output directory
while IFS= read -r line || [ "$line" ];  do
	sample_name=$(echo "${line}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${line}" | cut -d'/' -f1 | tr -d '[:space:]')
	cp ${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta ${outdir}
done < ${input}

# Call mashtree on all copied fasta
cd ${outdir}
mashtree --numcpus ${procs} *.fasta > "${outdir}/${output_file}.dnd";

ml -perl/5.16.1-MT

exit
