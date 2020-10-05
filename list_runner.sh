#!/bin/sh -l

#$ -o 	list_runner.out
#$ -e 	list_runner.err
#$ -N 	list_runner
#$ -cwd
#$ -q short.q

#
# Description: Quick and dirty way to perform something on a list of samples in the standard format (project/sample_name)
#
# Usage: ./list_runner.sh -l path_to_list [-c path_to_config_file]
#
# Output location: Varies on contents
#
# Modules required: Varies on contents but should be loaded by what it calls
#
# v1.0.2 (09/08/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./list_runner.sh -l path_to_list [-c path_to_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?l:c:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		l)
			echo "Option -l triggered, argument = ${OPTARG}"
			list=${OPTARG};;
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
if [[ ! -f "${list}" ]] || [[ -z "${list}" ]]; then
	echo "List empty or non-existent, exiting"
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

# Loop through and act on each sample name in the passed/provided list
while IFS= read -r var; do
	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
	if [[ -f "${processed}/${project}/${sample_name}/MLST/${sample_name}_Oxford.mlst" ]]; then
		echo echo "Found Oxford for ${project}/${sample_name} in the list"
		info=$(head -n 1 "${processed}/${project}/${sample_name}/MLST/${sample_name}_Oxford.mlst")
		assembly=$(echo "${info}" | cut -d'	' -f1)
		db=$(echo "${info}" | cut -d'	' -f2)
		st=$(echo "${info}" | cut -d'	' -f3)
		gltA=$(echo "${info}" | cut -d'	' -f4)
		gyrB=$(echo "${info}" | cut -d'	' -f5)
		gdhB=$(echo "${info}" | cut -d'	' -f6)
		recA=$(echo "${info}" | cut -d'	' -f7)
		cpn60=$(echo "${info}" | cut -d'	' -f8)
		gpi=$(echo "${info}" | cut -d'	' -f9)
		rpoD=$(echo "${info}" | cut -d'	' -f10)
		allEntries=( ${db} ${st} ${gltA} ${gyrB} ${gdhB} ${recA} ${cpn60} ${gpi} ${rpoD})

		output="${assembly}"
		for entry in ${allEntries[@]}; do
			output="${output}	${entry//\//|}"
		done

		echo "${output}" > "${processed}/${project}/${sample_name}/MLST/${sample_name}_Oxford.mlst"
	fi

done < "${list}"

# Send a completion email to whoever ran the script
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
#printf "%s %s" "Act_by_list_template.sh has completed check of snd MLSTs" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
