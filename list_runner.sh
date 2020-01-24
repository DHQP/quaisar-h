#!/bin/sh -l

#$ -o 	list_runner.out
#$ -e 	list_runner.err
#$ -N 	list_runner
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Quick and dirty way to perform something on a list of samples in the standard format (project/sample_name)
#
# Usage: ./list_runner.sh path_to_list
#
# Output location: Varies on contents
#
# Modules required: Varies on contents
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./act_by_list_template.sh path_for_list_file"
	exit 0
#elif [[ ! -f "${1}" ]]; then
#	echo "List file (${1}) does not exist, can not proceed"
#	exit 0
fi

sample_name="${1}"
project="${processed}/${2}"

# Loop through and act on each sample name in the passed/provided list
#while IFS= read -r var; do
#	sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
#	project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
#	echo "Found ${project}/${sample_name} in the list"
	if [[ -f "${project}/${sample_name}/MLST/${sample_name}_Oxford.mlst" ]]; then
		info=$(head -n 1 "${project}/${sample_name}/MLST/${sample_name}_Oxford.mlst")
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
	fi

	output="${assembly}"
	for entry in ${allEntries[@]}; do
		output="${output}	${entry//\//|}"
	done

	echo "${output}"
#done < "${1}"

# Send a completion email to whoever ran the script
echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
#Script exited gracefully (unless something else inside failed)
#printf "%s %s" "Act_by_list_template.sh has completed check of snd MLSTs" "${global_end_time}" | mail -s "act_by_list complete" nvx4@cdc.gov
exit 0
