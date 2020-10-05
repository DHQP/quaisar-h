#!/bin/sh -l

#$ -o run_sum.out
#$ -e run_sum.err
#$ -N run_sum
#$ -cwd
#$ -q short.q

#
# Description: Creates a summary file for the run and prints out a one word status and a short description of each step being reported
#
# Usage ./run_sum.sh -p run_ID [-c path_to_config_file]
#
# Output location: default_config.sh_output_location/run_ID/
#
# Modules required: None
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_sum.sh -p run_ID [-c path_to_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?c:p:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		p)
			echo "Option -p triggered, argument = ${OPTARG}"
			project=${OPTARG};;
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

# Checks for proper argumentation
if [ -z "${project}" ]; then
	echo "Empty project name given. Exiting"
	exit 1
fi

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${processed}/${project}"


echo "Checking for ${OUTDATADIR}/${project}_list(_ordered).txt"

# Checks for existence of list files in specific order
#if [[ -z ${2} ]]; then
	if [[ -f "${OUTDATADIR}/${project}_list_ordered.txt" ]]; then
		list="${OUTDATADIR}/${project}_list_ordered.txt"
	elif [[ -f "${OUTDATADIR}/${project}_list.txt" ]]; then
		list="${OUTDATADIR}/${project}_list.txt"
	else
		echo "No list file exists, cannot do a summary, unless I add in an automagic creator later"
		exit
	fi
#	type="project"
#else
#	type="list"
#	list=${project}
#fi

# Gets todays date to show when summary was run
runsumdate=$(date "+%Y_%m_%d_at_%Hh_%Mm")
echo "Creating run summary at ${runsumdate}"
# Status of each individual sample is updated in its own folder and the run_summary file
#if [[ "${type}" = "project" ]]; then
sum_name="${project}_run_summary_at_${runsumdate}.sum"
#	echo "named as project"
#else
#	sum_name="list_summary_at_${runsumdate}.sum"
#	echo "named as list"
#fi

echo "List=${list}"

# Run validate_piperun.sh on every sample in the list and cat output into one summary file
while IFS= read -r samples || [ -n "$samples" ]; do
	file=$(echo "${samples}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	proj=$(echo "${samples}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	echo ${file}
	"${shareScript}/validate_piperun.sh" -n "${file}" -p "${proj}" -c "${config}" > "${processed}/${proj}/${file}/${file}_pipeline_stats.txt"
	#if [[ "${type}" = "project" ]]; then
		cat "${processed}/${proj}/${file}/${file}_pipeline_stats.txt" >> "${processed}/${proj}/${sum_name}"
	#else
	#	cat "${processed}/${proj}/${file}/${file}_pipeline_stats.txt" >> "${3}/${sum_name}"
	#fi
done < ${list}
