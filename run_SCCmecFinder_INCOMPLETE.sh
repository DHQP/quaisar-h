#!/bin/sh -l

#$ -o run_SPAdes.out
#$ -e run_SPAdes.err
#$ -N run_SPAdes
#$ -cwd
#$ -q short.q

#
# Description: Runs SCCmesFinder on isolate to find mec types in staph samples
#
# Usage: ./run_SCCmecFinder.sh -n sample_name	-p run_ID [-c path_to_config_file]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/SCCmec
#
# Modules required: None
#
# v1.0.1 (09/08/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_SCCmecFinder.sh -n sample_name	-p run_ID [-c path_to_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?c:p:n:" option; do
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
		n)
			echo "Option -n triggered, argument = ${OPTARG}"
			sample_name=${OPTARG};;
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
if [[ -z "${sample_name}" ]]; then
	echo "Empty sample name supplied to run_kraken.sh, exiting"
	exit 1
elif [ -z "${project}" ]; then
	echo "Empty project name given. Exiting"
	exit 1
fi

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${processed}/${project}/${sample_name}"



if [[ ! -d "${OUTDATADIR}/SCCmec" ]]; then
	mkdir -p "${OUTDATADIR}/SCCmec"
fi

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${OUTDATADIR}/SCCmec"

python2 "${shareScript}/SCCmecFinder_v4.py"

#Script exited gracefully (unless something else inside failed)
exit 0
