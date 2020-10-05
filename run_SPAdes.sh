#!/bin/sh -l

#$ -o run_SPAdes.out
#$ -e run_SPAdes.err
#$ -N run_SPAdes
#$ -cwd
#$ -q short.q

#
# Description: Runs SPAdes on sample to align reads into best possible assembly
#
# Usage: ./run_spades.sh -n sample_name -t normal/continue -p run_ID	[-c path_config_file]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/Assembly
#
# Modules required: SPAdes/3.10.1
#
# v1.0.2 (08/18/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml SPAdes/3.13.0

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_spades.sh -n sample_name -t normal/continue -p run_ID	[-c path_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?c:p:n:t:" option; do
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
		t)
			echo "Option -t triggered, argument = ${OPTARG}"
			type=${OPTARG};;
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
elif [ -z "${type}" ]; then
	echo "Empty analysis type given. Exiting"
	exit 1
fi

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${processed}/${project}/${sample_name}"


#Calls spades depending on if it is supposed to look for plasmids or not, all other arguments are the same and pulled from config.sh
if [ "${type}" = "normal" ]; then
	spades.py --careful --only-assembler --pe1-1 "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" --pe1-2 "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" --pe1-s "${OUTDATADIR}/trimmed/${sample_name}.single.fq" -o "${OUTDATADIR}/Assembly" --phred-offset "${phred}" -t "${procs}"
elif [ "${type}" = "continue" ]; then
 	spades.py -o "${OUTDATADIR}/Assembly"
else
	echo "Unknown type requested...not running SPAdes"
fi

ml -SPAdes/3.13.0

#Script exited gracefully (unless something else inside failed)
exit 0
