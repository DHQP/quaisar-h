#!/bin/sh -l

#$ -o run_MLST.out
#$ -e run_MLST.err
#$ -N run_MLST
#$ -cwd
#$ -q short.q

#
#  Description: Script to find mlst profile using tseeman's mlst script and pubmlst
# 	The most similar match is identified and provided for confirmation
#
# Usage: ./run_MLST.sh -n sample_name -p run_ID [-f DB_to force_comparison_to] [-c path_to_config_file]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/MLST/
#
# Modules required: mlst/2.16, perl/5.16.1-MT
#
# v1.0.1 (09/08/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml mlst/2.16 perl/5.16.1-MT

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_MLST.sh -n sample_name -p run_ID [-f DB_to force_comparison_to] [-c path_to_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?c:p:n:f:" option; do
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
		f)
			echo "Option -f triggered, argument = ${OPTARG}"
			forced_scheme=${OPTARG};;
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

# Checks to see if a mlst folder already exists and creates it if not
if [ ! -d "${OUTDATADIR}/MLST" ]; then
	echo "Creating ${OUTDATADIR}/MLST"
	mkdir -p "${OUTDATADIR}/MLST"
fi

# Call mlst against the given mlst DB
#. "${shareScript}/module_changers/perl_5221_to_5161mt.sh"
if [[ ! -z "${forced_scheme}" ]]; then
	#echo "(${forced_scheme})-In forced"
	mlst --scheme "${forced_scheme}" "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta" > "${OUTDATADIR}/MLST/${sample_name}_${forced_scheme}.mlst"
# Call mlst using built autoidentifier
else
	mlst "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta" > "${OUTDATADIR}/MLST/${sample_name}.mlst"
fi
#. "${shareScript}/module_changers/perl_5161mt_to_5221.sh"
#Script exited gracefully (unless something else inside failed)

ml -mlst/2.16 -perl/5.16.1-MT

exit 0
