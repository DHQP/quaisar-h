#!/bin/sh -l

#$ -o run_prokka.out
#$ -e run_prokka.err
#$ -N run_prokka
#$ -cwd
#$ -q short.q

#
# Description: Runs prokka gene identifier on sample to discover all identifiable genes. Also necessary for downstream busco processing
#
# Usage ./run_prokka.sh -n sample_name -p run_ID [-c path_to_config_file]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/prokka
#
# Modules required: prokka/1.12, perl/5.12.3 java/jdk1.8.0_221
#
# v1.0.1 (09/08/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml prokka/1.12 perl/5.12.3 java/jdk1.8.0_221

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_prokka.sh -n sample_name -p run_ID [-c path_to_config_file]"
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

# Checks for existence of prokka output folder and deletes and recreates it if there
if [ -d "${OUTDATADIR}/prokka" ]; then  #removes old prokka results before continuing (it will complain otherwise)
	echo "Removing old prokka results ${OUTDATADIR}/prokka"
	rm -rf "${OUTDATADIR}/prokka"
fi

### Prokka to identify genes in ###
echo "Running Prokka for gene identification"
# Run prokka
prokka --outdir "${OUTDATADIR}/prokka" "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta"

#echo "About to rename files"
for pfile in ${OUTDATADIR}/prokka/*.*; do
	fullname=$(basename "${pfile}")
	ext="${fullname##*.}"
	echo "Renaming ${pfile} to ${OUTDATADIR}/prokka/${sample_name}_PROKKA.${ext}"
	mv "${pfile}" "${OUTDATADIR}/prokka/${sample_name}_PROKKA.${ext}"
done

#Script exited gracefully (unless something else inside failed)

ml -prokka/1.12 -perl/5.12.3

exit 0
