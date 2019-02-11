#!/bin/sh -l

#$ -o run_ANI_folder_3.out
#$ -e run_ANI_folder_3.err
#$ -N run_ANI_folder_3
#$ -cwd
#$ -q all.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp config_template.sh config.sh
fi
. ./config.sh
. ${mod_changers}/pipeline_mods

#
# Script to calculate the All vs. All average nucleotide identity of a folder of fasta files
#
# Usage ./run_ANI_folder.sh path_to_folder
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to run_ANI.sh, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_ANI.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_ANI_folder.sh path_to_folder"
	exit 0
fi

python -V
echo "Running ALL vs ALL aniM on ${1} and placing results in ${1}/aniM"
#python "${shareScript}/pyani/average_nucleotide_identity.py" -i "${OUTDATADIR}/ANI/localANIDB" -o "${OUTDATADIR}/ANI/aniM" --write_excel
if [[ -d "${1}/aniM" ]]; then
	mv "${1}/aniM" "${1}/aniM_023"
fi

average_nucleotide_identity.py -i "${1}" -o "${1}/aniM"