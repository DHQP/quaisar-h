#!/bin/sh -l

#$ -o run_sum.out
#$ -e run_sum.err
#$ -N run_sum
#$ -cwd
#$ -q short.q


#
# Description: Allows script to source this file to be able to pull out latest database filenames
#
# Usage ./get_latest_DBs.sh -d path_to_database_folder
#
# Output loction: std_out
#
# Modules required: None
#
# v1.0 (05/13/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./get_latest_DBs.sh -d path_to_database_folder"
}

# Parse command line options
options_found=0
while getopts ":h?d:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		d)
			echo "Option -d triggered, argument = ${OPTARG}"
			databases=${OPTARG};;
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
if [[ -z "${databases}" ]] || [[ ! -d "${databases}" ]]; then
	echo "Empty/non-existent database folder (${databases}) supplied, exiting"
	exit 1
fi


function get_ANI_REFSEQ {
	REFSEQ="NOT_FOUND"
	REFSEQ=$(find ${databases}/ANI -maxdepth 1 -name "REFSEQ_*.msh" -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
	echo "${REFSEQ}"
}

function get_ANI_REFSEQ_Date {
	REFSEQ="NOT_FOUND"
	REFSEQ=$(find ${databases}/ANI -maxdepth 1 -name "REFSEQ_*.msh" -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
	REFSEQ_date=$(echo "${REFSEQ}" | rev | cut -d'_' -f1 | rev | cut -d'.' -f1)
	echo "${REFSEQ_date}"
}

function get_srst2 {
	ResGANNCBI_srst2="NOT_FOUND"
	ResGANNCBI_srst2=$(find ${databases}/star -maxdepth 1 -name "ResGANNCBI_*_srst2.fasta" -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
	echo "${ResGANNCBI_srst2}"
}

function get_srst2_filename {
	ResGANNCBI_srst2="NOT_FOUND"
	ResGANNCBI_srst2=$(find ${databases}/star -maxdepth 1 -name "ResGANNCBI_*_srst2.fasta" -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
	ResGANNCBI_srst2_filename=$(echo "${ResGANNCBI_srst2}" | rev | cut -d'/' -f1 | rev | cut -d'_' -f1,2)
	echo "${ResGANNCBI_srst2_filename}"
}
