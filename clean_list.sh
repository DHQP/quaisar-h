#!/bin/sh -l

#$ -o clean_list.out
#$ -e clean_list.err
#$ -N clean_list
#$ -cwd
#$ -q short.q

#
# Description: Script to clean any list file of extra newlines and space
#
# Usage ./clean_list.sh -l path_to_list_file
#
# Output location: same folder as path_to_list
#
# Modules required: None
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#


show_help () {
	echo "Usage is ./clean_list.sh -l path_to_list_file"
}

# Parse command line options
options_found=0
while getopts ":h?l:" option; do
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

if [[ -z "${list}" ]] || [[ ! -f ${list} ]]; then
	echo "Empty list supplied to clean_list.sh or list file does not exist, exiting"
	exit 1
fi

# Makes a backup of original list file
cp -f ${list} ${list}.original

# Puts list through dos2unix to convert any windows line returns to unix
dos2unix ${list}.original

tr -d ' \t\r\f' < ${list}.original > ${list}
ex -s +'v/\S/d' -cwq ${list}

# Shows all old lines from list
while IFS= read -r line; do
	echo "O:${line}:"
done < "${list}.original"

# Shows all new lines in list
while IFS= read -r line; do
	echo "N:${line}:"
done < "${list}"

#Script exited gracefully (unless something else inside failed)
exit 0
