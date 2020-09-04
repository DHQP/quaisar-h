#!/bin/bash -l

#$ -o snv_OA.out
#$ -e snv_OA.err
#$ -N snv_OA
#$ -cwd
#$ -q short.q

#
# Description: Script that sorts a list of samples into similar MLST groups and runs snvPhyl on each group
#
# Usage ./SNVPhyl_OA.sh   -l list_file  -a Outbreak_Name -o Output_folder
#
# Output location: Parameter
#
# Modules required: None
#
# v1.0.2 (09/04/2020)
#
# Created by Rich Stanton (njr5@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./SNVPhyl_OA.sh   -l list_file  -a Outbreak_Name -o Output_folder"
	echo "Output is saved to output_directory/tree_output_namepath_to_folder"
}

options_found=0
while getopts ":h?l:a:o:" option; do
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
		a)
			echo "Option -a triggered, argument = ${OPTARG}"
			analysisName=${OPTARG};;
		o)
			echo "Option -o triggered, argument = ${OPTARG}"
			outFolder=${OPTARG};;
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

if [[ ! "${list}" ]]; then
	echo "List (${list}) does not exist"
	echo "EXITING..."
	exit 1
fi

if [[ ! -d ${outFolder}/${analysisName} ]]; then
	mkdir ${outFolder}/${analysisName}
fi

ml Python3/3.5.2

python MLST_compare.py -i ${list} -o ${outFolder}/${analysisName}/${analysisName}

for k in ${outFolder}/${analysisName}/*.samples
do
	sample=$(basename $k)
	echo "qSNVPhyl.sh -l $k -o ${outFolder}/${analysisName} -n ${sample:0: -8}"
	qsub qSNVPhyl.sh -l $k -o ${outFolder}/${analysisName} -n ${sample:0: -8}
done

short_list=$(basename ${list})
cp ${list} ${outFolder}/${analysisName}/${short_list}.original

ml -Python3/3.5.2
