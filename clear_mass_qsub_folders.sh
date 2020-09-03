#!/bin/sh -l

#$ -o clear_mqs.out
#$ -e clear_mqs.err
#$ -N clear_mqs
#$ -cwd
#$ -q short.q

#
# Description: Script to clean up script and mass qsub folders of proscripts made while mass submitting many parallele jobs
#
# Usage ./clear_mass_qsub_folders -l path_to_folder_to_clean -t type_of_cleaning (1=home_Script_folder=errs/outs 2=mass_qsub_folder=errs/outs/sh's...be VERY careful with #2)
#
# Output location: No output created
#
# Modules required: None
#
# v1.0.2 (08/18/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./clear_mass_qsub_folders -l path_to_folder_to_clean -t type_of_cleaning (1=home_script_folder=errs/outs 2=mass_qsub_folder=errs/outs/sh's...be VERY careful with #2)"
}

# Parse command line options
options_found=0
while getopts ":h?l:t:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		l)
			echo "Option -l triggered, argument = ${OPTARG}"
			location=${OPTARG};;
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

# Show help info for when no options are given
if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
fi

# Number regex to test max concurrent submission parametr
number='^[1-3]+$'

# Checks for proper argumentation
if [[ -z "${location}" ]] || [[ ! -d ${location} ]]; then
	echo "Empty or non-existent folder supplied, exiting"
	exit 1
elif ! [[ ${type} =~ $number ]]; then
	echo "Type is not a number or is empty. Please input 1 or 2 (1 for just out/err, 2 for out/err/sh)...exiting"
	exit 2
fi

# Clears out mass qsub script folder of all .sh, .err, and .out files and the complete folders within each qsub type folder
if [[ ${type} -eq 2 ]]; then
	for folder in ${location}/*
	do
		if [[ -d ${folder} ]]; then
			for folder1 in ${folder}
			do
				echo "Deleting every .sh script in ${folder1}"
				rm ${folder1}/*.sh
				rm ${folder1}/*.out
				rm ${folder1}/*.err

				echo "Deleting every .txt in ${folder1}/complete"
				rm ${folder1}/complete/*.txt
			done
		fi
	done
elif [[ ${1} = "-h" ]]; then
	echo "Will clean out script directory and mass qsub folders depending on parameter. 1 for qsub folder, 2, script folder, or 3 for both"
fi

# Deletes all straggling .err and .out files left in the home shareScript directory
if [[ ${type} -eq 1 ]]; then
	rm ${location}/blast16sID_*.out
	rm ${location}/blast16sID_*.err
	rm ${location}/blast16s_*.out
	rm ${location}/blast16s_*.err
	rm ${location}/ani_*.out
	rm ${location}/ani_*.err
	rm ${location}/ANI_*.out
	rm ${location}/ANI_*.err
	rm ${location}/BTQC_*.out
	rm ${location}/BTQC_*.err
	rm ${location}/BUSCO_*.out
	rm ${location}/BUSCO_*.err
	rm ${location}/getFASTQR1_*.out
	rm ${location}/getFASTQR1_*.err
	rm ${location}/getFASTQR2_*.out
	rm ${location}/getFASTQR2_*.err
	rm ${location}/csstn_*.out
	rm ${location}/csstn_*.err
	rm ${location}/csstp_*.out
	rm ${location}/csstp_*.err
	rm ${location}/kraka_*.out
	rm ${location}/kraka_*.err
	rm ${location}/krakr_*.out
	rm ${location}/krakr_*.err
	rm ${location}/gott_*.out
	rm ${location}/gott_*.err
	rm ${location}/mlst_*.out
	rm ${location}/mlst_*.err
	rm ${location}/pFinf_*.out
	rm ${location}/pFinf_*.err
	rm ${location}/pFinp_*.out
	rm ${location}/pFinp_*.err
	rm ${location}/SPAdn_*.out
	rm ${location}/SPAdn_*.err
	rm ${location}/SPAdp_*.out
	rm ${location}/SPAdp_*.err
	rm ${location}/plasFlow_*.out
	rm ${location}/plasFlow_*.err
	rm ${location}/pFlow_*.out
	rm ${location}/pFlow_*.err
	rm ${location}/pFinf_*.out
	rm ${location}/pFinf_*.err
	rm ${location}/PROKK_*.out
	rm ${location}/PROKK_*.err
	rm ${location}/PROKK_*.e*
	rm ${location}/srst2AR_*.out
	rm ${location}/srst2AR_*.err
	rm ${location}/srst2MLST_*.out
	rm ${location}/srst2MLST_*.err
	rm ${location}/srst22MLST_*.out
	rm ${location}/srst22MLST_*.err
	rm ${location}/QUAST_*.out
	rm ${location}/QUAST_*.err
	rm ${location}/QC_*.out
	rm ${location}/QC_*.err
	rm ${location}/MLST_*.out
	rm ${location}/MLST_*.err
	rm ${location}/taxID_*.out
	rm ${location}/taxID_*.err
	rm ${location}/validate_*.out
	rm ${location}/validate_*.err
	rm ${location}/sum_*.out
	rm ${location}/sum_*.err
	rm ${location}/pFn_*.out
	rm ${location}/pFn_*.err
	rm ${location}/pFp_*.out
	rm ${location}/pFp_*.err
	rm ${location}/aniB_*.out
	rm ${location}/aniB_*.err
	rm ${location}/aniM_*.out
	rm ${location}/aniM_*.err
	rm ${location}/node_*.out
	rm ${location}/node_*.err
	rm ${location}/core.*
	rm ${location}/quaisar_*.out
	rm ${location}/quaisar_*.err
	rm ${location}/SNVPhyl_*.out
	rm ${location}/SNVPhyl_*.err
fi
