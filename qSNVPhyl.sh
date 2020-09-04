#!/bin/sh -l

#$ -o qSNVPhyl.out
#$ -e qSNVPhyl.err
#$ -N qsnv
#$ -cwd
#$ -q short.q

#
# Description: The wrapper script that runs the run_SNVPhyl_template.sh script on the cluster, allowing multiple instances to run different sets of isolates
#
# Usage: ./qSNVPhyl -l path_to_list_file -o output_directory -p project_identifier [-c path_to_config_file]
#
# Output location: Parameter
#
# Modules required: None
#
# v1.0.1 (08/24/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./qSNVPhyl.sh -l path_to_list -o output_directory -n tree_output_name"
	echo "Output is saved to output_directory/tree_output_namepath_to_folder"
}

options_found=0
while getopts ":h?l:n:o:c:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		l)
			echo "Option -l triggered, argument = ${OPTARG}"
			input=${OPTARG};;
		o)
			echo "Option -o triggered, argument = ${OPTARG}"
			outdir=${OPTARG};;
		n)
			echo "Option -n triggered, argument = ${OPTARG}"
			output_file=${OPTARG};;
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

# Not being run on cluster=no run
if [[ ${host} != "cluster"* ]]; then
	echo "No scheduling system, can not run qSNVPhyl.sh"
	#exit 1
fi
sleep 5

# Loop through and act on each sample name in the passed/provided list
counter=0
#echo "Starting counting"
while [[ -f ${shareScript}/run_SNVPhyl_${counter}.sh ]]; do
	counter=$(( counter + 1 ))
	echo ${counter}
	# Do we need to set a global variable for max SNVPhyls??? (instead of 10)
	if [[ ${counter} -ge 10 ]]; then
		#echo "Showing all jobs by MMB team..."
		qstat -u nvx4 -u nvd4 -u njr5 -u xku6 -u kqj9 -u nyx3 -u kbi5 -u yer1 -u vif0 -u hex1
		echo "Too many (10) SNVPhyls running already; please try again later (also check that there arent unused scripts in ${shareScript}/run_SNVPhyl_*.sh)"
		exit
	fi
done

# Copy and convert snvphyl template to a specific numbered version
cp ${shareScript}/run_SNVPhyl_template.sh ${shareScript}/run_SNVPhyl_temp.sh
sed -i -e "s/run_SNVPhyl/run_SNVPhyl_${counter}/g" "${shareScript}/run_SNVPhyl_temp.sh"
sed -i -e "s/SNVPhyl_X/SNVPhyl_${counter}/g" "${shareScript}/run_SNVPhyl_temp.sh"
mv ${shareScript}/run_SNVPhyl_temp.sh ${shareScript}/run_SNVPhyl_${counter}.sh
echo "${shareScript}/run_SNVPhyl_${counter}.sh $@"
echo "Created and ran run_SNVPhyl_${counter}.sh"

# Submit the new snyphyl run file
qsub -sync y "${shareScript}/run_SNVPhyl_${counter}.sh" "$@"
rm ${shareScript}/run_SNVPhyl_${counter}.sh

# Sent an email showing that the script has completed
submitter=$(whoami)
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
printf "%s %s" "qSNVPhyl.sh has completed running run_SNVPhyl_${counter}.sh" "${global_end_time}" | mail -s "SNVPhyl ${counter} analysis complete" "${submitter}@cdc.gov"

#Script exited gracefully (unless something else inside failed)
exit 0
