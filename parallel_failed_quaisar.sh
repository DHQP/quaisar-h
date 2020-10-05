#!/bin/bash -l

#$ -o parfqua.out
#$ -e parfqua.err
#$ -N parqfua
#$ -cwd
#$ -q short.q

#
# Description: This version of the script on ly needs to know the location of the list to retry the pipeline on. All pieces should already be in place from the previous attempt.
# Usage: ./parallel_failed_quaisar.sh -l /list_of_samples_to_redo [-c path_to_config_file]
# Output location: default_config.sh_output_location
#
# Modules required: None
#		*script must be run on cluster or grid scheduler machine
#
# v1.0.2 (08/24/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#


#  Function to print out help blurb
show_help () {
	echo "Usage is ./parallel_failed_quaisar.sh -l path_to_list_of_samples_to_rerun [-c path_to_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?l:c:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		l)
			echo "Option -l triggered, argument = ${OPTARG}"
			list_path=${OPTARG};;
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

# Show help info for when no options are given
if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
fi

# Checks for proper argumentation
if [[ ! -f "${list_path}" ]] || [[ -z "${list_path}" ]]; then
	echo "List empty or non-existent, exiting"
	exit 1
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



# Creates a copy of config file to use for each run of the script, in case there is a change in output locations
echo "${shareScript}"
config_counter=0
while true
do
	if [[ "${config_counter}" -gt ${max_quaisars} ]]; then
		echo "Already ${max_quaisars} parallel quaisar sets running, please wait until one finishes (or check script directory for any straggling config_X.sh files that may not be being used anymore)...exiting"
		exit 324
	fi
	if [[ ! -f "${shareScript}/config_${config_counter}.sh" ]]; then
		if [[ -f "${shareScript}/config_template.sh" ]]; then
			echo "Trying to copy config_template.sh to config_${config_counter}"
			cp "${shareScript}/config_template.sh" "${shareScript}/config_${config_counter}.sh"
			break
		else
			echo "config_template.sh does not exist, cannot copy and must exit..."
			exit 333
		fi
	else
		config_counter=$(( config_counter + 1 ))
	fi
done

#Print out which type of machine the script is running on (Biolinux or Aspen as an interactive session or node based)
if [ "${host}" = "biolinux" ]; then
	echo "Running pipeline on Biolinux"
elif [ "${host}" = "aspen_login" ]; then
	echo "Running pipeline on Aspen interactive node"
elif [[ "${host}" = "cluster"* ]]; then
	echo "Running pipeline on Aspen ${host}"
fi

# Checking BASH version
if [ "${BASH_VERSINFO}" -lt 4 ];
then
	echo "Sorry, you need at least bash-4.0 to run this script." >&2;
	exit 1;
fi

# Checks the arguments (more to come)
nopts=$#
global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
requestor=$(whoami)



# Loops through list file to create an array of all isolates to run through pipeline
declare -a file_list=()
declare -a run_list=()
while IFS= read -r file || [ -n "$file" ]; do
	echo "Found: ${file}"
	sample_id=$(echo "${file}" | cut -d'/' -f2)
	run_id=$(echo ${file} | cut -d'/' -f1)
	file_list+=("${run_id}/${sample_id}")
	found=0
	for item in "${run_list[@]}"; do
		if [[ "${item}" == "${run_id}" ]]; then
			found=1
			break
		fi
	done
	if [[ "${found}" -eq 0 ]]; then
		run_list+=("${run_id}")
	fi
done < "${list_path}"

# Displays number and names of files found to analyze
if [[ ${#file_list[@]} -gt 1 ]]; then
	echo "Will analyze these ${#file_list[@]} files: ${file_list[*]}"
elif [[ ${#file_list[@]} -eq 1 ]]; then
	echo "Will analyze this file: ${file_list[0]}"
else
	echo "No files found in ${list_path}"
fi

run_start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

#Each file in the list is checked individually for successful completion and added then added to the log for the run
mkdir -p "${Quaisar_H_log_directory}/failed_redo_on_${run_start_time}"
log_dir="${Quaisar_H_log_directory}/failed_redo_on_${run_start_time}"

#Get the time the run started to use as the identifier
outarray=()
echo "Run started at ${run_start_time}; Log directory will be ${log_dir}/failed_redo_on_${run_start_time}"
echo "Run started at ${run_start_time}" > "${log_dir}/failed_redo_on_${run_start_time}.log"
outarray+=("failed_redo started at ${run_start_time} and saved to failed_redo_on_${run_start_time}.log")


#Each file in the list is put through the full pipeline
for projfile in "${file_list[@]}";
do
	echo "${projfile}"
	file=$(echo "${projfile}" | cut -d'/' -f2 | tr -d '[:space:]')
	proj=$(echo "${projfile}" | cut -d'/' -f1 | tr -d '[:space:]')
	echo "${file} ${proj}"
	if [[ -f "${processed}/${proj}/${file}/${file}_pipeline_stats.txt" ]]; then
		rm "${processed}/${proj}/${file}/${file}_pipeline_stats.txt"
	fi
	if [[ ! -f ${shareScript}/qfa_${file}.sh ]]; then
		cp ${shareScript}/quaisar_failed_assembly_template.sh ${shareScript}/qfa_${file}.sh
		sed -i -e "s/qfa_X/qfa_${file}/g" "${shareScript}/qfa_${file}.sh"
		echo "Entering ${shareScript}/qfa_${file}.sh" "${file}" "${proj}" "${shareScript}/config_${config_counter}.sh"
		qsub "${shareScript}/qfa_${file}.sh" -n "${file}" -p "${proj}" -c "${shareScript}/config_${config_counter}.sh"
		echo "Created and ran qfa_${file}.sh"
	else
		echo "${shareScript}/qfa_${file}.sh already exists, will resubmit"
		qsub "${shareScript}/qha_${file}.sh" -n "${file}" -p "${proj}" -c "${shareScript}/config_${config_counter}.sh"
	fi
done

# Hold for completion of all submited single quaisars
for run_sample in "${file_list[@]}"; do
	waiting_sample=$(echo "${run_sample}" | cut -d'/' -f2)
	proj=$(echo "${run_sample}" | cut -d'/' -f1)
	if [[ -f "${processed}/${proj}/${waiting_sample}/${waiting_sample}_pipeline_stats.txt" ]]; then
		echo "${waiting_sample} is complete"
		mv ${shareScript}/qfa_${waiting_sample}.sh ${log_dir}
		mv ${shareScript}/qfa_${waiting_sample}.err ${log_dir}
		mv ${shareScript}/qfa_${waiting_sample}.out ${log_dir}
	else
		while :
		do
				if [[ ${timer} -gt 1440 ]]; then
					echo "Timer exceeded limit of 86400 seconds(24 hours), Must complete other steps manually for ${1}"
					exit 1
				fi
				if [[ -f "${processed}/${proj}/${waiting_sample}/${waiting_sample}_pipeline_stats.txt" ]]; then
					echo "${waiting_sample} is complete"
					mv ${shareScript}/quaisar_${waiting_sample}.sh ${log_dir}
					mv ${shareScript}/quaisar_${waiting_sample}.err ${log_dir}
					mv ${shareScript}/quaisar_${waiting_sample}.out ${log_dir}
					break
				else
					timer=$(( timer + 1 ))
					if [[ $(( timer % 5 )) -eq 0 ]]; then
						echo "Slept for ${timer} minutes so far"
					fi
					sleep 60
				fi
		done
	fi
done

# Get run summary info to send in an email
# Hold for completion of all submited single quaisars
for run_sample in "${run_list[@]}"; do
	"${shareScript}/run_sum.sh" -p "${run_sample}" -c "${config}"
done

# Add print time the run completed in the text that will be emailed
global_end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "Finished with ${1} at ${global_end_time}"
outarray+=("
${1} finished at ${global_end_time}")
exit

#Send email to submitter and Nick with run status
if [ "${requestor}" != "nvx4" ]; then
	echo "Sending summary email to ${requestor}@cdc.gov & nvx4@cdc.gov"
	printf "%s\\n" "${outarray}" | mail -s "Run Status for ${1}_on_${run_start_time}_run.log" "nvx4@cdc.gov"
	printf "%s\\n" "${outarray}" | mail -s "Run Status for ${1}_on_${run_start_time}_run.log" "${requestor}@cdc.gov"
else
	echo "Sending summary email to nvx4@cdc.gov"
	printf "%s\\n" "${outarray}" | mail -s "Run Status for ${1}_on_${run_start_time}_run.log" "nvx4@cdc.gov"
fi

# One final check for any dump files
if [ -n "$(find "${shareScript}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
	echo "Found core dump files at end of processing ${1} and attempting to delete"
	find "${shareScript}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
fi

# Copy the config file to the log directory so as not to hold up any future quaisar runs that count the number of config files present, but for some reason does not remove from script folder
if [[ -f "${shareScript}/config_${config_counter}.sh" ]]; then
	echo "Supposedly moving config file(config_${config_counter}.sh) to log directory ($log_dir)"
	mv "${shareScript}/config_${config_counter}.sh" "${log_dir}/config_${PROJECT}.sh"
fi

end_date=$(date "+%m_%d_%Y_at_%Hh_%Mm")
echo "Run ended at ${end_date}" >> "${log_dir}/failed_redo_on_${run_start_time}/${PROJECT}_on_${run_start_time}.log"
