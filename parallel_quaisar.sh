#!/bin/bash -l

#$ -o parqua.out
#$ -e parqua.err
#$ -N parqua
#$ -cwd
#$ -q all.q

#
# Description: The wrapper script that runs the QuAISAR-H pipeline in a grid scheduler environment.
# There are 2 ways to call the script. 1. If there is a default location and format that reads are kept then set that in the config file and use -p and the subfolder name containing the fastq files,
# or if the location is not in a standard place then use -i and the format number matching the reads postfix and set output directory with -o as follows
#
# Usage 1: ./parallel_quaisar.sh -p name_of_subfolder_within_default_read_location
# Usage 2: ./parallel_quaisar.sh -i path_to_reads_folder -f 1|2|3|4 -o path_to_parent_output_folder_location -p name_of_output_folder [-c path_to_config_file] [-a assemblies_only]
#
# Output location: Parameter or default_config.sh_output_location if p flag is used
#
# Modules required: None
#		*script must be run on cluster or grid scheduler machine
#
# v1.0.1 (08/24/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo -e "Usage 1: ./parallel_quaisar.sh -p MiSeq_Run_ID\n Usage 2: ./parallel_quaisar.sh -i path_to_reads_folder -f 1|2|3|4 -o path_to_parent_output_folder_location -n name_of_output_folder [-c path_to_config_file] [-a assemblies_only]"
}

change_processed="false"

# Parse command line options
options_found=0
while getopts ":h?c:p:i:f:o:n:a" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		p)
			echo "Option -p triggered, argument = ${OPTARG}"
			PROJECT=${OPTARG};;
		n)
			echo "Option -n triggered, argument = ${OPTARG}"
			named_project=${OPTARG};;
		i)
			echo "Option -i triggered, argument = ${OPTARG}"
			INDATADIR=${OPTARG}
			# Checks for proper argumentation
			if [[ -d  ${INDATADIR} ]]; then
				do_download="true"
			else
					echo "${INDATADIR} does not exist...exiting"
					exit 1
			fi
			indir_set="true";;
		f)
			echo "Option -f triggered, argument = ${OPTARG}"
			postfix=${OPTARG};;
		o)
			echo "Option -o triggered, argument = ${OPTARG}"
			new_processed=${OPTARG}
			change_processed="true";;
		c)
			echo "Option -c triggered, argument = ${OPTARG}"
			config=${OPTARG};;
		a)
			echo "Option -a triggered, argument = ${OPTARG}"
			assemblies="true";;
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
	if [[ "${config_counter}" -gt "${max_quaisars}" ]]; then
		echo "Already ${max_quaisars} parallel quaisar sets running, please wait until one finishes (or check script directory for any straggling config_X.sh files that may not be being used anymore)...exiting"
		exit 324
	fi
	if [[ ! -f "${shareScript}/config_${config_counter}.sh" ]]; then
		if [[ -f "${shareScript}/config_template.sh" ]]; then
			echo "Trying to copy ${shareScript}/config_template.sh to ${shareScript}/config_${config_counter}.sh"
			cp "${shareScript}/config_template.sh" "${shareScript}/config_${config_counter}.sh"
			echo "CP-${change_processed}"
			if [[ "${change_processed}" == "true" ]]; then
				echo "Changing processed location in config file"
				echo "processed=${new_processed}" >> "${shareScript}/config_${config_counter}.sh"
			fi
			if [[ -f "${shareScript}/config_${config_counter}.sh" ]]; then
				. "${shareScript}/config_${config_counter}.sh"
				break
			else
				echo "New config file was NOT created....find out why, exiting"
				exit
			fi
		else
			echo "config_template.sh does not exist, cannot copy and must exit..."
			exit 333
		fi
	else
		config_counter=$(( config_counter + 1 ))
	fi
done

BASEDIR="${processed}"

# Apply default values to other variables if a
if [[ ! -z "${PROJECT}" ]]; then
	is_proj="true"
	for instrument in "${all_instruments[@]}"
	do
		# Goes through each subfolder of the current instrument
		for run_folder in "${instrument}"/*
		do
			# Gets folder names in current directory
			run_ID=${run_folder##*/}
			#echo "${run_ID} - ${PROJECT}"
			# If folder name matches project name
			if [[ "${run_ID}" == "${PROJECT}" ]]; then
				# Print that a match was found
				echo "Found project ${run_ID} in ${instrument}"
				# Go through every file in the Basecalls folder of the found folder (all files will only be fastq.gzs)
				INDATADIR="${instrument}/${run_ID}/Data/Intensities/BaseCalls"
				break
			fi
		done
	done
	if [[ ! -d "${INDATADIR}" ]]; then
		echo "No FOLDER ${INDATADIR} exists on any MiSeq instrument, exiting"
		exit 123
	fi
	postfix=1
else
		if [[ -z "${named_project}" ]]; then
			echo "No project name given for the alternate output location"
		elif [[ -z "${postfix}" ]] || [[ ${postfix} -lt 1 ]] || [[ ${postfix} -gt 4 ]]; then
			echo "No read postfix, or it is outside the 1-4 range, to know how to deal with reads...exiting"
		fi
		PROJECT="${named_project}"

fi

list_path="${processed}/${PROJECT}/${PROJECT}_list.txt"



#Print out which type of machine the script is running on (Biolinux or Aspen as an interactive session or node based)
if [ "${host}" = "biolinux" ];
then
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

# Short print out summary of run settings
echo -e "Source folder: ${INDATADIR}\\nOutput folder: ${BASEDIR}\\nList based analysis:  ${list_path}"


# Checks that a full FASTQ source path is given
if [ "${INDATADIR:0:1}" != "/" ] && [ "${INDATADIR:0:1}" != "$" ]; then
	echo "${INDATADIR}"
	echo 'ERROR: The full path was not specified.' >&2
	exit 1
fi

# Copies reads from source location to working directory and creates a list of IDs
if [[ "${assemblies}" == "true" ]]; then
	"${shareScript}/get_Assemblies_from_folder.sh" -p "${PROJECT}" -i "${INDATADIR}" -o "${processed}" -c "${shareScript}/config_${config_counter}.sh"
else
	"${shareScript}/get_Reads_from_folder.sh" -p "${PROJECT}" -i "${INDATADIR}" -f "${postfix}" -o "${processed}" -c "${shareScript}/config_${config_counter}.sh"
fi

# Loops through list file to create an array of all isolates to run through pipeline
declare -a file_list=()
while IFS= read -r file || [ -n "$file" ]; do
	echo "Found: ${file}"
	file=$(echo "${file}" | tr -d '[:space:]')
	file_list+=("${file}")
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
mkdir -p "${Quaisar_H_log_directory}/${PROJECT}_on_${run_start_time}"
log_dir="${Quaisar_H_log_directory}/${PROJECT}_on_${run_start_time}"

#Get the time the run started to use as the identifier
outarray=()
echo "Run started at ${run_start_time}; Log directory will be ${Quaisar_H_log_directory}/${PROJECT}_on_${run_start_time}"
echo "Run started at ${run_start_time}" > "${log_dir}/${PROJECT}_on_${run_start_time}/${PROJECT}_on_${run_start_time}.log"
outarray+=("${PROJECT} started at ${run_start_time} and saved to ${PROJECT}_on_${run_start_time}.log\n")


#Each file in the list is put through the full pipeline
if [[ "${assemblies}" == "true" ]]; then
	for projfile in "${file_list[@]}";
	do
		echo "${projfile}"
		file=$(echo "${projfile}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
		proj=$(echo "${projfile}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
		echo "${file} ${proj} ${BASEDIR}"
		if [[ -f "${processed}/${proj}/${file}/${file}_pipeline_stats.txt" ]]; then
			rm "${processed}/${proj}/${file}/${file}_pipeline_stats.txt"
		fi
		if [[ ! -f ${shareScript}/quaisar_on_assembly_${file}.sh ]]; then
			cp ${shareScript}/quaisar_on_assembly_template.sh ${shareScript}/quaisar_on_assembly_${file}.sh
			sed -i -e "s/qoa_X/qoa_${file}/g" "${shareScript}/quaisar_on_assembly_${file}.sh"
			sed -i -e "s/qoaX/qoa_${file}/g" "${shareScript}/quaisar_on_assembly_${file}.sh"
			echo "Entering ${shareScript}/quaisar_on_assembly_${file}.sh" "${file}" "${proj}" "${shareScript}/config_${config_counter}.sh"
			qsub "${shareScript}/quaisar_on_assembly_${file}.sh" -n "${file}" -p "${proj}" -c "${shareScript}/config_${config_counter}.sh"
			echo "Created and ran quaisar_on_assembly_${file}.sh"
		else
			echo "${shareScript}/quaisar_on_assembly_${file}.sh already exists, will resubmit"
			qsub "${shareScript}/quaisar_on_assembly_${file}.sh" -n "${file}" -p "${proj}" -c "${shareScript}/config_${config_counter}.sh"
		fi
	done
else
	for projfile in "${file_list[@]}";
	do
		echo "${projfile}"
		file=$(echo "${projfile}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
		proj=$(echo "${projfile}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
		echo "${file} ${proj} ${BASEDIR}"
		if [[ -f "${processed}/${proj}/${file}/${file}_pipeline_stats.txt" ]]; then
			rm "${processed}/${proj}/${file}/${file}_pipeline_stats.txt"
		fi
		if [[ ! -f ${shareScript}/quaisar_${file}.sh ]]; then
			cp ${shareScript}/quaisar_template.sh ${shareScript}/quaisar_${file}.sh
			sed -i -e "s/quaisar_X/quaisar_${file}/g" "${shareScript}/quaisar_${file}.sh"
			sed -i -e "s/quasX/quasp_${file}/g" "${shareScript}/quaisar_${file}.sh"
			echo "Entering ${shareScript}/quaisar_${file}.sh" "${file}" "${proj}" "${shareScript}/config_${config_counter}.sh"
			qsub "${shareScript}/quaisar_${file}.sh" -n "${file}" -p "${proj}" -c "${shareScript}/config_${config_counter}.sh"
			echo "Created and ran quaisar_${file}.sh"
		else
			echo "${shareScript}/quaisar_${file}.sh already exists, will resubmit"
			qsub "${shareScript}/quaisar_${file}.sh" -n "${file}" -p "${proj}" -c "${shareScript}/config_${config_counter}.sh"
		fi
	done
fi

# Hold for completion of all submited single quaisars
for run_sample in "${file_list[@]}"; do
	waiting_sample=$(echo "${run_sample}" | cut -d'/' -f2)
	if [[ -f "${processed}/${PROJECT}/${waiting_sample}/${waiting_sample}_pipeline_stats.txt" ]]; then
		echo "${waiting_sample} is complete"
		mv ${shareScript}/quaisar_${waiting_sample}.sh ${log_dir}
		mv ${shareScript}/quaisar_${waiting_sample}.err ${log_dir}
		mv ${shareScript}/quaisar_${waiting_sample}.out ${log_dir}
	else
		while :
		do
				if [[ ${timer} -gt 1440 ]]; then
					echo "Timer exceeded limit of 86400 seconds(24 hours), Must complete other steps manually for ${PROJECT}"
					exit 1
				fi
				if [[ -f "${processed}/${PROJECT}/${waiting_sample}/${waiting_sample}_pipeline_stats.txt" ]]; then
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

# Concatenates lists if this run was an addition to an already processed folder
if [[ -f "${processed}/${PROJECT}/${PROJECT}_list_original.txt" ]]; then
	cat "${processed}/${PROJECT}/${PROJECT}_list.txt" >> "${processed}/${PROJECT}/${PROJECT}_list_original.txt"
	rm "${processed}/${PROJECT}/${PROJECT}_list.txt"
	mv "${processed}/${PROJECT}/${PROJECT}_list_original.txt" "${processed}/${PROJECT}/${PROJECT}_list.txt"
fi

# Run the Seqlog creator on the proper file
if [ "${is_proj}" = "true" ]; then
	"${shareScript}/make_Seqlog_from_log.sh" -p "${PROJECT}" -c "${shareScript}/config_${config_counter}.sh"
else
	"${shareScript}/make_Seqlog_from_list.sh" -l "${processed}/${PROJECT}/${PROJECT}_list.txt"  -c "${shareScript}/config_${config_counter}.sh"
fi

# Get run summary info to send in an email
runsumdate=$(date "+%m_%d_%Y_at_%Hh_%Mm")
${shareScript}/run_sum.sh -p ${PROJECT} -c "${shareScript}/config_${config_counter}.sh"
runsum=$(${shareScript}/view_sum.sh ${PROJECT})
outarray+="${runsum}"

# Add print time the run completed in the text that will be emailed
global_end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "Finished with run ${PROJECT} at ${global_end_time}"
outarray+=("
${PROJECT} finished at ${global_end_time}")

#Send email to submitter and Nick with run status
if [ "${requestor}" != "nvx4" ]; then
	echo "Sending summary email to ${requestor}@cdc.gov & nvx4@cdc.gov & kbi5@cdc.gov"
	printf "%s\\n" "${outarray}" | mail -s "Run Status for ${PROJECT}_on_${run_start_time}_run.log" "nvx4@cdc.gov"
	printf "%s\\n" "${outarray}" | mail -s "Run Status for ${PROJECT}_on_${run_start_time}_run.log" "kbi5@cdc.gov"
	printf "%s\\n" "${outarray}" | mail -s "Run Status for ${PROJECT}_on_${run_start_time}_run.log" "${requestor}@cdc.gov"
else
	echo "Sending summary email to nvx4@cdc.gov and kbi5@cdc.gov"
	printf "%s\\n" "${outarray}" | mail -s "Run Status for ${PROJECT}_on_${run_start_time}_run.log" "nvx4@cdc.gov"
	printf "%s\\n" "${outarray}" | mail -s "Run Status for ${PROJECT}_on_${run_start_time}_run.log" "kbi5@cdc.gov"
fi

# One final check for any dump files
if [ -n "$(find "${shareScript}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
	echo "Found core dump files at end of processing ${PROJECT} and attempting to delete"
	find "${shareScript}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
fi

# Copy the config file to the log directory so as not to hold up any future quaisar runs that count the number of config files present, but for some reason does not remove from script folder
if [[ -f "${shareScript}/config_${config_counter}.sh" ]]; then
	echo "Supposedly moving config file(config_${config_counter}.sh) to log directory ($log_dir)"
	cp "${shareScript}/config_${config_counter}.sh" "${processed}/${PROJECT}/config_${PROJECT}.sh"
	#mv "${shareScript}/config_${config_counter}.sh" "${log_dir}/config_${PROJECT}.sh"
fi

end_date=$(date "+%m_%d_%Y_at_%Hh_%Mm")
echo "Run ended at ${end_date}" >> "${log_dir}/${PROJECT}_on_${run_start_time}/${PROJECT}_on_${run_start_time}.log"
