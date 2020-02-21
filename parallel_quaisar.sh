#!/bin/bash -l

#$ -o parqua.out
#$ -e parqua.err
#$ -N parqua
#$ -cwd
#$ -q short.q

# Sets the sharescript variable temporarily to the current working directory, allowing it to find the original config.sh file
shareScript=$(pwd)
#Import the config file with shortcuts and settings
. ${shareScript}/config.sh

#
# Description: The wrapper script that runs the QuAISAR-H pipeline in a grid scheduler environment.
# There are 2 ways to call the script. 1. If there is a default location and format that reads are kept then set that in the config file and use -p and the subfolder name containing the fastq files,
# or if the location is not in a standard place then use -i and the format number matching the reads post_fix and set output directory with -o as follows
#
# Usage 1: ./parallel_quaisar.sh -p name_of_subfolder_within_default_read_location
# Usage 2: ./parallel_quaisar.sh -i path_to_reads_folder 1|2|3|4 -o path_to_parent_output_folder_location name_of_output_folder
#
# Output location: Parameter or default_config.sh_output_location if p flag is used
#
# Modules required: None
#		*script must be run on cluster or grid scheduler machine
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#


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
if [ "${host}" = "biolinux" ];
then
	echo "Running pipeline on Biolinux"
elif [ "${host}" = "aspen_login" ];
then
	echo "Running pipeline on Aspen interactive node"
elif [[ "${host}" = "cluster"* ]];
then
	echo "Running pipeline on Aspen ${host}"
fi

# Checking BASH version
if [ "${BASH_VERSINFO}" -lt 4 ];
then
	echo "Sorry, you need at least bash-4.0 to run this script." >&2;
	exit 1;
fi

# Checking for proper number of arguments from command line
if [[ $# -lt 1  || $# -gt 7 ]]; then
	echo "If reads are in default location set in config file then"
  echo "Usage: ./parallel_quaisar.sh -p project_name"
	echo "else if you are running it on reads not in the default location or format"
	echo "Usage: ./parallel_quaisar.sh -i location_of_reads 1|2|3|4 -o path_to_parent_output_folder_location name_of_output_folder [-a]"
	echo "filename postfix numbers are as follows 1:_SX_L001_RX_00X.fastq.gz 2: _(R)X.fastq.gz 3: _RX_00X.fastq.gz 4: _SX_RX_00X.fastq.gz"
  echo "You have used $# args"
  exit 3
fi

# Checks the arguments (more to come)
nopts=$#
do_download=false
global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
requestor=$(whoami)
PROJECT="${requestor}_${global_time}"
BASEDIR="${processed}"

for ((i=1 ; i <= nopts ; i++)); do
	#echo "${1} ${2}"
	case "${1}" in
		#Help/Usage section
		-h | --help)
			echo -e "\\n\\n\\n"
			echo "If reads are in default location set in config file then"
		  echo "Usage: ./parallel_quaisar.sh -p project_name"
			echo "else if you are running it on reads not in the default location or format"
			echo "Usage: ./parallel_quaisar.sh -i location_of_reads 1|2|3|4 -o path_to_parent_output_folder_location name_of_output_folder [-a]"
			echo "filename postfix numbers are as follows 1:_SX_L001_RX_00X.fastq.gz 2: _(R)X.fastq.gz 3: _RX_00X.fastq.gz 4: _SX_RX_00X.fastq.gz"
			echo -e "\\n\\n\\n"
			exit 0
			;;
		#Gets name of folder that FASTA files will be in
		-i | --in-dir)
			INDATADIR="$2"
			if [[ -d  ${INDATADIR} ]]; then
				do_download="true"
				list_path="${BASEDIR}/${PROJECT}/${PROJECT}_list.txt"
			else
					echo "${INDATADIR} does not exist...exiting"
					exit 1
			fi
			indir_set="true"
			postfix="$3"
			#is_full_run="false"
			#echo "$INDATADIR $2"
			shift 3
			;;
		#Gets output directory name of folder that all output files will be stored
		-o | --out-dir)
			BASEDIR="$2"
			PROJECT="$3"
			shift 3

			echo "processed=${BASEDIR}" >> "${shareScript}/config.sh"
			. ${shareScript}/config.sh
			echo "${processed}"
			list_path="${BASEDIR}/${PROJECT}/${PROJECT}_list.txt"
			if [[ ! -d ${BASEDIR} ]]; then
				mkdir -p ${BASEDIR}
			fi
			;;
		#Checks for (project) name of folder that all output files will be stored
		-p | --project-name)
			PROJECT="$2"
			is_proj="true"
			#Copies over and unzips all 'Determined' zipped FASTQ files from the sequencing instrument for the requested project_id
			for instrument in "${all_instruments[@]}"
			do
				# Goes through each subfolder of the current instrument
				for run_folder in "${instrument}"/*
				do
					# Gets folder names in current directory
					run_ID=${run_folder##*/}
					#echo "${run_ID} - ${PROJECT}"
					# If folder name matches project name
					if [[ "${run_ID}" = "${PROJECT}" ]]; then
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
			if [[ ${BASEDIR} = "${requestor}_${global_time}" ]]; then
				BASEDIR="${processed}"
			fi
			list_path="${processed}/${PROJECT}/${PROJECT}_list.txt"
			postfix=1
			shift 2
			;;
			# If the target files assemblies only
			-a | --assemblies)
				assemblies="true"
				shift
				;;
		#Captures any other characters in the args
		\?)
			echo "ERROR: ${BOLD}$2${NORM} is not a valid argument" >&2
			usage
			exit 1
			;;
	esac
done

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
	"${shareScript}/get_assemblies_from_folder.sh" "${PROJECT}" "${INDATADIR}"
else
	"${shareScript}/get_Reads_from_folder.sh" "${PROJECT}" "${INDATADIR}" "${postfix}"
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
outarray+=("${PROJECT} started at ${run_start_time} and saved to ${PROJECT}_on_${run_start_time}.log")


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
			qsub "${shareScript}/quaisar_on_assembly_${file}.sh" "${file}" "${proj}" "${shareScript}/config_${config_counter}.sh"
			echo "Created and ran quaisar_on_assembly_${file}.sh"
		else
			echo "${shareScript}/quaisar_on_assembly_${file}.sh already exists, will resubmit"
			qsub "${shareScript}/quaisar_on_assembly_${file}.sh" "${file}" "${proj}" "${shareScript}/config_${config_counter}.sh"
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
			qsub "${shareScript}/quaisar_${file}.sh" "${file}" "${proj}" "${shareScript}/config_${config_counter}.sh"
			echo "Created and ran quaisar_${file}.sh"
		else
			echo "${shareScript}/quaisar_${file}.sh already exists, will resubmit"
			qsub "${shareScript}/quaisar_${file}.sh" "${file}" "${proj}" "${shareScript}/config_${config_counter}.sh"
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
					echo "Timer exceeded limit of 86400 seconds(24 hours), Must complete other steps manually for ${project}"
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
	"${shareScript}/make_Seqlog_from_log.sh" "${PROJECT}"
else
	"${shareScript}/make_Seqlog_from_list.sh" "${processed}/${PROJECT}/${PROJECT}_list.txt"
fi

# Get run summary info to send in an email
runsumdate=$(date "+%m_%d_%Y_at_%Hh_%Mm")
${shareScript}/run_sum.sh ${PROJECT}
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
	mv "${shareScript}/config_${config_counter}.sh" "${log_dir}/config_${PROJECT}.sh"
fi

end_date=$(date "+%m_%d_%Y_at_%Hh_%Mm")
echo "Run ended at ${end_date}" >> "${log_dir}/${PROJECT}_on_${run_start_time}/${PROJECT}_on_${run_start_time}.log"
