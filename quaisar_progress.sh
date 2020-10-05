#!/bin/sh -l

#$ -o quapro.out
#$ -e quapro.err
#$ -N quapro
#$ -cwd
#$ -q short.q

#
# Description: Helper for the main quaisar script to allow easy visualization of progress of run
#
# Usage ./quaisar-progress.sh -i path_to_run_folder [-c path_to_config_file]
#
# Output loction: screen
#
# Modules required: None
#
# v1.0.2 (09/08/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

run_tasks=11
tasks_per_isolate=29

#Change window size to match progress bars (assuming this will be running for status only)
printf '\e[8;6;140t'
printf '\e[2t' && sleep 1 && printf '\e[1t'

#  Function to print out help blurb
show_help () {
	echo "Usage: ./quaisar-progress.sh -i path_to_run_folder [-c path_to_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?c:i:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		i)
			echo "Option -i triggered, argument = ${OPTARG}"
			run_path=${OPTARG};;
		c)
			echo "Option -v triggered, argument = ${OPTARG}"
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
if [ ! -d "${run_path}" ]; then
	echo "Run path (${run_path}) does not exist. Exiting"
	exit 1
fi

run_name=$(echo "${run_path}" | rev | cut -d'/' -f1 | rev)
BAR_length=100
BAR_character='#'
BAR=$(printf "[%${BAR_length}s]" | tr ' ' $BAR_character)

# run_task AA
declare -a run_AA=('Copying Reads/Assemblies to project directory' 'Inverting list' 'Listing all isolates' 'Displaying isolates' 'Creating unique run Identifier' 'Catting list' 'running isolates' 'Creating bug array' 'Creating Seqlog' 'Creating run summary' 'Copying config and closing out run')
# isolate task AA
declare -a iso_AA=('Prepping FASTQ folder' 'Raw Read Quality count' 'BBDUK PhiX' 'Trimmomatic' 'Trimmed Read Quality Count' 'Kraken on reads' 'GOTTCHA' 'SRST2 AR' 'SPAdes Assembling' 'Trimming Assemmbly' 'Kraken on Assembly' '16s Identification' 'Assembly QC' 'PROKKA' 'Rename Contig Headers' 'ANI' 'Taxon classification' 'BUSCO' 'c-SSTAR' 'GAMA' 'MLST' 'plasmidFinder' 'plasFlow' 'Check plasFlow assembly' 'c-SSTAR on plasFlow' 'plasmidFinder on PlasFlow' 'GAMA on plasFlow' 'Summarize isolate' 'Cleaning isolate')

while true; do
	pro_run_task_id=$(head -n1 ${run_path}/progress.txt | cut -d':' -f2)
	pro_Isolate_count=$(head -n2 ${run_path}/progress.txt | tail -n1 | cut -d':' -f2)
	current_Isolate_number=$(head -n3 ${run_path}/progress.txt | tail -n1 | cut -d':' -f2)
	isolate_index=$((current_Isolate_number + 1))
	current_Isolate_name=$(head -n${isolate_index} ${run_path}/${run_name}_list.txt | tail -n1 | cut -d'/' -f2)
	pro_Isolate_task_number=$(tail -n1 ${run_path}/progress.txt | cut -d':' -f2)
	total_jobs=$(( run_tasks + pro_Isolate_count * tasks_per_isolate ))
	#echo -e "${pro_Isolate_task_number}	${tasks_per_isolate}\n\n\n"
	current_Isolate_progress=$(( 100 * pro_Isolate_task_number / tasks_per_isolate ))
	jobs_completed=$(( current_Isolate_number * tasks_per_isolate + pro_run_task_id + pro_Isolate_task_number))
	total_progress=$(( 100 * jobs_completed / total_jobs ))
	#echo -e "${pro_run_task_id}	${pro_Isolate_count}	${current_Isolate_number}	${pro_Isolate_task_number}	${total_jobs}	${jobs_completed}\n\n\n"
	#echo "${current_Isolate_progress}"
	#echo "${total_progress}"
	isolate_incomplete_percent=$(( 100 - current_Isolate_progress ))
	total_incomplete_percent=$(( 100 - total_progress ))
	isolate_completed_string=$(printf "%0.s=" $(seq 1 ${current_Isolate_progress})) # Fill $variable with $n periods
	isolate_incomplete_string=$(printf "%0.s " $(seq 1 ${isolate_incomplete_percent})) # Fill $variable with $n periods
	total_completed_string=$(printf "%0.s=" $(seq 1 ${total_progress})) # Fill $variable with $n periods
	total_incomplete_string=$(printf "%0.s " $(seq 1 ${total_incomplete_percent}))
	isolate_progress="${isolate_completed_string}${isolate_incomplete_string}"
	run_progress="${total_completed_string}${total_incomplete_string}"
	#echo -e "${current_Isolate_progress}+${isolate_incomplete_percent}=100?"
	#echo -e "${total_progress}+${total_incomplete_percent}=100?"
	clear
	echo -en "\n\nProgress for run $run_name\n[${isolate_progress}]\t${current_Isolate_progress}%-${current_Isolate_name}-${iso_AA[${pro_Isolate_task_number}]}\n[${run_progress}]\t${total_progress}%-${run_AA[${pro_run_task_id}]}\n\n"

	#echo -ne "\r${BAR:0:$current_Isolate_progress}(${current_Isolate_progress}%-${current_Isolate_name}-${iso_AA[${pro_Isolate_task_number}]})"
	#	echo -ne "\r${BAR:0:$total_progress}(${total_progress}%-${run_AA[${pro_run_task_id}]})"
	sleep 1
	if [[ "${total_progress}" -eq 100 ]] || [[ ${pro_run_task_id} -eq 12 ]]; then
		echo "Run is complete!!!"
		printf '\e[8;24;80t'
		printf '\e[2t' && sleep 1 && printf '\e[1t'
		exit
	fi
done
