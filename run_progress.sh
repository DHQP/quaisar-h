#!/bin/sh -l

#$ -o quapro.out
#$ -e quapro.err
#$ -N quapro
#$ -cwd
#$ -q short.q

#
# Description: Helper for the main quaisar script to allow easy visualization of progress of run
#
# Usage ./quaisar-progress.sh path_to_run_folder
#
# Output loction: screen
#
# Modules required: None
#
# v1.0.1 (05/14/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

run_tasks=11
tasks_per_isolate=29

#Change window size to match progress bars (assuming this will be running for status only)
printf '\e[8;6;140t'
printf '\e[2t' && sleep 1 && printf '\e[1t'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./quaisar_progress.sh path_to_run_folder"
	echo "Output is saved to path_to_run_folder"
	exit 0
elif [[ ! -d "${1}" ]]; then
	echo "Path (${1}) does not exist, exiting quaisar-progress.sh"
	exit 2
fi

run_path=${1}
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
	echo -en "\n\nProgress for run $1\n[${isolate_progress}]\t${current_Isolate_progress}%-${current_Isolate_name}-${iso_AA[${pro_Isolate_task_number}]}\n[${run_progress}]\t${total_progress}%-${run_AA[${pro_run_task_id}]}\n\n"

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
