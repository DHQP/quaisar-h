#!/bin/sh -l

#$ -o ablmq-cs.out
#$ -e ablmq-cs.err
#$ -N ablmq-cs
#$ -cwd
#$ -q all.q

#
# Description: A script to submit a list of isolates to the cluster to perform csstar on many isolates in parallel
#
# Usage is ./abl_mass_qsub_csstar.sh -l path_to_list -m max_concurrent_submissions -o output_folder_for_scripts -k clobberness[keep|clobber] -s %ID(optional)[80|95|98|99|100] -c path_to_config_file(optional)
#
# Output location: default_config.sh_output_location/run_ID/sample_name/c-sstar/. Temp scripts will be in default_mass_qsubs_folder_from_config.sh/c-sstar_subs
#
# Modules required: None, run_c-sstar.sh will load Python3/3.5.2 and ncbi-blast+/LATEST
#
# v1.0.1 (08/18/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./abl_mass_qsub_csstar.sh -l path_to_list -m max_concurrent_submissions -o output_folder_for_scripts -k clobberness[keep|clobber] -s %ID(optional)[80|95|98|99|100] -c path_to_config_file(optional)"
	echo "Output is saved to ${processed}/run_ID/sample_name/csstar where processed is retrieved from config file, either default or imported"
}

# Parse command line options
options_found=0
while getopts ":h?l:o:m:s:k:c:" option; do
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
		m)
			echo "Option -m triggered, argument = ${OPTARG}"
			max_subs=${OPTARG};;
		o)
			echo "Option -o triggered, argument = ${OPTARG}"
			script_output=${OPTARG};;
		k)
			echo "Option -k triggered, argument = ${OPTARG}"
			clobberness=${OPTARG};;
		s)
			echo "Option -s triggered, argument = ${OPTARG}"
			percent_ID=${OPTARG};;
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

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

if [[ ! -f "${list}" ]]; then
	echo "${list} (list) does not exist...exiting"
	exit 1
elif ! [[ ${max_subs} =~ $number ]] || [[ -z "${max_subs}" ]]; then
	echo "${max_subs} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
elif [[ -z "${script_output}" ]]; then
	echo "Script output directory parameter is empty...exiting"
	exit 1
elif [[ -z "${cloberness}" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 4th parameter...exiting"
	exit 1
fi

# Check that clobberness is a valid option
if [[ "${clobberness}" != "keep" ]] && [[ "${cloberness}" != "clobber" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 4th parameter...exiting"
	exit 1
fi

# Checks that value given for % Identity is one of the presets for csstar
if [[ "${percent_ID}" != 80 ]] && [[ "${percent_ID}" != 95 ]] && [[ "${percent_ID}" != 98 ]] && [[ "${percent_ID}" != 99 ]] && [[ "${percent_ID}" != 100 ]]; then
	echo "Identity is not one of the presets for csstar and therefore will fail, defaulting to 98..."
	sim="h"
	simnum=98
else
	if [ "${percent_ID}" == 98 ]; then
		sim="h"
	elif [ "${percent_ID}" == 80 ]; then
		sim="l"
	elif [ "${percent_ID}" == 99 ]; then
		sim="u"
	elif [ "${percent_ID}" == 95 ]; then
		sim="m"
	elif [ "${percent_ID}" == 100 ]; then
		sim="p"
	elif [ "${percent_ID}" == 40 ]; then
		sim="o"
	fi
	simnum=${percent_ID}
fi

if [[ -f "${config}" ]];
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

# create an array of all samples in the list
arr=()
while IFS= read -r line || [ "$line" ];  do
	if [[ ! -z "${line}" ]]; then
		line=$(echo ${line} | tr -d '\n' | tr -d '\r')
		arr+=($line)
	fi
done < ${list}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"


# Create direcory to hold all temporary qsub scripts
counter=0

# format name being extracted from alt database
main_dir="${script_output}/csstar_subs"
#cp ./config ${main_dir}
if [[ ! -d "${script_output}/csstar_subs" ]]; then
	mkdir -p "${script_output}/csstar_subs/complete"
elif [[ ! -d  "${script_output}/csstar_subs/complete" ]]; then
	mkdir "${script_output}/csstar_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Create and submit scripts to run default csstar on all samples on the list
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ "${clobberness}" = "clobber" ]]; then
		echo "Removing old c-sstar for ${project}/${sample} and ${ResGANNCBI_srst2_filename}"
		rm ${processed}/${project}/${sample}/c-sstar/${sample}.${ResGANNCBI_srst2_filename}.gapped_${simnum}_sstar_summary.txt
		rm -r ${processed}/${project}/${sample}/c-sstar/${ResGANNCBI_srst2_filename}_gapped/
	fi
	#echo ${counter}-${project}-${sample}
	# Check if sample has a usable assembly file
	if [[ -s "${processed}/${project}/${sample}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
		#echo "Test"
		# Check if counter is below max number of concurrent submissions
		if [[ ${counter} -lt ${max_subs} ]]; then
			# Check if old data exists, skip if so
			if [[ ! -f "${processed}/${project}/${sample}/c-sstar/${sample}.${ResGANNCBI_srst2_filename}.gapped_${simnum}_sstar_summary.txt" ]]; then
				echo  "Index is below max submissions, submitting"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -o csstn_${sample}.out" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -e csstn_${sample}.err" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -N csstn_${sample}"   >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				# Defaulting to gapped/98, change if you want to include user preferences
				echo -e "cd ${shareScript}" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "\"${shareScript}/run_c-sstar.sh\" -n \"${sample}\" -g gapped -s \"${sim}\" -p \"${project}\" -c \"${config}\""  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarn_complete.txt\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
				cd "${main_dir}"
				echo "submitting ${main_dir}/csstn_${sample}_${start_time}.sh"
				qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
			# Old data exists
			else
				echo "${project}/${sample} already has the newest ResGANNCBI (${ResGANNCBI_srst2_filename})"
				echo -e "$(date)" > "${main_dir}/complete/${sample}_csstarn_complete.txt"
			fi
		# Counter is above max number of submissions
		else
			waiting_for_index=$(( counter - max_subs ))
			waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
			timer=0
			echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
			while :
			do
				# Check that the timer has not exceeded max amount of time to wait
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 1800 seconds 30 minutes"
					break
				fi
				# Check if usable assembly exists for current sample or that one does not exist for the waiting sample (therefore there would be no need to wait on it)
				if [[ -f "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
					# Check if old data exists, skip if so
					if [[ ! -f "${processed}/${project}/${sample}/c-sstar/${sample}.${ResGANNCBI_srst2_filename}.gapped_${simnum}_sstar_summary.txt" ]]; then
						echo  "Index is below max submissions, submitting"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -o csstn_${sample}.out" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -e csstn_${sample}.err" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -N csstn_${sample}"   >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "module load Python/3.6.1\n" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						# Defaulting to gapped/98, change if you want to include user preferences
						echo -e "cd ${shareScript}" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "\"${shareScript}/run_c-sstar.sh\" -n \"${sample}\" -g gapped -s \"${sim}\" -p \"${project}\" -c \"${config}\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_csstarn_complete.txt\"" >> "${main_dir}/csstn_${sample}_${start_time}.sh"
						cd ${main_dir}
						qsub "${main_dir}/csstn_${sample}_${start_time}.sh"
					# Skipping because old data exists
					else
						echo "${project}/${sample} already has the newest ResGANNCBI (${ResGANNCBI_srst2_filename})"
						echo -e "$(date)" > "${main_dir}/complete/${sample}_csstarn_complete.txt"
					fi
					break
				# If waiting sample has not completed, wait 5 more seconds and try again
				else
					timer=$(( timer + 5 ))
					echo "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt not ready, sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
			done
		fi
	fi
	counter=$(( counter + 1 ))
done

# Loop to ensure all samples are complete (or time runs) before allowing the script to exit
timer=0
ptimer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt" ]] || [[ ! -s "${processed}/${project}/${waiting_sample}/Assembly/${waiting_sample}_scaffolds_trimmed.fasta" ]]; then
		echo "${item} is complete normal"
		if [[ -f "${shareScript}/csstn_${waiting_sample}.out" ]]; then
			mv "${shareScript}/csstn_${waiting_sample}.out" "${main_dir}"
		fi
		if [[ -f "${shareScript}/csstn_${waiting_sample}.err" ]]; then
			mv "${shareScript}/csstn_${waiting_sample}.err" "${main_dir}"
		fi
	else
		# Check every 5 seconds to see if the sample has completed normal csstar analysis
		while :
		do
				if [[ ${timer} -gt 1800 ]]; then
					echo "Timer exceeded limit of 1800 seconds = 30 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_csstarn_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/csstn_${waiting_sample}.out" ]]; then
						mv "${shareScript}/csstn_${waiting_sample}.out" "${main_dir}"
					fi
					if [[ -f "${shareScript}/csstn_${waiting_sample}.err" ]]; then
						mv "${shareScript}/csstn_${waiting_sample}.err" "${main_dir}"
					fi
					break
				else
					timer=$(( timer + 5 ))
					echo "sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
		done
	fi
done

echo "All isolates completed"
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
printf "%s %s" "abl_mass_qsub_csstar.sh has completed" "${global_end_time}" | mail -s "abl_mass_qsub_csstar.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
