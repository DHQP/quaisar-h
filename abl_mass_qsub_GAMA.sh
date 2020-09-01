#!/bin/sh -l

#$ -o ablmq_GAMA.out
#$ -e ablmq_GAMA.err
#$ -N ablmq_GAMA
#$ -cwd
#$ -q all.q

#
# Description: A script to submit a list of isolates to the cluster to perform GAMA on many isolates in parallel
#
# ./abl_mass_qsub_GAMA.sh -l path_to_list -m max_concurrent_submissions -o output_folder_for_scripts -k clobberness[keep|clobber] -c path_to_config_file(optional)
#
# Output location: default_config.sh_output_location/run_ID/sample_name/GAMA/. Temp scripts will be in default_mass_qsubs_folder_from_config.sh/GAMA_subs
#
# Modules required: None, run_GAMA.sh will load Python3/3.5.2 and blat/35
#
# v1.0.1 (08/18/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./abl_mass_qsub_GAMA.sh -l path_to_list -m max_concurrent_submissions -o output_folder_for_scripts -k clobberness[keep|clobber] -c path_to_config_file(optional)"
	echo "Output is saved to ${processed}/run_ID/sample_name/GAMA where processed is retrieved from config file, either default or imported"
}

# Parse command line options
options_found=0
while getopts ":h?l:s:m:s:k:c:" option; do
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
	echo "${script_output} directory parameter for script output is empty...exiting"
	exit 1
elif [[ -z "${cloberness}" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 4th parameter...exiting"
	exit 1
fi

# Check that clobberness is a valid option
if [[ "${clobberness}" != "keep" ]] && [[ "${clobberness}" != "clobber" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 5th parameter...exiting"
	exit 1
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
  arr+=("$line")
done < ${list}

arr_size="${#arr[@]}"
last_index=$(( arr_size -1 ))
echo "-${arr_size}:${arr[@]}-"

# Create counter and set max number of concurrent submissions
counter=0

"${shareScript}/clean_list.sh" "${list}"

# Set script directory
main_dir="${script_output}/GAMA_subs"
if [[ ! -d "${script_output}/GAMA_subs" ]]; then
	mkdir "${script_output}/GAMA_subs"
	mkdir "${script_output}/GAMA_subs/complete"
elif [[ ! -d  "${script_output}/GAMA_subs/complete" ]]; then
	mkdir "${script_output}/GAMA_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Creates and submits qsub scripts to check all isolates on the list against the newest ResGANNCBI DB
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ "${clobberness}" = "clobber" ]]; then
		rm ${processed}/${project}/${sample}/GAMA/${sample}_${ResGANNCBI_srst2_filename}.GAMA
	fi
	echo ${counter}
	# Check if counter is below max number of concurrent submissions
	if [ ${counter} -lt ${max_subs} ]; then
		# Check if the output file of GAMA exist, skip submission if so
		if [[ ! -f "${processed}/${project}/${sample}/GAMA/${sample}_${ResGANNCBI_srst2_filename}.GAMA" ]]; then
			echo  "Index is below max submissions, submitting"
			echo -e "#!/bin/bash -l\n" > "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "#$ -o GAMAAR_${sample}.out" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "#$ -e GAMAAR_${sample}.err" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "#$ -N GAMAAR_${sample}"   >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "cd ${shareScript}" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "\"${shareScript}/run_GAMA.sh\" -n \"${sample}\" -p \"${project}\" -k -c \"${config.sh}\"" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_GAMAAR_complete.txt\"" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			if [[ "${counter}" -lt "${last_index}" ]]; then
				qsub "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			else
				qsub -sync y "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
			fi
			mv "${shareScript}/GAMAAR_${sample}.err" ${main_dir}
			mv "${shareScript}/GAMAAR_${sample}.out" ${main_dir}
		# Old data existed, skipping
		else
			echo -e $(date) > "${main_dir}/complete/${sample}_GAMAAR_complete.txt"
			echo "${project}/${sample} already has newest GAMA ResGANNCBI ${ResGANNCBI_srst2_filename}"
		fi
	# Counter is above max submission, must wait for previous ones to finish before moving on
	else
		waiting_for_index=$(( counter - max_subs ))
		waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
		timer=0
		echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
		while :
		do
			# Check if timer is above max time allowed
			if [[ ${timer} -gt 1800 ]]; then
				echo "Timer exceeded limit of 1800 seconds 30 minutes"
				break
			fi
			# Check if waiting sample is finished
			if [ -f "${main_dir}/complete/${waiting_sample}_GAMAAR_complete.txt" ]; then
				# Check if current sample has etiher one of the output files from GAMA, skip analysis if so
				if [[ ! -f "${processed}/${project}/${sample}/GAMA/${sample}_${ResGANNCBI_srst2_filename}.GAMA" ]]; then
					echo  "Index is below max submissions, submitting"
					echo -e "#!/bin/bash -l\n" > "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "#$ -o GAMAAR_${sample}.out" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "#$ -e GAMAAR_${sample}.err" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "#$ -N GAMAAR_${sample}"   >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "cd ${shareScript}" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "\"${shareScript}/run_GAMA.sh\" -n \"${sample}\" -p \"${project}\" -k -c \"${config}\"" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_GAMAAR_complete.txt\"" >> "${main_dir}/GAMAAR_${sample}_${start_time}.sh"

					#cd "${main_dir}"
					if [[ "${counter}" -lt "${last_index}" ]]; then
						qsub "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					else
						qsub -sync y "${main_dir}/GAMAAR_${sample}_${start_time}.sh"
					fi
					mv "${shareScript}/GAMAAR_${sample}.err" ${main_dir}
					mv "${shareScript}/GAMAAR_${sample}.out" ${main_dir}
				# Old data existed, skipping
				else
					echo -e $(date) > "${main_dir}/complete/${sample}_GAMAAR_complete.txt"
					echo "${project}/${sample} already has newest GAMA ResGANNCBI ${ResGANNCBI_srst2_filename}"
				fi
				break
			# Wait 5 seconds and then check if "waiting" sample is complete
			else
				timer=$(( timer + 5 ))
				echo "sleeping for 5 seconds, so far slept for ${timer}"
				sleep 5
			fi
		done
	fi
	counter=$(( counter + 1 ))
done

# Check for completion of all samples
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_GAMAAR_complete.txt" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/GAMAAR_${waiting_sample}.out" ]]; then
			mv "${shareScript}/GAMAAR_${waiting_sample}.out" "${main_dir}"
		fi
		if [[ -f "${shareScript}/GAMAAR_${waiting_sample}.err" ]]; then
			mv "${shareScript}/GAMAAR_${waiting_sample}.err" "${main_dir}"
		fi
	else
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_GAMAAR_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/GAMAAR_${waiting_sample}.out" ]]; then
						mv "${shareScript}/GAMAAR_${waiting_sample}.out" "${main_dir}"
					fi
					if [[ -f "${shareScript}/GAMAAR_${waiting_sample}.err" ]]; then
						mv "${shareScript}/GAMAAR_${waiting_sample}.err" "${main_dir}"
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

exit 0
