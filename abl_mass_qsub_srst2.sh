#!/bin/sh -l

#$ -o ablmq_srst2.out
#$ -e ablmq_srst2.err
#$ -N ablmq_srst2
#$ -cwd
#$ -q all.q

#
# Description: A script to submit a list of isolates to the cluster to perform SRST2 AR on many isolates in parallel
#
# Usage: Usage is ./abl_mass_qsub_srst2.sh -l path_to_list -m max_concurrent_submissions -o output_folder_for_scripts -k clobberness[keep|clobber] [-c path_to_config_file] [-d path_to_alt_DB]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/c-sstar/. Temp scripts will be in default_mass_qsubs_folder_from_config.sh/srst2_subs
#
# Modules required: None, run_srst2AR.sh will load srst2/0.2.0
#
# v1.0.2 (09/02/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./abl_mass_qsub_srst2.sh -l path_to_list -m max_concurrent_submissions -o output_folder_for_scripts -k clobberness[keep|clobber] [-c path_to_config_file] [-d path_to_alt_DB]"
	echo "Output is saved to *processed/run_ID/sample_name/srst2 where processed is retrieved from config file, either default or imported"
}

# Parse command line options
options_found=0
while getopts ":h?l:p:m:s:b:c:d:" option; do
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
		d)
			echo "Option -d triggered, argument = ${OPTARG}"
			alt_DB=${OPTARG};;
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

# Number regex to test max concurrent submission parameter
number='^[0-9]+$'

# Checks for proper argumentation
if [[ ! -f "${list}" ]]; then
	echo "${list} (list) does not exist...exiting"
	exit 1
elif ! [[ ${max_subs} =~ $number ]] || [[ -z "${max_subs}" ]]; then
	echo "${max_subs} is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
elif [[ -z "${script_output}" ]]; then
	echo "${script_output}, script output directory parameter is empty...exiting"
	exit 1
elif [[ -z "${clobberness}" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 4th parameter...exiting"
	exit 1
elif [[ ! -z "${alt_db}" ]]; then
	if [[ ! -f "${alt_db}" ]]; then
		echo " No or empty alternate database location supplied to run_c-sstar_altDB.sh, exiting"
		exit 39
	else
		use_alt_db="true"
		database_path="${alt_DB}"
	fi
fi

# Check that clobberness is a valid option
if [[ "${clobberness}" != "keep" ]] && [[ "${clobberness}" != "clobber" ]]; then
	echo "Clobberness was not input, be sure to add keep or clobber as 5th parameter...exiting"
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

# Set script directory
main_dir="${script_output}/srst2_subs"
if [[ ! -d "${script_output}/srst2_subs" ]]; then
	mkdir "${script_output}/srst2_subs"
	mkdir "${script_output}/srst2_subs/complete"
elif [[ ! -d  "${script_output}/srst2_subs/complete" ]]; then
	mkdir "${script_output}/srst2_subs/complete"
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Creates and submits qsub scripts to check all isolates on the list against the newest ResGANNCBI DB
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)
	if [[ "${clobberness}" = "clobber" ]]; then
		rm ${processed}/${project}/${sample}/srst2/${sample}__fullgenes__${ResGANNCBI_srst2_filename}_srst2__results.txt
		rm ${processed}/${project}/${sample}/srst2/${sample}__genes__${ResGANNCBI_srst2_filename}_srst2__results.txt
	fi
	echo ${counter}
	# Check if counter is below max number of concurrent submissions
	if [ ${counter} -lt ${max_subs} ]; then
		#echo "if [[ ! -f ${processed}/${project}/${sample}/srst2/${sample}__genes__${ResGANNCBI_srst2_filename}_srst2__results.txt ]] || [[ ! -f ${processed}/${project}/${sample}/srst2/${sample}__fullgenes__${ResGANNCBI_srst2_filename}_srst2__results.txt ]]; then"
		# Check if either one of the output files of srst2 files exist, skip submission if so
		if [[ ! -f "${processed}/${project}/${sample}/srst2/${sample}__genes__${ResGANNCBI_srst2_filename}_srst2__results.txt" ]] || [[ ! -f "${processed}/${project}/${sample}/srst2/${sample}__fullgenes__${ResGANNCBI_srst2_filename}_srst2__results.txt" ]]; then
			echo  "Index is below max submissions, submitting"
			echo -e "#!/bin/bash -l\n" > "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -o srst2AR_${sample}.out" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -e srst2AR_${sample}.err" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -N srst2AR_${sample}"   >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -cwd"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			echo -e "#$ -q short.q\n"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			# Can we somehow consolidate into one srst2 analysis to do MLST/AR/SEROTYPE
			echo -e "cd ${shareScript}" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			if [[ "${use_alt_db}" == "true" ]]; then
				echo -e "\"${shareScript}/run_srst2AR.sh\" -n \"${sample}\" -p \"${project}\" -c \"${config}\" -d \"${database_path}\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			else
				echo -e "\"${shareScript}/run_srst2AR.sh\" -n \"${sample}\" -p \"${project}\" -c \"${config}\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			fi
			echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_srst2AR_complete.txt\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"

			#cd "${main_dir}"
			if [[ "${counter}" -lt "${last_index}" ]]; then
				qsub "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			else
				qsub -sync y "${main_dir}/srst2AR_${sample}_${start_time}.sh"
			fi
		# Old data existed, skipping
		else
			echo -e $(date) > "${main_dir}/complete/${sample}_srst2AR_complete.txt"
			echo "${project}/${sample} already has newest srst2 ResGANNCBI ${ResGANNCBI_srst2_filename}"
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
			if [ -f "${main_dir}/complete/${waiting_sample}_srst2AR_complete.txt" ]; then
				# Check if current sample has etiher one of the output files from srst2, skip analysis if so
				if [[ ! -f "${processed}/${project}/${sample}/srst2/${sample}__genes__${ResGANNCBI_srst2_filename}_srst2__results.txt" ]] && [[ ! -f "${processed}/${project}/${sample}/srst2/${sample}__fullgenes__${ResGANNCBI_srst2_filename}_srst2__results.txt" ]]; then
					echo "${waiting_sample} has completed, starting ${sample}"
					echo -e "#!/bin/bash -l\n" > "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -o srst2AR_${sample}.out" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -e srst2AR_${sample}.err" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -N srst2AR_${sample}"   >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -cwd"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					echo -e "#$ -q short.q\n"  >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					# Can we somehow consolidate into one srst2 analysis to do MLST/AR/SEROTYPE
					echo -e "cd ${shareScript}" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					if [[ "${use_alt_db}" == "true" ]]; then
						echo -e "\"${shareScript}/run_srst2AR.sh\" -n \"${sample}\" -p \"${project}\" -c \"${config}\" -d \"${database_path}\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					else
						echo -e "\"${shareScript}/run_srst2AR.sh\" -n \"${sample}\" -p \"${project}\" -c \"${config}\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					fi
					echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_srst2AR_complete.txt\"" >> "${main_dir}/srst2AR_${sample}_${start_time}.sh"
					cd "${main_dir}"
					qsub "${main_dir}/srst2AR_${sample}_${start_time}.sh"
				# Old data existed, skipping
				else
					echo -e $(date) > "${main_dir}/complete/${sample}_srst2AR_complete.txt"
					echo "${project}/${sample} already has newest srst2 ResGANNCBI ${ResGANNCBI_srst2_filename}"
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
	if [[ -f "${main_dir}/complete/${waiting_sample}_srst2AR_complete.txt" ]]; then
		echo "${item} is complete"
	else
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_srst2AR_complete.txt" ]]; then
					echo "${item} is complete"
					break
				else
					timer=$(( timer + 5 ))
					echo "sleeping for 5 seconds, so far slept for ${timer}"
					sleep 5
				fi
		done
	fi
done

rm "${main_dir}"/complete/*.txt

exit 0
