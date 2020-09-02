#!/bin/sh -l

#$ -o ablmq-animrq.out
#$ -e ablmq-animrq.err
#$ -N ablmq-animrq
#$ -cwd
#$ -q short.q


#
# Description: A script to submit a list of isolates to the cluster to perform csstar on many isolates in parallel
#
# Usage is ./abl_mass_qsub_ANI_REFSEQ.sh -l path_to_list -m max_concurrent_submissions -o output_folder_for_scripts -k clobberness[keep|clobber] [-c path_to_config_file] [-d path_to_alt_DB]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/ANI/. Temp scripts will be in default_mass_qsubs_folder_from_config.sh/c-sstar_subs
#
# Modules required: None, run_ANI_REFSEQ.sh will load Python3/3.5.2 and ncbi-blast+/LATEST
#
# v1.0.2 (09/02/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./abl_mass_qsub_csstar.sh -l path_to_list -m max_concurrent_submissions -o output_folder_for_scripts -k clobberness[keep|clobber] [-c path_to_config_file] [-d path_to_alt_DB]"
	echo "Output is saved to *processed/run_ID/sample_name/csstar where processed is retrieved from config file, either default or imported"
}

# Parse command line options
options_found=0
while getopts ":h?l:o:m:s:k:c:d:" option; do
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
		d)
			echo "Option -d triggered, argument = ${OPTARG}"
			alt_db=${OPTARG};;
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




# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Check that clobberness is a valid option
if [[ "${clobberness}" != "keep" ]] && [[ "${clobberness}" != "clobber" ]]; then
	echo "Clobberness was not input correctly [keep|clobber]...exiting"
	exit 1
fi

# Create an array of all samples in the list
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
main_dir="${script_output}/ANI_REFSEQ_subs"

if [[ ! -d "${main_dir}" ]]; then
	mkdir -p "${main_dir}/complete"
elif [[ ! -d "${main_dir}/complete" ]]; then
	mkdir "${main_dir}/complete"
fi
cp ./config.sh ${main_dir}

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

OUTDATADIR="${processed}/${project}/${sample}"

# Create and submit qsub scripts to get ANI for all isolates
while [ ${counter} -lt ${arr_size} ] ; do
	sample=$(echo "${arr[${counter}]}" | cut -d'/' -f2 | cut -d':' -f1)
	project=$(echo "${arr[${counter}]}" | cut -d'/' -f1)

	if [[ "${clobberness}" == "clobber" ]]; then
		echo "Trying to remove ${OUTDATADIR}/ANI/aniM_REFSEQ"
		rm -r "${OUTDATADIR}/ANI/aniM_REFSEQ"
		rm -r "${OUTDATADIR}/ANI/localANIDB_REFSEQ"
		rm -r "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample}_vs_REFSEQ_"*").txt"
	fi

	# Check if counter is below max submission limit
	if [[ ${counter} -lt ${max_subs} ]]; then
		# Check if old data exists, skip if so
		if [[ -s "${OUTDATADIR}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
			if [[ ! -f "${OUTDATADIR}/ANI/best_anim_hits_ordered(${sample}_vs_${database_and_version})" ]]; then
				echo  "Index is below max submissions, submitting"
				echo "Going to make ${main_dir}/anim_refseq_${sample}_${start_time}.sh"
				echo -e "#!/bin/bash -l\n" > "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
				echo -e "#$ -o anim_refseq_${sample}.out" >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
				echo -e "#$ -e anim_refseq_${sample}.err" >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
				echo -e "#$ -N anim_refseq_${sample}"   >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
				echo -e "#$ -cwd"  >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
				echo -e "#$ -q short.q\n"  >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
				echo -e "cd ${shareScript}" >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
				if [[ "${use_alt_db}" == "true" ]]; then
					echo -e "\"${shareScript}/run_ANI_REFSEQ.sh\" -n \"${sample}\" -p \"${project}\" -c \"${config}\" -d \"${database_path}\"" >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
				else
					echo -e "\"${shareScript}/run_ANI_REFSEQ.sh\" -n \"${sample}\" -p \"${project}\" -c \"${config}\"" >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
				fi
				echo -e "\"${shareScript}/determine_taxID.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
				echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_anim_refseq_complete.txt\"" >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
				cd "${main_dir}"
				qsub "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
			# Old data exists, skipping
			else
				echo "${project}/${sample} already has ANI summary"
				echo "$(date)" > "${main_dir}/complete/${sample}_anim_refseq_complete.txt"
			fi
		# No Assembly file to run ANI on, skipping
		else
			echo "${project}/${sample} does not have assembly"
			echo "$(date)" > "${main_dir}/complete/${sample}_anim_refseq_complete.txt"
		fi
	# Counter is above limit, wait until "slot" opens up"
	else
		waiting_for_index=$(( counter - max_subs ))
		waiting_sample=$(echo "${arr[${waiting_for_index}]}" | cut -d'/' -f2)
		timer=0
		echo "Index is above max submissions, waiting for index ${waiting_for_index}:${waiting_sample} to complete"
		# Loop to wait until "waiting" sample is complete
		while :
		do
			# If timer exceeeds limit then exit
			if [[ ${timer} -gt 1800 ]]; then
				echo "Timer exceeded limit of 1800 seconds 30 minutes"
				break
			fi
			# Check if "waiting" sample is complete
			if [[ -f "${main_dir}/complete/${waiting_sample}_anim_refseq_complete.txt" ]]; then
				# Check if an assembly exists to run ANI on
				if [[ -s "${OUTDATADIR}/Assembly/${sample}_scaffolds_trimmed.fasta" ]]; then
					# Check if old data exists, skip if so
					if [[ ! -f "${OUTDATADIR}/ANI/best_anim_hits_ordered(${sample}_vs_${genus,})" ]]; then
						echo  "${waiting_sample}(${waiting_for_index}) is not complete, submitting ${sample} ($counter)"
						echo -e "#!/bin/bash -l\n" > "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
						echo -e "#$ -o anim_refseq_${sample}.out" >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
						echo -e "#$ -e anim_refseq_${sample}.err" >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
						echo -e "#$ -N anim_refseq_${sample}"   >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
						echo -e "#$ -cwd"  >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
						echo -e "#$ -q short.q\n"  >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
						echo -e "cd ${shareScript}" >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
						echo -e "\"${shareScript}/run_ANI_REFSEQ.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
	  				echo -e "\"${shareScript}/determine_taxID.sh\" \"${sample}\" \"${project}\"" >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
						echo -e "echo \"$(date)\" > \"${main_dir}/complete/${sample}_anim_refseq_complete.txt\"" >> "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
						cd "${main_dir}"
						qsub "${main_dir}/anim_refseq_${sample}_${start_time}.sh"
					# Old data exists, skipping
					else
						echo "${project}/${sample} already has ANI summary"
						echo "$(date)" > "${main_dir}/complete/${sample}_anim_refseq_complete.txt"
					fi
					# No Assembly file to run ANI on, skipping
				else
					echo "${project}/${sample} does not have assembly"
					echo "$(date)" > "${main_dir}/complete/${sample}_anim_refseq_complete.txt"
				fi
				break
			# Wait 5 seconds before checking if "waiting" sample is complete
			else
				timer=$(( timer + 5 ))
				echo "sleeping for 5 seconds, so far slept for ${timer}"
				sleep 5
			fi
		done
	fi
	counter=$(( counter + 1 ))
done

# Loop to ensure all samples are complete (or time runs) before allowing the script to exit
timer=0
for item in "${arr[@]}"; do
	waiting_sample=$(echo "${item}" | cut -d'/' -f2)
	if [[ -f "${main_dir}/complete/${waiting_sample}_anim_refseq_complete.txt" ]]; then
		echo "${item} is complete"
		if [[ -f "${shareScript}/anim_refseq_${waiting_sample}.out" ]]; then
			mv "${shareScript}/anim_refseq_${waiting_sample}.out" ${main_dir}
		fi
		if [[ -f "${shareScript}/anim_refseq_${waiting_sample}.err" ]]; then
			mv "${shareScript}/anim_refseq_${waiting_sample}.err" ${main_dir}
		fi
	else
		# Check every 5 seconds to see if the sample has completed normal csstar analysis
		while :
		do
				if [[ ${timer} -gt 3600 ]]; then
					echo "Timer exceeded limit of 3600 seconds = 60 minutes"
					exit 1
				fi
				if [[ -f "${main_dir}/complete/${waiting_sample}_anim_refseq_complete.txt" ]]; then
					echo "${item} is complete"
					if [[ -f "${shareScript}/anim_refseq_$waiting_{sample}.out" ]]; then
						mv "${shareScript}/anim_refseq_${waiting_sample}.out" ${main_dir}
					fi
					if [[ -f "${shareScript}/anim_refseq_${waiting_sample}.err" ]]; then
						mv "${shareScript}/anim_refseq_${waiting_sample}.err" ${main_dir}
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
printf "%s %s" "abl_mass_qsub_ANIm_REFSEQ.sh has completed" "${global_end_time}" | mail -s "abl_mass_qsub_ANIm_REFSEQ.sh complete" nvx4@cdc.gov

#Script exited gracefully (unless something else inside failed)
exit 0
