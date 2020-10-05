#!/bin/sh -l

#$ -o get_Reads_from_Instruments.out
#$ -e get_Reads_from_Instruments.err
#$ -N get_Reads_from_Instruments
#$ -cwd
#$ -q short.q

#
# Description: Will find all fastq.gz files within the given folder. It will move and rename them to the location that the pipeline will expect
#
# Usage: ./get_Reads_from_folder.sh -p run_ID -i folder_with_fastqs -f postfix_for_reads(1:_SX_L001_RX_00X.fastq.gz 2: _(R)X.fastq.gz 3: _RX_00X.fastq.gz 4: _SX_RX_00X.fastq.gz) -o output_directory [-c path_to_config_file]
#
# Output location: output_directory/run_ID
#
# Modules required: None
#
# v1.0.1 (08/18/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./get_reads_from_folder.sh -i input_folder -o output folder -f post-fix_ID -p project_ID [-c path_to_config_file]"
	echo "Output is saved to ${processed}/run_ID/ where processed is retrieved from config file, either default or imported"
}

# Parse command line options
options_found=0
while getopts ":h?i:p:o:c:f:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		i)
			echo "Option -i triggered, argument = ${OPTARG}"
			input_folder=${OPTARG};;
		p)
			echo "Option -p triggered, argument = ${OPTARG}"
			project=${OPTARG};;
		o)
			echo "Option -o triggered, argument = ${OPTARG}"
			output_folder=${OPTARG};;
		c)
			echo "Option -c triggered, argument = ${OPTARG}"
			config=${OPTARG};;
		f)
			echo "Option -f triggered, argument = ${OPTARG}"
			postfix_index=${OPTARG};;
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

number='^[0-9]+$'

# Checks for proper argumentation
if [[ -z "${input_folder}" ]] || [[ ! -d "${input_folder}" ]]; then
	echo "Empty/non-existent input folder (${input_folder}) supplied, exiting"
	exit 1
elif [[ -z "${output_folder}" ]] || [[ ! -d "${output_folder}" ]]; then
	echo "Empty/non-existent output folder (${output}) supplied, exiting"
	exit 1
elif [[ -z "${project}" ]]; then
	echo "Empty project ID supplied, can not place files correctly, must exit"
	exit 1
elif ! [[ ${postfix_index} =~ $number ]] || [[ -z "${postfix_index}" ]]; then
	echo "postfix index is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
fi

if [[ "${postfix_index}" -gt 4 ]] || [[ "${postfix_index}" -lt 1 ]]; then
	echo "postfix for reads is TOO high, only 4 options...1:_SX_L001_RX_00X.fastq.gz 2: _(R)X.fastq.gz 3: _RX_00X.fastq.gz 4: _SX_RX_00X.fastq.gz , exiting"
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

# Sets folder to where files will be downloaded to
OUTDATADIR="${output_folder}/${project}"
if [ ! -d "${OUTDATADIR}" ]; then
	echo "Creating $OUTDATADIR"
	mkdir -p "${OUTDATADIR}"
fi
if [ -f "${OUTDATADIR}/${project}_list.txt" ]; then
	mv "${OUTDATADIR}/${project}_list.txt" "${OUTDATADIR}/${project}_list_original.txt"
fi
out_list="${OUTDATADIR}/${project}_list.txt"

####### Set trailing match pattern (everything after R in filename, to call unzip once for each pair) ######
match="${postfix_index}"

# Goes through given folder
echo "${input_folder}"
for file in ${input_folder}/*
do
	# Check if file is a zipped reads file
	if [[ "${file}" = *.gz ]] || [[ "${file}" = *.fastq ]] && [[ "${file}" != *_L001_I1_001.fastq.gz ]]; then
		echo "filename: ${file}"
		# Gets full file name from path
		full_sample_name=${file##*/}
		if [[ "${full_sample_name}" = *"_R1_"* ]]; then
			full_sample_name_pair=${full_sample_name/_R1_/_R2_}
		elif [[ "${full_sample_name}" = *"_R2_"* ]]; then
			full_sample_name_pair="${full_sample_name/_R2_/_R1_}"
		elif [[ "${full_sample_name}" = *"_1.fast"* ]]; then
			full_sample_name_pair=${full_sample_name/_1.fast/_2.fast}
		elif [[ "${full_sample_name}" = *"_2.fast"* ]]; then
			full_sample_name_pair="${full_sample_name/_2.fast/_1.fast}"
		fi
		# gets path from file
		source_path=$(dirname "${file}")
		# Extracts filename keeping only isolate ID, if it matches standard miseq naming
		echo "${match}:${full_sample_name}"
		if [[ "${match}" -eq 1 ]]; then
			short_name=$(echo "${full_sample_name}" | rev | cut -d'_' -f5- | rev)
			postfix=$(echo "${full_sample_name}" | rev | cut -d'_' -f1,2,3 | rev)
			# Create an array out of the full sample name, delimited by _
			#IFS='_' read -r -a name_array <<< "${full_sample_name}"
			#long_name=${name_array[0]}_${name_array[1]}_${name_array[2]}_${name_array[3]}
		elif [[ "${match}" -eq 4 ]]; then
			short_name=$(echo "${full_sample_name}" | rev | cut -d'_' -f4- | rev)
			postfix=$(echo "${full_sample_name}" | rev | cut -d'_' -f1,2,3 | rev)
		elif [[ "${match}" -eq 3 ]]; then
			short_name=$(echo "${full_sample_name}" | rev | cut -d'_' -f3- | rev)
			postfix=$(echo "${full_sample_name}" | rev | cut -d'_' -f1,2 | rev)
		elif [[ "${match}" -eq 2 ]]; then
			short_name=$(echo "${full_sample_name}" | rev | cut -d'_' -f2- | rev)
			postfix=$(echo "${full_sample_name}" | rev | cut -d'_' -f1 | rev)
		else
			echo "Magic - should have never gotten here as this number does not match any of the input numbers... 1:_SX_L001_RX_00X.fastq.gz 2: _(R)X.fastq.gz 3: _RX_00X.fastq.gz 4: _SX_RX_00X.fastq.gz , exiting"
			exit
		fi

		#long_name=$(echo "${full_sample_name}" | cut -d'_' -f1,2,3)
		echo "Short: ${short_name}"
		#echo "Does ${full_sample_name} match *${match}"

		# Skip file if it happens to be undetermined
    if [[ "${short_name}" == "Undetermined" ]]; then
			echo "found undetermined (${file})"
			continue
		# If the file matches the postfix given in the arguments proceed with moving and unzipping to the output directory
		else
			# Creates output folder
			if [ ! -d "${OUTDATADIR}/${short_name}" ]; then
				echo "Creating $OUTDATADIR/${short_name}"
				mkdir -p "${OUTDATADIR}/${short_name}"
				echo "Creating $OUTDATADIR/${short_name}/FASTQs"
				mkdir -p "${OUTDATADIR}/${short_name}/FASTQs"
			fi
			# Announces name of file being unzipped and then unzips it to the FASTQs folder for the matching sample name. Files are shortened to just name_R1_001.fastq or name_R2_001.fastq
			echo "Retrieving ${source_path}/${full_sample_name} and ${full_sample_name_pair}"
			#if [[ "${match}" -eq 1 ]] || [[ "${match}" -eq 4 ]]; then
				if [[ "${postfix}" = *"R1_001.fast"* ]] || [[ "${postfix}" = *"R1.fast"* ]] || [[ "${postfix}" = *"1.fast"* ]]; then
					if [[ "${full_sample_name}" = *".fastq.gz" ]]; then
						if [[ -f "${source_path}/${full_sample_name_pair}" ]]; then
							if [[ -f "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" ]] && [[ -f "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" ]]; then
								echo "${short_name} already has both zipped FASTQs (Probably done when R2 was found, this is the R1 tester)"
							else
								echo "Moving ${source_path}/${full_sample_name} and ${source_path}/${full_sample_name_pair} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz and ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
								#clumpify.sh in1="${source_path}/${full_sample_name}" in2="${source_path}/${full_sample_name_pair}" out1="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" out2="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" reorder=p
								cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
								cp "${source_path}/${full_sample_name_pair}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
							fi
						else
							echo "No R2 found for ${source_path}/${full_sample_name}"
							#clumpify.sh in="${source_path}/${full_sample_name}" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
							cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
						fi
					elif [[ "${full_sample_name}" = *".fastq" ]]; then
						if [[ -f "${source_path}/${full_sample_name_pair}" ]]; then
							if [[ -f "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" ]] && [[ -f "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" ]]; then
								echo "${short_name} already has both unzipped FASTQs (Probably done when R2 was found, this is the R1 tester)"
							else
								#if [[ ! -d "${OUTDATADIR}/${short_name}/FASTQs/temp" ]]; then
								#	mkdir -p "${OUTDATADIR}/${short_name}/FASTQs/temp"
								#fi
								echo "Zipping and not clumping ${source_path}/${full_sample_name} and ${source_path}/${full_sample_name_pair} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz and ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
								gzip -c "${source_path}/${full_sample_name}" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
								gzip -c "${source_path}/${full_sample_name_pair}" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
								#clumpify.sh in1="${OUTDATADIR}/${short_name}/FASTQs/temp/${short_name}_R1_001.fastq.gz" in2="${OUTDATADIR}/${short_name}/FASTQs/temp/${short_name}_R2_001.fastq.gz" out1="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" out2="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" reorder=p
								#rm -r "${OUTDATADIR}/${short_name}/FASTQs/temp"
							fi
						else
							echo "No R2 found for ${source_path}/${full_sample_name}"
							gzip -c "${source_path}/${full_sample_name}" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
							#clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/temp/${short_name}_R1_001.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
							#rm -r "${OUTDATADIR}/${short_name}/FASTQs/temp"
						fi
					fi
				fi
				if [[ "${postfix}" = *"R2_001.fast"* ]] || [[ "${postfix}" = *"R2.fast"* ]] || [[ "${postfix}" = *"2.fast"* ]]; then
					if [[ "${full_sample_name}" = *".fastq.gz" ]]; then
						if [[ -f "${source_path}/${full_sample_name_pair}" ]]; then
							if [[ -f "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" ]] && [[ -f "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" ]]; then
								echo "${short_name} already has both zipped FASTQs (Probably done when R1 was found, this is the R2 tester)"
							else
								echo "Not Clumping ${source_path}/${full_sample_name} and ${source_path}/${full_sample_name_pair} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz and ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
								#clumpify.sh in1="${source_path}/${full_sample_name_pair}" in2="${source_path}/${full_sample_name}" out1="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" out2="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" reorder=p
								cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
								cp "${source_path}/${full_sample_name_pair}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
							fi
						else
							echo "No R1 found for ${source_path}/${full_sample_name}"
							#clumpify.sh in="${source_path}/${full_sample_name}" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
							cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
						fi
					elif [[ "${full_sample_name}" = *".fastq" ]]; then
						if [[ -f "${source_path}/${full_sample_name_pair}" ]]; then
							if [[ -f "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" ]] && [[ -f "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" ]]; then
								echo "${short_name} already has both zipped FASTQs (Probably done when R1 was found, this is the R2 tester)"
							else
								#if [[ ! -d "${OUTDATADIR}/${short_name}/FASTQs/temp" ]]; then
								#	mkdir -p "${OUTDATADIR}/${short_name}/FASTQs/temp"
								#fi
									echo "Zipping and not clumping ${source_path}/${full_sample_name} and ${source_path}/${full_sample_name_pair} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz and ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
									gzip -c "${source_path}/${full_sample_name}" > "${OUTDATADIR}/${short_name}/FASTQs/temp/${short_name}_R1_001.fastq.gz"
									gzip -c "${source_path}/${full_sample_name_pair}" > "${OUTDATADIR}/${short_name}/FASTQs/temp/${short_name}_R2_001.fastq.gz"
								#	clumpify.sh in1="${OUTDATADIR}/${short_name}/FASTQs/temp/${short_name}_R1_001.fastq.gz" in2="${OUTDATADIR}/${short_name}/FASTQs/temp/${short_name}_R2_001.fastq.gz" out1="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" out2="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" reorder=p
								#	rm -r "${OUTDATADIR}/${short_name}/FASTQs/temp"
							fi
						else
							echo "No R1 found for ${source_path}/${full_sample_name}"
							gzip -c "${source_path}/${full_sample_name}" > "${OUTDATADIR}/${short_name}/FASTQs/temp/${short_name}_R2_001.fastq.gz"
							#clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/temp/${short_name}_R2_001.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
							#rm -r "${OUTDATADIR}/${short_name}/FASTQs/temp"
						fi
					fi
				fi
				if grep -Fxq "${project}/${short_name}" "${out_list}"
				then
					echo -e "${project}/${short_name} already on list ${out_list}, not adding again"
				else
					echo -e "${project}/${short_name}" >> "${out_list}"
				fi
		fi
	else
		echo "${file} is not a FASTQ(.gz) read file, not acting on it"
	fi
done

# Invert list so that the important isolates (for us at least) get run first
if [[ -f "${out_list}" ]]; then
	sort -k2,2 -t'/' -r "${out_list}" -o "${out_list}"
fi

#Script exited gracefully (unless something else inside failed)
exit 0
