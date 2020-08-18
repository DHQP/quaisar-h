#!/bin/sh -l

#$ -o get_Reads_from_Instruments.out
#$ -e get_Reads_from_Instruments.err
#$ -N get_Reads_from_Instruments
#$ -cwd
#$ -q short.q

# Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Will find all fastq.gz files within the given folder. It will move and rename them to the location that the pipeline will expect
#
# Usage: ./get_Reads_from_folder.sh full_path_Output_directory folder_with_fastqs postfix_for_reads(1:_SX_L001_RX_00X.fastq.gz 2: _(R)X.fastq.gz 3: _RX_00X.fastq.gz 4: _SX_RX_00X.fastq.gz)
#
# Output location: output_location
#
# Modules required: None
#
# v1.0.1 (08/18/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

number='^[0-9]+$'

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty project name supplied to $0, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./get_Reads_from_folder.sh  full_path_Output_directory location_of_fastqs postfix_for_reads( 1: _SX_L001_RX_00X.fastq.gz 2: _(R)X.fastq.gz 3: _RX_00X.fastq.gz 4: _SX_RX_00X.fastq.gz)"
	echo "Output by default is downloaded to ${processed}/run_ID and extracted to ${processed}/run_ID/sample_name/FASTQs"
	exit 0
elif [[ -z "${2}" ]]; then
	echo "Empty folder supplied to $0, exiting"
	exit 1
elif ! [[ ${3} =~ $number ]] || [[ -z "${3}" ]]; then
	echo "postfix is not a number or is empty. Please input max number of concurrent qsub submissions...exiting"
	exit 2
fi

if [[ "${3}" -gt 4 ]]; then
	echo "postfix for reads is TOO high, only 4 options...1:_SX_L001_RX_00X.fastq.gz 2: _(R)X.fastq.gz 3: _RX_00X.fastq.gz 4: _SX_RX_00X.fastq.gz , exiting"
	exit
fi


# Sets folder to where files will be downloaded to
OUTDATADIR="${1}"
if [ ! -d "${OUTDATADIR}" ]; then
	echo "Creating $OUTDATADIR"
	mkdir -p "${OUTDATADIR}"
fi
if [ -f "${OUTDATADIR}/${1}_list.txt" ]; then
	mv "${OUTDATADIR}/${1}_list.txt" "${OUTDATADIR}/${1}_list_original.txt"
fi
out_list="${OUTDATADIR}/${1}_list.txt"

####### Set trailing match pattern (everything after R in filename, to call unzip once for each pair) ######
match="${3}"

# Goes through given folder
echo "${2}"
for file in ${2}/*
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
				if grep -Fxq "${1}/${short_name}" "${out_list}"
				then
					echo -e "${1}/${short_name} already on list ${out_list}, not adding again"
				else
					echo -e "${1}/${short_name}" >> "${out_list}"
				fi
			# elif [[ "${match}" -eq 3 ]]; then
			# 	if [[ "${postfix}" = *"1.fast"* ]]; then
			# 		if [[ "${full_sample_name}" = *".fastq.gz" ]]; then
			# 			echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
			# 			cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
			# 			#clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
			# 		elif [[ "${full_sample_name}" = *".fastq" ]]; then
			# 			echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
			# 			gzip -c "${source_path}/${full_sample_name}" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
			# 			#clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
			# 		fi
			# 		echo -e "${1}/${short_name}" >> "${out_list}"
			# 	elif [[ "${postfix}" = *"2.fast"* ]]; then
			# 		if [[ "${full_sample_name}" = *".fastq.gz" ]]; then
			# 			echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
			# 			cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
			# 			#clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
			# 		elif [[ "${full_sample_name}" = *".fastq" ]]; then
			# 			echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
			# 			gzip -c "${source_path}/${full_sample_name}" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
			# 			#clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
			# 		fi
			# 	fi
			# elif [[ "${match}" -eq 2 ]]; then
			# 	if [[ "${postfix}" = *"R1.fast"* ]]; then
			# 		if [[ "${full_sample_name}" = *".fastq.gz" ]]; then
			# 			echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
			# 			cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
			# 			#clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
			# 		elif [[ "${full_sample_name}" = *".fastq" ]]; then
			# 			echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
			# 			gzip -c "${source_path}/${full_sample_name}" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
			# 			#clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R1_001.fastq.gz"
			# 		fi
			# 		echo -e "${1}/${short_name}" >> "${out_list}"
			# 	elif [[ "${postfix}" = *"R2.fast"* ]]; then
			# 		if [[ "${full_sample_name}" = *".fastq.gz" ]]; then
			# 			echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
			# 			cp "${source_path}/${full_sample_name}" "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
			# 			#clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
			# 		elif [[ "${full_sample_name}" = *".fastq" ]]; then
			# 			echo "${source_path}/${full_sample_name} to ${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
			# 			gzip -c "${source_path}/${full_sample_name}" > "${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
			# 			#clumpify.sh in="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz" out="${OUTDATADIR}/${short_name}/FASTQs/${short_name}_R2_001.fastq.gz"
			# 		fi
			# 	fi
			# else
			# 	echo "Unrecognized postfix type, but how did it get this far?"
			# fi
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
