#!/bin/bash -l

#$ -o sample_cleaner.out
#$ -e sample_cleaner.err
#$ -N sample_cleaner
#$ -cwd
#$ -q short.q

#
# Description: Script uses Gulviks SPAdes cleaner along with general folder cleanup to decrease footprint of samples after processing
#
# Usage ./sample_cleaner.sh -n sample_name -p run_ID [-c path_config_file]
#
# Output location: No output created
#
# Modules required: None
#
# v1.0.2 (09/08/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage: ./sample_cleaner.sh -n sample_name -p run_ID [-c path_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?c:p:n:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		p)
			echo "Option -p triggered, argument = ${OPTARG}"
			project=${OPTARG};;
		n)
			echo "Option -n triggered, argument = ${OPTARG}"
			sample_name=${OPTARG};;
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
if [[ -z "${sample_name}" ]]; then
	echo "Empty sample name supplied to run_kraken.sh, exiting"
	exit 1
elif [ -z "${project}" ]; then
	echo "Empty project name given. Exiting"
	exit 1
fi


# Set main sample folder to clean
sample_folder="${processed}/${project}/${sample_name}"
echo "Cleaning ${sample_folder}"
echo "Source - ${sample_folder}"
sample_name=$(echo "${sample_folder}" | rev | cut -d'/' -f1 | rev)
echo "Sample_ID - ${sample_name}"
# Remove the localANIDB from ANI output folder, if found
echo "Cleaning srst2"
if [[ -f "${sample_folder}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz" ]]; then
	rm -r "${sample_folder}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz"
fi
if [[ -f "${sample_folder}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz" ]]; then
	rm -r "${sample_folder}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz"
fi

echo "Cleaning ANI"
if [ -d "${sample_folder}/ANI/localANIDB" ]; then
	echo "removing localANIDb"
	rm -r "${sample_folder}/ANI/localANIDB"
fi
if [ -d "${sample_folder}/ANI/localANIDB_REFSEQ" ]; then
	echo "removing localANIDB_REFSEQ"
	rm -r "${sample_folder}/ANI/localANIDB_REFSEQ"
fi
if [ -d "${sample_folder}/ANI/localANIDB_full" ]; then
	echo "removing localANIDB_full"
	rm -r "${sample_folder}/ANI/localANIDB_full"
fi
if [ -d "${sample_folder}/ANI/temp" ]; then
	echo "removing temp"
	rm -r "${sample_folder}/ANI/temp"
fi
# Remove the hmmer output from the BUSCO folder
echo "Cleaning BUSCO"
if [ -d "${sample_folder}/BUSCO/hmmer_output" ]; then
	echo "removing hmmer output"
	rm -r "${sample_folder}/BUSCO/hmmer_output"
fi
# Use Gulviks cleaner script on regular SPAdes output
echo "Cleaning Assembly Folder"
if [ -d "${sample_folder}/Assembly" ]; then
	echo "Using Gulviks SPAdes cleaner on Assembly"
	${shareScript}/gulvic_SPAdes_cleaner.sh "${sample_folder}/Assembly"
fi
# Use Gulviks cleaner script on regular SPAdes output
echo "Cleaning plasmidAssembly Folder"
if [ -d "${sample_folder}/plasmidAssembly" ]; then
	echo "Using Gulviks SPAdes cleaner on plasmidAssembly"
	${shareScript}/gulvic_SPAdes_cleaner.sh "${sample_folder}/plasmidAssembly"
fi
echo "Cleaning GOTTCHA Folder"
# Remove splitrim fodler from gottcha output, if found
if [ -d "${sample_folder}/gottcha/gottcha_S/${sample_name}_temp/splitrim" ]; then
	echo "Deleting splitrim folder"
	rm -r "${sample_folder}/gottcha/gottcha_S/${sample_name}_temp/splitrim"
fi
# Removed intermediate folder that has reads with no adapters, but have not been trimmed yet
echo "Cleaning Adapter Folder"
if [ -d "${sample_folder}/removedAdapters" ]; then
	echo "Deleting adapterless reads"
	rm "${sample_folder}/removedAdapters/*.fsq"
	rm "${sample_folder}/removedAdapters/*.fsq.gz"
fi
# Clean trimmed folder of catted and unpaired reads, leaving only R1 and R2
echo "Cleaning Trimmed Folder"
if [ -d "${sample_folder}/trimmed" ]; then
	echo "Deleting extraneous reads"
	if [ -f "${sample_folder}/trimmed/${sample_name}.paired.fq" ]; then
		if [ ! -f "${sample_folder}/trimmed/${sample_name}.paired.fq.gz" ]; then
			gzip "${sample_folder}/trimmed/${sample_name}.paired.fq"
		fi
		echo "Deleting catted paired reads"
		rm "${sample_folder}/trimmed/${sample_name}.paired.fq"
	fi
	#if [ -f "${sample_folder}/trimmed/${sample_name}.single.fq" ]; then
	#	echo "Deleting catted single reads"
	#	rm "${sample_folder}/trimmed/${sample_name}.single.fq"
	#fi
	if [ ! -f "${sample_folder}/trimmed/${sample_name}.single.fq.gz" ]; then
		if  [ ! -f "${sample_folder}/trimmed/${sample_name}.single.fq" ]; then
			if [ -f "${sample_folder}/trimmed/${sample_name}_R1_001.unpaired.fq" ] && [ -f "${sample_folder}/trimmed/${sample_name}_R2_001.unpaired.fq" ]; then
				cat "${sample_folder}/trimmed/${sample_name}_R1_001.unpaired.fq"  "${sample_folder}/trimmed/${sample_name}_R2_001.unpaired.fq" > "${sample_folder}/trimmed/${sample_name}.single.fq"
			fi
			gzip "${sample_folder}/trimmed/${sample_name}.single.fq"
		fi
	fi
	if [ -f "${sample_folder}/trimmed/${sample_name}_R1_001.unpaired.fq" ]; then
		echo "Deleting unpaired R1 reads"
		rm "${sample_folder}/trimmed/${sample_name}_R1_001.unpaired.fq"
	fi
	if [ -f "${sample_folder}/trimmed/${sample_name}_R2_001.unpaired.fq" ]; then
		echo "Deleting unpaired R2 reads"
		rm "${sample_folder}/trimmed/${sample_name}_R2_001.unpaired.fq"
	fi
fi
# Clean FASTQ folder by zipping any unzipped reads
	# if [ -s "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq.gz" ] && [ -s "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq.gz"]; then
	# 	echo "Clumping R1 and R2"
	#
	# 	if [ -f "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq" ]; then
	# 		rm "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
	# 	fi
	# 	if [ -f "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq" ]; then
	# 		rm "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
	# 	fi
	if [ -f "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq" ]; then
		#echo "Found unzipped FASTQ"
		if [ ! -f "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq.gz" ]; then
			echo "Zipping R1"
			gzip -c "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
			if [ -s "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq.gz" ]; then
				echo "Zipping seems to have been successful, deleting ${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
				rm -r "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
			fi
		else
			echo "Zipped file found, checking for substance"
			if [ -s "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq.gz" ]; then
				echo "Zipped file not zero, deleting ${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
				rm -r "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
			else
				echo "Zip file was empty, trying to rezip"
				gzip -c "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
				if [ -s "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq.gz" ]; then
					echo "Rezip worked,Deleting ${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
					rm -r "${sample_folder}/FASTQs/${sample_name}_R1_001.fastq"
				else
					echo "Rezip did not work, what to do now???"
				fi
			fi
		fi
	fi
	if [ -f "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq" ]; then
	#echo "Found unzipped FASTQ"
	if [ ! -f "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq.gz" ]; then
		echo "Zipping R2"
		gzip -c "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
		if [ -s "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq.gz" ]; then
			echo "Zipping seems to have been successful, deleting ${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
			rm -r "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
		fi
	else
		echo "Zipped file found, checking for substance"
		if [ -s "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq.gz" ]; then
			echo "Zipped file not zero, deleting ${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
			rm -r "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
		else
			echo "Zip file was empty, trying to rezip"
			gzip -c "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
			if [ -s "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq.gz" ]; then
				echo "Rezip worked,Deleting ${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
				rm -r "${sample_folder}/FASTQs/${sample_name}_R2_001.fastq"
			else
				echo "Rezip did not work, what to do now???"
				fi
			fi
		fi
	fi
	# Zip trimmed R1 read, if not already done
	if [ -f "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq" ]; then
		#echo "Found unzipped FASTQ"
		if [ ! -f "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]; then
			echo "Zipping paired R1"
			gzip "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
			if [ -s "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]; then
				echo "Zipping seems to have been successful, deleting ${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
				rm -r "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
			fi
		else
			echo "Zipped file found, checking for substance"
			if [ -s "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]; then
				echo "Zipped file not zero, deleting ${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
				rm -r "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
			else
				echo "Zip file was empty, trying to rezip"
				gzip "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
				if [ -s "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]; then
					echo "Zipping seems to have been successful, deleting ${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
					rm -r "${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq"
				else
					echo "Rezip did not work, what to do now???"
				fi
			fi
		fi
	fi
	# Zip trimmed R2 read, if not already done
	if [ -f "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq" ]; then
		#echo "Found unzipped FASTQ"
		if [ ! -f "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]; then
			echo "Zipping paired R2"
			gzip "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
			if [ -s "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]; then
				echo "Zipping seems to have been successful, deleting ${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
				rm -r "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
			fi
		else
			echo "Zipped file found, checking for substance"
			if [ -s "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]; then
				echo "Zipped file not zero, deleting ${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
				rm -r "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
			else
				echo "Zip file was empty, trying to rezip"
				gzip "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
				if [ -s "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]; then
					echo "Zipping seems to have been successful, deleting ${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
					rm -r "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq"
				else
					echo "Rezip did not work, what to do now???"
				fi
			fi
		fi
	fi
	# Zip single fq files if present
	if [ -f "${sample_folder}/trimmed/${sample_name}.single.fq" ]; then
		if [ ! -f "${sample_folder}/trimmed/${sample_name}.single.fq.gz" ]; then
			gzip "${sample_folder}/trimmed/${sample_name}.single.fq.gz"
		fi
	fi
	#if [[ -f "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]] && [[ -f "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]]; then
	#	clumpify in1="${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq.gz" in2="${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" out1="${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq.clumped.gz" out2="${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.clumped.gz" reorder
	#fi
#echo "Sample ${project}/${sample_name} should now be clean" >> "${processed}/cleaned_sample_list.txt"
