#!/bin/bash -l

#$ -o sample_cleaner.out
#$ -e sample_cleaner.err
#$ -N sample_cleaner
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Script uses Gulviks SPAdes cleaner along with general folder cleanup to decrease footprint of samples after processing
#
# Usage ./sample_cleaner.sh   sample_name   run_ID
#
# Output location: No output created
#
# Modules required: None
#
# v1.0.1 (11/18/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#


# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./sample_cleaner.sh  sample_name MiSeq_Run_ID"
	echo "Will clean ${processed}/${2}/${1} folder"
	exit 0
elif [[ ! -d "${processed}/${2}" ]] || [[ ! -d "${processed}/${2}/${1}" ]]; then
	if [ ! -d "${processed}/${2}" ]; then
		echo "Project ${processed}/${2} does not exist"
	elif [ ! -d "${processed}/${2}/${1}" ]; then
		echo "Isolate ${processed}/${2}/${1} does not exist"
	fi
	echo "EXITING..."
	exit 1
else
	echo "Cleaning ${processed}/${2}/${1}"
fi

# Set main sample folder to clean
sample_folder="${processed}/${2}/${1}"
echo "Source - ${sample_folder}"
sample_name=$(echo "${sample_folder}" | rev | cut -d'/' -f1 | rev)
echo "Sample_ID - ${sample_name}"
# Remove the localANIDB from ANI output folder, if found
echo "Cleaning srst2"
if [[ -f "${sample_folder}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz" ]]; then
	rm -r "${sample_folder}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz"
fi
if [[ -f "${sample_folder}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz"]]; then
	rm -r "${sample_folder}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz"
fi

echo "Cleaning ANI"
if [ -d "${sample_folder}/ANI/localANIDB" ]; then
	echo "removing localANIDb"
	rm -r "${sample_folder}/ANI/localANIDB"
fi
if [ -d "${sample_folder}/ANI/localANIDB_full" ]; then
	echo "removing localANIDb_full"
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
	rm -r "${sample_folder}/removedAdapters"
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
	if [[ -f "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]] && [[ -f "${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]]; then
		clumpify in1="${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq.gz" in2="${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.gz" out1="${sample_folder}/trimmed/${sample_name}_R1_001.paired.fq.clumped.gz" out2="${sample_folder}/trimmed/${sample_name}_R2_001.paired.fq.clumped.gz" reorder
	fi
echo "Sample ${2}/${1} should now be clean" >> "${processed}/cleaned_sample_list.txt"
