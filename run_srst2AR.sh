#!/bin/sh -l

#$ -o srst2.out
#$ -e srst2.err
#$ -N srst2
#$ -cwd
#$ -q short.q

#
# Description: Script to use srst2 to attempt to find AR genes in parllel with assembly searches by other tools. This uses an alternate DB of genes
#
# Usage: ./run_srst2AR.sh -n sample_name -p MiSeq_Run_ID [-c	path_to_config_file] [-d path_to_alternate_DB]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/srst2/
#
# Modules required: srst2/0.2.0 bowtie2/2.2.4(?) Python2/2.7.13
#
# v1.0.3 (9/2/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml srst2 bowtie2/2.2.4


#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_srst2AR.sh -n sample_name -p MiSeq_Run_ID [-c	path_to_config_file] [-d path_to_alternate_DB]"
}

# Parse command line options
options_found=0
while getopts ":h?c:p:n:d:" option; do
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

if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
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

database_path="${ResGANNCBI_srst2}"
#database_and_version="${ResGANNCBI_srst2_filename}"

if [[ -z "${project}" ]]; then
	echo "No Project/Run_ID supplied to run_c-sstar_altDB.sh, exiting"
	exit 33
elif [[ -z "${sample_name}" ]]; then
	echo "No sample name supplied to run_c-sstar_altDB.sh, exiting"
	exit 34
elif [[ ! -z "${alt_db}" ]]; then
	if [[ ! -f "${alt_db}" ]]; then
		echo " No or empty alternate database location supplied to run_c-sstar_altDB.sh, exiting"
		exit 39
	else
		database_path="${alt_DB}"
#		database_basename=$(basename -- "${alt_db}")
#		database_basename2=$(echo ${database_basename##*/} | cut -d'.' -f1)
#		database_and_version=${database_basename2//_srst2/}
	fi
fi


OUTDATADIR="${processed}/${project}/${sample_name}"

# Create output directory
if [[ ! -d "${OUTDATADIR}/srst2" ]]; then
	echo "Creating ${OUTDATADIR}/srst2"
	mkdir "${OUTDATADIR}/srst2"
fi


# Prep any read files that have not been trimmed yet
if [ ! -f "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz" ]; then
	if [ -f "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]; then
		#echo "1"
		cp "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq.gz" "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz"
	elif [ -f "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" ]; then
		#echo "2"
		gzip -c "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" > "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz"
		#gzip < "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" > "${OUTDATADIR}/trimmed/${sample_name}_S1_L001_R1_001.fastq.gz"
	elif [[ ! -d "${OUTDATADIR}/trimmed" ]]; then
		#echo "5"
		if [[ -f "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq.gz" ]] && [[ ! -f "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq" ]]; then
			gunzip -c "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq.gz" > "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq"
		fi
		if [[ -f "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq.gz" ]] && [[ ! -f "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq" ]]; then
			gunzip -c "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq.gz" > "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq"
		fi
		ml BBMap/38.26 trimmomatic/0.35
		bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq" in2="${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq" out="${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R1.fsq" out2="${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
		trimmomatic "${trim_endtype}" -"${trim_phred}" -threads "${procs}" "${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R1.fsq" "${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R2.fsq" "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R1_001.unpaired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R2_001.unpaired.fq" ILLUMINACLIP:"${trim_adapter_location}:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome}" SLIDINGWINDOW:"${trim_window_size}:${trim_window_qual}" LEADING:"${trim_leading}" TRAILING:"${trim_trailing}" MINLEN:"${trim_min_length}"
		ml -BBMap/38.26 -trimmomatic/0.35
		cat "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" > "${OUTDATADIR}/trimmed/${sample_name}.paired.fq"
		cat "${OUTDATADIR}/trimmed/${sample_name}_R1_001.unpaired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R2_001.unpaired.fq" > "${OUTDATADIR}/trimmed/${sample_name}.single.fq"
		gzip -c "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" > "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz"
		gzip -c "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" > "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz"
	fi
fi
if [ ! -f "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz" ]; then
	if [ -f "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]; then
		#echo "3"
		cp "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq.gz" "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz"
	elif [ -f "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" ]; then
		#echo "4"
		gzip -c "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" > "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz"
		#gzip < "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" > "${OUTDATADIR}/trimmed/${sample_name}_S1_L001_R2_001.fastq.gz"
	fi
fi

# Prints the command that will be submitted to use srst2 to find AR genes
echo "--input_pe ${OUTDATADIR}/trimmed/${sample_name}_S1_L001_R1_001.fastq.gz ${OUTDATADIR}/trimmed/${sample_name}_S1_L001_R2_001.fastq.gz --output ${OUTDATADIR}/srst2/${sample_name}_ResGANNCBI --gene_db ${database_path}"

# Calls srst2 with the options for AR discovery
#python2 ${shareScript}/srst2-master/scripts/srst2.py --input_pe "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz" "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz" --output "${OUTDATADIR}/srst2/${sample_name}_ResGANNCBI" --gene_db "${ResGANNCBI_srst2}"
srst2 --input_pe "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz" "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz" --output "${OUTDATADIR}/srst2/${sample_name}_ResGANNCBI" --gene_db "${database_path}"

# Cleans up leftover files
rm -r "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz"
rm -r "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz"
rm -r "${OUTDATADIR}/srst2/"*".bam"
rm -r "${OUTDATADIR}/srst2/"*".pileup"

# Removes the extra ResGANNCBI__ from all files created
find ${OUTDATADIR}/srst2 -type f -name "*ResGANNCBI__*" | while read FILE ; do
  dirname=$(dirname $FILE)
	filename=$(basename $FILE)
	filename="${filename/_ResGANNCBI__/__}"
	#echo "Found-${FILE}"
	#echo "${filename}"
    mv "${FILE}" "${dirname}/${filename}"
done

ml -srst2 #-bowtie2/2.2.4
