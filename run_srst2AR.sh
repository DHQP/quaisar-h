#!/bin/sh -l

#$ -o srst2.out
#$ -e srst2.err
#$ -N srst2
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Script to use srst2 to attempt to find AR genes in parllel with assembly searches by other tools. This uses an alternate DB of genes
#
# Usage: ./run_srst2AR_altDB.sh   sample_name   MiSeq_Run_ID	path_to_alternate_DB
#
# Output location: default_config.sh_output_location/run_ID/sample_name/srst2/
#
# Modules required: srst2/0.2.0 bowtie2/2.2.4(?) Python2/2.7.13
#
# v1.0.1 (1/6/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml srst2 bowtie2/2.2.4

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./run_srst2.sh  sample_name run_ID"
	echo "Output location is ${processed}/run_ID/sample_name/srst2"
	exit 0
fi

# Create output directory
mkdir "${processed}/${2}/${1}/srst2"


# Prep any read files that have not been trimmed yet
if [ ! -f "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" ]; then
	if [ -f "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" ]; then
		#echo "1"
		cp "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
	elif [ -f "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" ]; then
		#echo "2"
		gzip -c "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
		#gzip < "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/trimmed/${1}_S1_L001_R1_001.fastq.gz"
	elif [[ ! -d "${processed}/${2}/${1}/trimmed" ]]; then
		#echo "5"
		if [[ -f "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq.gz" ]] && [[ ! -f "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq" ]]; then
			gunzip -c "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq.gz" > "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq"
		fi
		if [[ -f "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq.gz" ]] && [[ ! -f "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq" ]]; then
			gunzip -c "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq.gz" > "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq"
		fi
		ml BBMap/38.26 trimmomatic/0.35
		bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq" in2="${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq" out="${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R1.fsq" out2="${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
		trimmomatic "${trim_endtype}" -"${trim_phred}" -threads "${procs}" "${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R1.fsq" "${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R2.fsq" "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" "${processed}/${2}/${1}/trimmed/${1}_R1_001.unpaired.fq" "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" "${processed}/${2}/${1}/trimmed/${1}_R2_001.unpaired.fq" ILLUMINACLIP:"${trim_adapter_location}:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome}" SLIDINGWINDOW:"${trim_window_size}:${trim_window_qual}" LEADING:"${trim_leading}" TRAILING:"${trim_trailing}" MINLEN:"${trim_min_length}"
		ml -BBMap/38.26 -trimmomatic/0.35
		cat "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/trimmed/${1}.paired.fq"
		cat "${processed}/${2}/${1}/trimmed/${1}_R1_001.unpaired.fq" "${processed}/${2}/${1}/trimmed/${1}_R2_001.unpaired.fq" > "${processed}/${2}/${1}/trimmed/${1}.single.fq"
		gzip -c "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
		gzip -c "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
	fi
fi
if [ ! -f "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" ]; then
	if [ -f "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq.gz" ]; then
		#echo "3"
		cp "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
	elif [ -f "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" ]; then
		#echo "4"
		gzip -c "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
		#gzip < "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/trimmed/${1}_S1_L001_R2_001.fastq.gz"
	fi
fi

# Prints the command that will be submitted to use srst2 to find AR genes
echo "--input_pe ${processed}/${2}/${1}/trimmed/${1}_S1_L001_R1_001.fastq.gz ${processed}/${2}/${1}/trimmed/${1}_S1_L001_R2_001.fastq.gz --output ${processed}/${2}/${1}/srst2/${1}_ResGANNCBI --gene_db ${ResGANNCBI_srst2}"

# Calls srst2 with the options for AR discovery
#python2 ${shareScript}/srst2-master/scripts/srst2.py --input_pe "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" --output "${processed}/${2}/${1}/srst2/${1}_ResGANNCBI" --gene_db "${ResGANNCBI_srst2}"
srst2 --input_pe "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" --output "${processed}/${2}/${1}/srst2/${1}_ResGANNCBI" --gene_db "${ResGANNCBI_srst2}"

# Cleans up leftover files
rm -r "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
rm -r "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
rm -r "${processed}/${2}/${1}/srst2/"*".bam"
rm -r "${processed}/${2}/${1}/srst2/"*".pileup"

# Removes the extra ResGANNCBI__ from all files created
find ${processed}/${2}/${1}/srst2 -type f -name "*ResGANNCBI__*" | while read FILE ; do
  dirname=$(dirname $FILE)
	filename=$(basename $FILE)
	filename="${filename/_ResGANNCBI__/__}"
	#echo "Found-${FILE}"
	#echo "${filename}"
    mv "${FILE}" "${dirname}/${filename}"
done

ml -srst2 #-bowtie2/2.2.4
