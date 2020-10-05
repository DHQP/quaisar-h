#!/bin/sh -l

#$ -o srst2.out
#$ -e srst2.err
#$ -N srst2
#$ -cwd
#$ -q short.q

#
# Description: Script to use srst2 to attempt to find AR genes in parllel with assembly searches by other tools. This uses an alternate DB of genes
#
# Usage: ./run_srst2AR_altDB.sh   sample_name   MiSeq_Run_ID	path_to_alternate_DB path_to_config_file(optional)
#
# Output location: default_config.sh_output_location/run_ID/sample_name/srst2/
#
# Modules required: srst2/0.2.0 bowtie2/2.2.4(?) Python2/2.7.13
#
# v1.0.2 (8/18/2020)
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
	echo "Usage is ./run_srst2AR_altDB.sh.sh  sample_name MiSeq_Run_ID path_to_alt_DB"
	echo "Output location is ${processed}/run_ID/srst2"
	exit 0
elif [[ ! -f ${3} ]]; then
	echo "Alternate database ($3) does not exist, exiting"
fi

if [[ -f "${4}" ]];
	echo "Loading special config file - ${4}"
	. "${4}"
else
	echo "Loading default config file"
	if [[ ! -f "./config.sh" ]]; then
		cp ./config_template.sh ./config.sh
	fi
	. ./config.sh
fi

alt_DB_path=${3}
alt_DB=$(echo ${alt_DB_path##*/} | cut -d'.' -f1)
alt_DB=${alt_DB//_srst2/}

# Create output folder
if [[ ! -d "${processed}/${2}/${1}/srst2" ]]; then
	mkdir "${processed}/${2}/${1}/srst2"
fi

# Prep any reads that have not been trimmed
if [ ! -f "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" ]; then
	if [ -f "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" ]; then
		#echo "1"
		cp "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
	elif [ -f "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" ]; then
		#echo "2"
		$(gzip -c "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz")
		#gzip < "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/trimmed/${1}_S1_L001_R1_001.fastq.gz"
	fi
else
	echo "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz already exists"
fi
if [ ! -f "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" ]; then
	if [ -f "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq.gz" ]; then
		#echo "3"
		cp "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
	elif [ -f "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" ]; then
		#echo "4"
		$(gzip -c "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz")
		#gzip < "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/trimmed/${1}_S1_L001_R2_001.fastq.gz"
	fi
else
	echo "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz already exists"
fi

# Display command that will be submitted
echo "--input_pe ${processed}/${2}/${1}/trimmed/${1}_S1_L001_R1_001.fastq.gz ${processed}/${2}/${1}/trimmed/${1}_S1_L001_R2_001.fastq.gz --output ${processed}/${2}/${1}/srst2/${alt_DB} --gene_db ${alt_DB_path}"

# Submit command to run srst2 on alternate DB
#python2 ${shareScript}/srst2-master/scripts/srst2.py --input_pe "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" --output "${processed}/${2}/${1}/srst2/${1}_${alt_DB}" --gene_db "${alt_DB_path}"
srst2 --input_pe "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" --output "${processed}/${2}/${1}/srst2/${1}_${alt_DB}" --gene_db "${alt_DB_path}"

# Clean up extra files
rm -r "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
rm -r "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
rm -r "${processed}/${2}/${1}/srst2/"*".bam"
rm -r "${processed}/${2}/${1}/srst2/"*".pileup"

# Rename files to reduce redundancy
find ${processed}/${2}/${1}/srst2 -type f -name "*_${alt_DB}__*" | while read FILE ; do
	#echo "Doing - $FILE"
  dirname=$(dirname $FILE)
	filename=$(basename $FILE)
	filename="${filename/_${alt_DB}__/__}"
	#echo "${filename}"
  mv "${FILE}" "${dirname}/${filename}"
done

ml -srst2 -bowtie2/2.2.4
