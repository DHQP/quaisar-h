#!/bin/sh -l

#$ -o srst2.out
#$ -e srst2.err
#$ -N srst2
#$ -cwd
#$ -q short.q

# Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Script to use srst2 to attempt to find mlst profile on reads. Used if standard mlst profiling fails
#
# Usage: ./run_srst2_mlst.sh   sample_name   MiSeq_Run_ID	Genus	species
#
# Output location: default_config.sh_output_location/run_ID/sample_name/MLST/
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
	echo "Usage is ./run_srst2.sh  sample_name MiSeq_Run_ID Genus species"
	echo "Output location is ${processed}/${2}/srst2"
	exit 0
fi

species="${4,,}"
genus="${3,,}"
genus="${genus^}"

# Creates folder for output
if [[ ! -d "${processed}/${2}/${1}/srst2" ]]; then
	mkdir "${processed}/${2}/${1}/srst2"
fi

# Preps reads if there are no current trimmed reads in the folder
if [ ! -f "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" ]; then
	if [ -f "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" ]; then
		#echo "1"
		cp "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
	elif [ -f "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" ]; then
		#echo "2"
		gzip -c "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
	else
		#echo "3"
		if [[ -f "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq.gz" ]] && [[ ! -f "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq" ]]; then
			gunzip -c "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq.gz" > "${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq"
		fi
		if [[ -f "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq.gz" ]] && [[ ! -f "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq" ]]; then
			gunzip -c "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq.gz" > "${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq"
		fi
		echo "Running BBDUK and trimmomatic"
		ml BBMAP/38.26 trimmomatic/0.35
		bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${processed}/${2}/${1}/FASTQs/${1}_R1_001.fastq" in2="${processed}/${2}/${1}/FASTQs/${1}_R2_001.fastq" out="${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R1.fsq" out2="${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
		trimmomatic "${trim_endtype}" -"${trim_phred}" -threads "${procs}" "${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R1.fsq" "${processed}/${2}/${1}/removedAdapters/${1}-noPhiX-R2.fsq" "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" "${processed}/${2}/${1}/trimmed/${1}_R1_001.unpaired.fq" "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" "${processed}/${2}/${1}/trimmed/${1}_R2_001.unpaired.fq" ILLUMINACLIP:"${trim_adapter_location}:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome}" SLIDINGWINDOW:"${trim_window_size}:${trim_window_qual}" LEADING:"${trim_leading}" TRAILING:"${trim_trailing}" MINLEN:"${trim_min_length}"
		ml -BBMAP/38.26 -trimmomatic/0.35
		#cat "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/trimmed/${1}.paired.fq"
		#cat "${processed}/${2}/${1}/trimmed/${1}_R1_001.unpaired.fq" "${processed}/${2}/${1}/trimmed/${1}_R2_001.unpaired.fq" > "${processed}/${2}/${1}/trimmed/${1}.single.fq"
		gzip -c "${processed}/${2}/${1}/trimmed/${1}_R1_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz"
		gzip -c "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
	fi
fi
if [ ! -f "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" ]; then
	if [ -f "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq.gz" ]; then
		#echo "4"
		cp "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
	elif [ -f "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" ]; then
		#echo "5"
		gzip -c "${processed}/${2}/${1}/trimmed/${1}_R2_001.paired.fq" > "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz"
	fi
fi

if [ ! -d "${processed}/${2}/${1}/MLST/srst2" ]; then
	mkdir -p "${processed}/${2}/${1}/MLST/srst2"
fi

cd "${processed}/${2}/${1}/MLST/srst2"


#python2 ${shareScript}/srst2-master/scripts/getmlst.py --species "${genus} ${species}" > "${processed}/${2}/${1}/MLST/srst2/getmlst.out"
getmlst.py --species "${genus} ${species}" > "${processed}/${2}/${1}/MLST/srst2/getmlst.out"

db_name="Pasteur"
# Checks for either of the 2 databases that have multiple scheme types and runs both
if [[ "${genus}" == "Acinetobacter" ]]; then
	echo "${processed}/${2}/${1}/MLST/srst2/${genus}_${species}.fasta"
	if [[ "${species}" == "baumannii#1" ]]; then
		sed -i -e 's/Oxf_//g' "${processed}/${2}/${1}/MLST/srst2/${genus}_${species}.fasta"
		sed -i -e 's/Oxf_//g' "${processed}/${2}/${1}/MLST/srst2/abaumannii.txt"
		db_name="Oxford"
	elif [[ "${species}" == "baumannii#2" ]]; then
		sed -i -e 's/Pas_//g' "${processed}/${2}/${1}/MLST/srst2/${genus}_${species}.fasta"
		sed -i -e 's/Pas_//g' "${processed}/${2}/${1}/MLST/srst2/abaumannii_2.txt"
		db_name="Pasteur"
	else
		echo "Unknown species in Acinetobacter MLST lookup"
	fi
elif [[ "${genus}" == "Escherichia" ]]; then
	echo "${processed}/${2}/${1}/MLST/srst2/${genus}_${species}.fasta"
	if [[ "${species}" == "coli#1" ]]; then
		db_name="Achtman"
	elif [[ "${species}" == "coli#2" ]]; then
		db_name="Pasteur"
	else
		echo "Unknown species in Escherichia MLST lookup"
	fi
fi


# Pulls suggested command info from the getmlst script
suggested_command=$(tail -n2 "${processed}/${2}/${1}/MLST/srst2/getmlst.out" | head -n1)
mlst_db=$(echo "${suggested_command}" | cut -d' ' -f11)
mlst_defs=$(echo "${suggested_command}" | cut -d' ' -f13)
mlst_delimiter=$(echo "${suggested_command}" | cut -d' ' -f15)
#echo "${mlst_db}"
#echo "${mlst_defs}"
#echo "${mlst_delimiter}"

if [[ "${mlst_delimiter}" != "'_'" ]]; then
	echo "Unknown delimiter - \"${mlst_delimiter}\""
	exit
else
	mlst_delimiter="_"
	#echo "Delimiter is OK (${mlst_delimiter})"
fi

# Print out what command will be submitted
echo "--input_pe ${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz ${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz --output ${processed}/${2}/${1}/MLST/srst2 --mlst_db ${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}"
# Run the srst2 command to find MLST types
#python2 ${shareScript}/srst2-master/scripts/srst2.py --input_pe "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" --output "${processed}/${2}/${1}/srst2/${1}_ResGANNCBI" --gene_db "${ResGANNCBI_srst2}"
srst2 --input_pe "${processed}/${2}/${1}/srst2/${1}_S1_L001_R1_001.fastq.gz" "${processed}/${2}/${1}/srst2/${1}_S1_L001_R2_001.fastq.gz" --output "${processed}/${2}/${1}/MLST/srst2/${1}" --mlst_db "${mlst_db}" --mlst_definitions "${mlst_defs}" --mlst_delimiter "${mlst_delimiter}"

today=$(date "+%Y-%m-%d")

# Cleans up extra files and renames output file
mv "${processed}/${2}/${1}/MLST/srst2/${1}__mlst__${genus}_${species}__results.txt" "${processed}/${2}/${1}/MLST/${1}_srst2_${genus}_${species}-${db_name}.mlst"
mv "${processed}/${2}/${1}/MLST/srst2/mlst_data_download_${genus}_${species}_${today}.log" "${processed}/${2}/${1}/MLST/"
rm -r "${processed}/${2}/${1}/MLST/srst2"

if [[ -f "${processed}/${2}/${1}/MLST/srst2/${1}__${1}.${genus}_${species}.pileup" ]]; then
	rm -r "${processed}/${2}/${1}/MLST/srst2/${1}__${1}.${genus}_${species}.pileup"
fi
if [[ -f "${processed}/${2}/${1}/MLST/${1}__${1}.${genus}_${species}.sorted.bam" ]]; then
	rm -r "${processed}/${2}/${1}/MLST/${1}__${1}.${genus}_${species}.sorted.bam"
fi

ml -srst2 -bowtie2/2.2.4
