#!/bin/sh -l

#$ -o run_BBDUK.out
#$ -e run_BBDUK.err
#$ -N run_BBDUK
#$ -cwd
#$ -q short.q

# Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Runs BBDUK on sample to align reads into best possible assembly
#
# Usage ./run_BBDUK.sh sample_name run_ID
#
# Modules required BBDUK/3.10.1
#

ml BBMap/38.26 java/latest

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_BBDUK.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_BBDUK.sh sample_name run_ID"
	echo "Output by default is sent to ${processed}/miseq_run_ID/sample_name/removedAdapters"
	exit 0
elif [ -z "${2}" ]; then
	echo "Empty project id supplied to run_BBDUK.sh, exiting"
	exit 1
fi

# Sets the output folder to the sample_name folder in processed samples
sample="${1}"
project="${2}"
OUTDATADIR="${processed}/${project}/${sample}"

if [[ -f "${OUTDATADIR}/FASTQs/${sample}_R1_001.fastq" ]]; then
	file_in_1="${OUTDATADIR}/FASTQs/${sample}_R1_001.fastq"
elif [[ -f "${OUTDATADIR}/FASTQs/${sample}_R1_001.fastq.gz" ]]; then
	file_in_1="${OUTDATADIR}/FASTQs/${sample}_R1_001.fastq.gz"
else
	echo "No input R1 fastq found (${OUTDATADIR}/FASTQs/${sample}_R1_001.fastq(.gz)"
	exit
fi

if [[ -f "${OUTDATADIR}/FASTQs/${sample}_R2_001.fastq" ]]; then
	file_in_2="${OUTDATADIR}/FASTQs/${sample}_R2_001.fastq"
elif [[ -f "${OUTDATADIR}/FASTQs/${sample}_R2_001.fastq.gz" ]]; then
	file_in_2="${OUTDATADIR}/FASTQs/${sample}_R2_001.fastq.gz"
else
	echo "No input R2 fastq found (${OUTDATADIR}/FASTQs/${sample}_R2_001.fastq(.gz)"
	exit
fi




bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${file_in_1}" in2="${file_in_2}" out="${OUTDATADIR}/removedAdapters/${sample}-noPhiX-R1.fsq" out2="${OUTDATADIR}/removedAdapters/${sample}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
remAdapt_length_R1=$(cat ${OUTDATADIR}/removedAdapters/${sample}-noPhiX-R1.fsq | paste - - - - | cut -f2 |tr -d '\n' | wc -c)
remAdapt_length_R2=$(cat ${OUTDATADIR}/removedAdapters/${sample}-noPhiX-R2.fsq | paste - - - - | cut -f2 |tr -d '\n' | wc -c)
echo -e "R1:	${remAdapt_length_R1}\nR2:	${remAdapt_length_R2}" > "${OUTDATADIR}/removedAdapters/no_PhiX_total_lengths.txt"

ml -BBMap/38.26

#Script exited gracefully (unless something else inside failed)
exit 0
