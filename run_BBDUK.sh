#!/bin/sh -l

#$ -o run_BBDUK.out
#$ -e run_BBDUK.err
#$ -N run_BBDUK
#$ -cwd
#$ -q short.q

#
# Runs BBDUK on sample to align reads into best possible assembly
#
# Usage ./run_BBDUK.sh -n sample_name -p run_ID [-c path_to_config_file]
#
# Modules required BBDUK/3.10.1
#

ml BBMap/38.26 java/latest

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_Assembly_Quality_Check.sh -n sample_name -p run_ID [-c path_to_config_file]"
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

if [[ -z "${project}" ]]; then
	echo "No Project/Run_ID supplied to retry_ANI_best_hit.sh, exiting"
	exit 33
elif [[ -z "${sample_name}" ]]; then
	echo "No sample name supplied to retry_ANI_best_hit.sh, exiting"
	exit 34
fi


OUTDATADIR="${processed}/${project}/${sample_name}"

if [[ -f "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq" ]]; then
	file_in_1="${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq"
elif [[ -f "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq.gz" ]]; then
	file_in_1="${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq.gz"
else
	echo "No input R1 fastq found (${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq(.gz)"
	exit
fi

if [[ -f "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq" ]]; then
	file_in_2="${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq"
elif [[ -f "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq.gz" ]]; then
	file_in_2="${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq.gz"
else
	echo "No input R2 fastq found (${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq(.gz)"
	exit
fi




bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${file_in_1}" in2="${file_in_2}" out="${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R1.fsq" out2="${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
remAdapt_length_R1=$(cat ${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R1.fsq | paste - - - - | cut -f2 |tr -d '\n' | wc -c)
remAdapt_length_R2=$(cat ${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R2.fsq | paste - - - - | cut -f2 |tr -d '\n' | wc -c)
echo -e "R1:	${remAdapt_length_R1}\nR2:	${remAdapt_length_R2}" > "${OUTDATADIR}/removedAdapters/no_PhiX_total_lengths.txt"

ml -BBMap/38.26

#Script exited gracefully (unless something else inside failed)
exit 0
