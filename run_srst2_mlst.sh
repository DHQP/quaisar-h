#!/bin/sh -l

#$ -o srst2.out
#$ -e srst2.err
#$ -N srst2
#$ -cwd
#$ -q short.q

#
# Description: Script to use srst2 to attempt to find mlst profile on reads. Used if standard mlst profiling fails
#
# Usage: ./run_srst2_mlst.sh -n sample_name -p MiSeq_Run_ID	-g Genus -s species [-c path_to_config_file]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/MLST/
#
# Modules required: srst2/0.2.0 bowtie2/2.2.4(?) Python2/2.7.13
#
# v1.0.3 (9/02/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml srst2 bowtie2/2.2.4

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_srst2_mlst.sh -n sample_name -p MiSeq_Run_ID	-g Genus -s species [-c path_to_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?c:p:n:g:s:" option; do
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
		g)
			echo "Option -g triggered, argument = ${OPTARG}"
			genus=${OPTARG,,}
			genus=${genus^};;
		s)
			echo "Option -s triggered, argument = ${OPTARG}"
			species=${OPTARG,,};;
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
	echo "No Project/Run_ID supplied to run_srst2_mlst.sh, exiting"
	exit 33
elif [[ -z "${sample_name}" ]]; then
	echo "No sample name supplied to run_srst2_mlst.sh, exiting"
	exit 34
fi


OUTDATADIR="${processed}/${project}/${sample_name}"

# Create output directory
if [[ ! -d "${OUTDATADIR}/srst2" ]]; then
	echo "Creating ${OUTDATADIR}/srst2"
	mkdir "${OUTDATADIR}/srst2"
fi

# Preps reads if there are no current trimmed reads in the folder
if [ ! -f "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz" ]; then
	if [ -f "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]; then
		#echo "1"
		cp "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq.gz" "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz"
	elif [ -f "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" ]; then
		#echo "2"
		gzip -c "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" > "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz"
	else
		#echo "3"
		if [[ -f "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq.gz" ]] && [[ ! -f "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq" ]]; then
			gunzip -c "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq.gz" > "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq"
		fi
		if [[ -f "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq.gz" ]] && [[ ! -f "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq" ]]; then
			gunzip -c "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq.gz" > "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq"
		fi
		echo "Running BBDUK and trimmomatic"
		ml BBMAP/38.26 trimmomatic/0.35
		bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq" in2="${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq" out="${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R1.fsq" out2="${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
		trimmomatic "${trim_endtype}" -"${trim_phred}" -threads "${procs}" "${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R1.fsq" "${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R2.fsq" "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R1_001.unpaired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R2_001.unpaired.fq" ILLUMINACLIP:"${trim_adapter_location}:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome}" SLIDINGWINDOW:"${trim_window_size}:${trim_window_qual}" LEADING:"${trim_leading}" TRAILING:"${trim_trailing}" MINLEN:"${trim_min_length}"
		ml -BBMAP/38.26 -trimmomatic/0.35
		#cat "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" > "${OUTDATADIR}/trimmed/${sample_name}.paired.fq"
		#cat "${OUTDATADIR}/trimmed/${sample_name}_R1_001.unpaired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R2_001.unpaired.fq" > "${OUTDATADIR}/trimmed/${sample_name}.single.fq"
		gzip -c "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" > "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz"
		gzip -c "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" > "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz"
	fi
fi
if [ ! -f "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz" ]; then
	if [ -f "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]; then
		#echo "4"
		cp "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq.gz" "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz"
	elif [ -f "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" ]; then
		#echo "5"
		gzip -c "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" > "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz"
	fi
fi

if [ ! -d "${OUTDATADIR}/MLST/srst2" ]; then
	mkdir -p "${OUTDATADIR}/MLST/srst2"
fi

cd "${OUTDATADIR}/MLST/srst2"


#python2 ${shareScript}/srst2-master/scripts/getmlst.py --species "${genus} ${species}" > "${OUTDATADIR}/MLST/srst2/getmlst.out"
getmlst.py --species "${genus} ${species}" > "${OUTDATADIR}/MLST/srst2/getmlst.out"

db_name="Pasteur"
# Checks for either of the 2 databases that have multiple scheme types and runs both
if [[ "${genus}" == "Acinetobacter" ]]; then
	echo "${OUTDATADIR}/MLST/srst2/${genus}_${species}.fasta"
	if [[ "${species}" == "baumannii#1" ]]; then
		sed -i -e 's/Oxf_//g' "${OUTDATADIR}/MLST/srst2/${genus}_${species}.fasta"
		sed -i -e 's/Oxf_//g' "${OUTDATADIR}/MLST/srst2/abaumannii.txt"
		db_name="Oxford"
	elif [[ "${species}" == "baumannii#2" ]]; then
		sed -i -e 's/Pas_//g' "${OUTDATADIR}/MLST/srst2/${genus}_${species}.fasta"
		sed -i -e 's/Pas_//g' "${OUTDATADIR}/MLST/srst2/abaumannii_2.txt"
		db_name="Pasteur"
	else
		echo "Unknown species in Acinetobacter MLST lookup"
	fi
elif [[ "${genus}" == "Escherichia" ]]; then
	echo "${OUTDATADIR}/MLST/srst2/${genus}_${species}.fasta"
	if [[ "${species}" == "coli#1" ]]; then
		db_name="Achtman"
	elif [[ "${species}" == "coli#2" ]]; then
		db_name="Pasteur"
	else
		echo "Unknown species in Escherichia MLST lookup"
	fi
fi


# Pulls suggested command info from the getmlst script
suggested_command=$(tail -n2 "${OUTDATADIR}/MLST/srst2/getmlst.out" | head -n1)
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
echo "--input_pe ${OUTDATADIR}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz ${OUTDATADIR}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz --output ${OUTDATADIR}/MLST/srst2 --mlst_db ${mlst_db} --mlst_definitions ${mlst_defs} --mlst_delimiter ${mlst_delimiter}"
# Run the srst2 command to find MLST types
#python2 ${shareScript}/srst2-master/scripts/srst2.py --input_pe "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz" "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz" --output "${OUTDATADIR}/srst2/${sample_name}_ResGANNCBI" --gene_db "${ResGANNCBI_srst2}"
srst2 --input_pe "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R1_001.fastq.gz" "${OUTDATADIR}/srst2/${sample_name}_S1_L001_R2_001.fastq.gz" --output "${OUTDATADIR}/MLST/srst2/${sample_name}" --mlst_db "${mlst_db}" --mlst_definitions "${mlst_defs}" --mlst_delimiter "${mlst_delimiter}"

today=$(date "+%Y-%m-%d")

# Cleans up extra files and renames output file
mv "${OUTDATADIR}/MLST/srst2/${sample_name}__mlst__${genus}_${species}__results.txt" "${OUTDATADIR}/MLST/${sample_name}_srst2_${genus}_${species}-${db_name}.mlst"
mv "${OUTDATADIR}/MLST/srst2/mlst_data_download_${genus}_${species}_${today}.log" "${OUTDATADIR}/MLST/"
rm -r "${OUTDATADIR}/MLST/srst2"

if [[ -f "${OUTDATADIR}/MLST/srst2/${sample_name}__${sample_name}.${genus}_${species}.pileup" ]]; then
	rm -r "${OUTDATADIR}/MLST/srst2/${sample_name}__${sample_name}.${genus}_${species}.pileup"
fi
if [[ -f "${OUTDATADIR}/MLST/${sample_name}__${sample_name}.${genus}_${species}.sorted.bam" ]]; then
	rm -r "${OUTDATADIR}/MLST/${sample_name}__${sample_name}.${genus}_${species}.sorted.bam"
fi

ml -srst2 -bowtie2/2.2.4
