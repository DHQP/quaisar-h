#!/bin/sh -l

#$ -o run_plasFlow.out
#$ -e run_plasFlow.err
#$ -N run_plasFlow
#$ -cwd
#$ -q short.q

#
# Description: Will attempt to find any plasmids in sample using plasFlow methods
#
# Usage: ./run_plasFlow.sh -n sample_name -p run_ID [-c path_to_config_file]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/plasFlow/
#
# Modules required: PlasFlow/1.1
#
# v1.0.2 (09/08/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_plasFlow.sh -n sample_name -p run_ID [-c path_to_config_file]"
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

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${processed}/${project}/${sample_name}"



ml PlasFlow/1.1 Python3/3.5.4

# Create output directory
if [[ ! -d "${OUTDATADIR}/plasFlow" ]]; then
	mkdir "${OUTDATADIR}/plasFlow"
fi

# Prep any samples that don't have the paired.fq reads
if [ -f "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" ] && [ -f "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" ]; then
	echo "1"
	#gzip -c "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" > "${OUTDATADIR}/trimmed/${sample_name}_S1_L001_R1_001.fastq.gz"
	#gzip -c "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" > "${OUTDATADIR}/trimmed/${sample_name}_S1_L001_R2_001.fastq.gz"
elif [ -f "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]; then
	echo "2"
	gunzip < "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq.gz" > "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq"
	if [[ -f "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]]; then
	echo "2A"
		gunzip < "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq.gz" > "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq"
	fi
elif [ -f "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]; then
	echo "3"
	gunzip < "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq.gz" > "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq"
else
	echo "4"
	if [[ -f "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq.gz" ]] && [[ ! -f "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq" ]]; then
		gunzip -c "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq.gz" > "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq"
	fi
	if [[ -f "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq.gz" ]] && [[ ! -f "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq" ]]; then
		gunzip -c "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq.gz" > "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq"
	fi
	echo "Running BBDUK and trimmomatic"
	bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq" in2="${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fastq" out="${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R1.fsq" out2="${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
	trimmomatic "${trim_endtype}" -"${trim_phred}" -threads "${procs}" "${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R1.fsq" "${OUTDATADIR}/removedAdapters/${sample_name}-noPhiX-R2.fsq" "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R1_001.unpaired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R2_001.unpaired.fq" ILLUMINACLIP:"${trim_adapter_location}:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome}" SLIDINGWINDOW:"${trim_window_size}:${trim_window_qual}" LEADING:"${trim_leading}" TRAILING:"${trim_trailing}" MINLEN:"${trim_min_length}"
fi

# Check if sample has original assembly to process plasflow from
if [[ -s "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]]; then
	# Trim contigs a little to 2000 and larger and put through plasflow.
	python3 "${shareScript}/removeShortContigs.py" -i "${OUTDATADIR}/Assembly/scaffolds.fasta" -t 2000 -s "normal_SPAdes"
	# Renames headers of fasta files accordingly
	python3 "${shareScript}/fasta_headers.py" -i "${OUTDATADIR}/Assembly/scaffolds.fasta.TRIMMED.fasta" -o "${OUTDATADIR}/plasFlow/${sample_name}_scaffolds_trimmed_2000.fasta"
	# Removes intermeidate fasta file
	rm -r "${OUTDATADIR}/Assembly/scaffolds.fasta.TRIMMED.fasta"
	# Run plasflow on newly trimmed assembly file
	PlasFlow.py --input "${OUTDATADIR}/plasFlow/${sample_name}_scaffolds_trimmed_2000.fasta" --output "${OUTDATADIR}/plasFlow/${sample_name}_plasFlow_results.tsv" --threshold 0.7
	# Load all necessary modules to complete the realignment portion of analysis

	ml -Python3/3.5.4 bowtie2/2.2.9 samtools/1.4.1 bam2fastq/1.1.0 Unicycler/0.4.4 #SPAdes/3.13.0 racon/1.3.1

	mkdir ${OUTDATADIR}/plasFlow/bowtie2-index/
	bowtie2-build -f "${OUTDATADIR}/plasFlow/${sample_name}_plasFlow_results.tsv_chromosomes.fasta" "${OUTDATADIR}/plasFlow/bowtie2-index/bowtie2_${sample_name}_chr"
	mkdir ${OUTDATADIR}/plasFlow/filtered_reads_70/
	bowtie2 -x "${OUTDATADIR}/plasFlow/bowtie2-index/bowtie2_${sample_name}_chr" -1 "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" -2 "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" -S "${OUTDATADIR}/plasFlow/filtered_reads_70/${sample_name}.sam" -p 12 --local
	samtools view -bS "${OUTDATADIR}/plasFlow/filtered_reads_70/${sample_name}.sam" > "${OUTDATADIR}/plasFlow/filtered_reads_70/${sample_name}.bam"
	bam2fastq --no-aligned -o "${OUTDATADIR}/plasFlow/filtered_reads_70/${sample_name}_R#_bacterial.fastq" "${OUTDATADIR}/plasFlow/filtered_reads_70/${sample_name}.bam"
	unicycler -1 "${OUTDATADIR}/plasFlow/filtered_reads_70/${sample_name}_R_1_bacterial.fastq" -2 "${OUTDATADIR}/plasFlow/filtered_reads_70/${sample_name}_R_2_bacterial.fastq" -o "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly"
	mv "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/assembly.fasta" "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_original.fasta"
	mv "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/assembly.gfa" "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_assembly.gfa"
	python3 ${shareScript}/fasta_headers_plasFlow.py -i "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_original.fasta" -o "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly.fasta"
	python3 "${shareScript}/removeShortContigs.py" -i "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly.fasta" -t 500 -s "plasFlow"
	mv "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly.fasta.TRIMMED.fasta" "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta"
	rm "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly.fasta"
else
	echo "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta (Assembly) not found, cant do anything"
fi

ml -Python3/3.5.4 -bowtie2/2.2.9 -samtools/1.4.1 -bam2fastq/1.1.0 -Unicycler/0.4.4 -SPAdes/3.13.0 -racon/1.3.1
