#!/bin/sh -l

#$ -o 16s_blast.out
#$ -e 16s_blast.err
#$ -N 16s_blast
#$ -cwd
#$ -q short.q

#
# Description: Creates a species prediction based on blasting the largest and also 'best' hit of the suggested 16s sequences found using barrnap
# NB: Hard coding thread ${procs} count usage as 4, until we put in a way to guage available threads
#
# Usage: ./16s_blast.sh -i path_to_sample_folder
#
# Output location: path_to_sample_folder/sample_name/16s/ path_to_scripts
#
# Modules required: barrnap/0.8, ncbi-blast+/LATEST, perl/5.12.3. hmmer/3.1b2 loaded by barrnap/0.8
#
# v1.0.3 (05/14/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml perl/5.12.3 barrnap/0.8 ncbi-blast+/LATEST Python3/3.5.2

#  Function to print out help blurb
show_help () {
	echo "Usage is ./16s_blast.sh -i path_to_sample_folder"
	echo "Output is saved to path_to_sample_folder/run_ID/sample_name/16s"
}

# Parse command line options
options_found=0
while getopts ":h?n:p:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		i)
			echo "Option -i triggered, argument = ${OPTARG}"
			SAMPDATADIR=${OPTARG};;
		s)
			echo "Option -s triggered, argument = ${OPTARG}"
			scripts=${OPTARG};;
		:)
			echo "Option -${OPTARG} requires as argument";;
		h)
			show_help
			exit 0
			;;
	esac
done

# Show help info for when no options are given
if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
fi

# Sets hard-coded processors variable, until a better way is found
procs=4

# Creates new output folder based on the universally set processed location from config.sh
if [[ ! -d "${SAMPDATADIR}" ]]; then
	echo "Sample folder (${SAMPDATADIR}) does not exist, exiting"
	exit 1
elif [[ ! -d "${scripts}" ]]; then
	echo "Script folder (${scripts}) does not exist, exiting"
	exit 2
fi

sample_name=$(echo ${SAMPDATADIR} | rev | cut -d'/' -f1 | rev)

if [ ! -d "${SAMPDATADIR}/16s" ]; then
	echo "Creating $SAMPDATADIR/16s"
	mkdir "${SAMPDATADIR}/16s"
fi

# Get original working directory so that it can return to it after running (may not need to do this but havent tested it out yet)
owd=$(pwd)
cd ${SAMPDATADIR}/16s

# Run barrnap to discover ribosomal sequences
barrnap --kingdom bac --threads ${procs} "${SAMPDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta" > ${SAMPDATADIR}/16s/${sample_name}_rRNA_finds.txt

# Checks for successful output from barrnap, *rRNA_seqs.fasta
if [[ ! -s ${SAMPDATADIR}/16s/${sample_name}_rRNA_finds.txt ]]; then
	echo "rNA_seqs.fasta does NOT exist"
	exit 1
fi

# Checks barrnap output and finds all 16s hits and creates a multi-fasta file to list all possible matches
lines=-1
found_16s="false"
if [[ -f "${SAMPDATADIR}/16s/${sample_name}_rRNA_finds.txt" ]]; then
	while IFS='' read -r line; do
		contig=$(echo ${line} | cut -d' ' -f1)
		cstart=$(echo ${line} | cut -d' ' -f4)
		cstop=$(echo ${line} | cut -d' ' -f5)
		ribosome=$(echo ${line} | cut -d' ' -f9 | cut -d'=' -f3)
		if [ "${ribosome}" = "16S" ]; then
			# Replace with subsequence once it can handle multi-fastas
			#make_fasta $1 $2 $contig $cstart $cstop
			lines=$((lines + 1))
			echo "About to make ${sample_name}_16s_rna_seq_${lines}.fasta"
			python3 ${scripts}/get_subsequence.py -i "${SAMPDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta" -s ${cstart} -e ${cstop} -t ${contig} -o "${contig} 16s-${lines}" >> ${SAMPDATADIR}/16s/${sample_name}_16s_rna_seq_${lines}.fasta
			found_16s="true"
		fi
	done < "${SAMPDATADIR}/16s/${sample_name}_rRNA_finds.txt"
else
	echo "Barrnap Failed maybe??!??"
fi

# Adds No hits found to output file in the case where no 16s ribosomal sequences were found
if [[ "${found_16s}" == "false" ]]; then
	echo -e "best_hit	${sample_name}	No_16s_sequences_found" > "${SAMPDATADIR}/16s/${sample_name}_16s_blast_id.txt"
	echo -e "largest_hit	${sample_name}	No_16s_sequences_found" >> "${SAMPDATADIR}/16s/${sample_name}_16s_blast_id.txt"
	exit
fi

# Blasts the NCBI database to find the closest hit to every entry in the 16s fasta list
###### MAX_TARGET_SEQS POSSIBLE ERROR
line_inc=0
while [[ -f ${sample_name}_16s_rna_seq_${line_inc}.fasta ]]; do
	echo "Blasting ${sample_name}_16s_rna_seq_${line_inc}.fasta"
	blastn -word_size 10 -task blastn -remote -db nt -max_hsps 1 -max_target_seqs 10 -query ${processed}/${project}/${sample_name}/16s/${sample_name}_16s_rna_seq_${line_inc}.fasta -out ${SAMPDATADIR}/16s/${sample_name}.nt.RemoteBLASTN_${line_inc} -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen ssciname";
	line_inc=$(( line_inc + 1 ))
done

cat ${SAMPDATADIR}/16s/${sample_name}.nt.RemoteBLASTN_* > ${SAMPDATADIR}/16s/${sample_name}.nt.RemoteBLASTN_all

# Sorts the list based on sequence match length to find the largest hit
sort -k4 -n "${SAMPDATADIR}/16s/${sample_name}.nt.RemoteBLASTN_all" --reverse > "${SAMPDATADIR}/16s/${sample_name}.nt.RemoteBLASTN_all.sorted"

# Gets taxon info from the best bitscore (literal top) hit from the blast list
if [[ -s "${SAMPDATADIR}/16s/${sample_name}.nt.RemoteBLASTN_all" ]]; then
	me=$(whoami)
	accessions=$(head -n 1 "${SAMPDATADIR}/16s/${sample_name}.nt.RemoteBLASTN_all")
	hits=$(echo "${SAMPDATADIR}/16s/${sample_name}.nt.RemoteBLASTN_all" | wc -l)
#	echo ${accessions}
	gb_acc=$(echo "${accessions}" | cut -d' ' -f2 | cut -d'|' -f4)
	echo ${gb_acc}
	attempts=0
	# Will try getting info from entrez up to 5 times, as it has a higher chance of not finishing correctly on the first try
	while [[ ${attempts} -lt 5 ]]; do
		blast_id=$(python ${scripts}/entrez_get_taxon_from_accession.py -a "${gb_acc}" -e "${me}@cdc.gov")
		if [[ ! -z ${blast_id} ]]; then
			break
		else
			attempts=$(( attempts + 1 ))
		fi
		sleep 1
	done
	echo ${blast_id}
	if [[ -z ${blast_id} ]]; then
		blast_id="No_16s_matches_found"
	fi
	#blast_id=$(echo ${blast_id} | tr -d '\n')
	echo -e "best_hit	${sample_name}	${blast_id}" > "${SAMPDATADIR}/16s/${sample_name}_16s_blast_id.txt"
	if [[ "${hits}" -eq 1 ]]; then
		echo -e "largest	${sample_name}	${blast_id}" >> "${SAMPDATADIR}/16s/${sample_name}_16s_blast_id.txt"
		skip_largest="true"
	fi
else
	echo "No remote blast file"
fi

best_blast_id=${blast_id}

# Gets taxon info from the largest hit from the blast list
if [[ ${skip_largest} != "true" ]]; then
	# Gets taxon info from the largest hit from the blast list
	if [[ -s "${SAMPDATADIR}/16s/${sample_name}.nt.RemoteBLASTN_all.sorted" ]]; then
		me=$(whoami)
		accessions=$(head -n 1 "${SAMPDATADIR}/16s/${sample_name}.nt.RemoteBLASTN_all.sorted")
		gb_acc=$(echo "${accessions}" | cut -d' ' -f2 | cut -d'|' -f4)
		attempts=0
		# Will try getting info from entrez up to 5 times, as it has a higher chance of not finishing correctly on the first try
		while [[ ${attempts} -lt 5 ]]; do
			blast_id=$(python ${scripts}/entrez_get_taxon_from_accession.py -a "${gb_acc}" -e "${me}@cdc.gov")
			if [[ ! -z ${blast_id} ]]; then
				break
			else
				attempts=$(( attempts + 1 ))
			fi
		done
		echo ${blast_id}
		if [[ -z ${blast_id} ]]; then
			blast_id="No_16s_matches_found"
		fi
		#	blast_id$(echo ${blast_id} | tr -d '\n')
		if [[ "${hits}" -eq 1 ]] && [[ "${best_blast_id}" == "No_16s_matches_found" ]]; then
			echo -e "best_hit	${sample_name}	${blast_id}" > "${SAMPDATADIR}/16s/${sample_name}_16s_blast_id.txt"
		fi
		echo -e "largest	${sample_name}	${blast_id}" >> "${SAMPDATADIR}/16s/${sample_name}_16s_blast_id.txt"
	else
		echo "No sorted remote blast file"
	fi
fi

# Go back to original working directory
cd ${owd}

# Return modules to original versions
ml -perl/5.12.3 -barrnap/0.8 -ncbi-blast+/LATEST

#Script exited gracefully (unless something else inside failed)
exit 0
