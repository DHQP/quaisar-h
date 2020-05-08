#!/bin/sh -l

#$ -o make_Seqlog_from_list.out
#$ -e make_Seqlog_from_list.err
#$ -N Seqlog_list
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Creates a tsv file that matches the order of samples on the this when the run was created in QuAISAR (instead of matching the MMB_Seq log)
#
# Usage: ./make_Seqlog_from_list.sh path_to_list
#
# Output location: same folder as path_to_list input parameter
#
# Modules required: None
#
# v1.0.2 (05/08/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to $0, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./make_Seqlog_from_list.sh list_of_samples"
	temp_dir=$(dirname ${1})
	echo "Output is saved to ${temp_dir}/Seqlog_output.txt"
	exit 0
fi

# Creates a dictionary of commonly found bugs to use when looking up sizes and assembly ratios later
declare -A mmb_bugs
while IFS= read -r bug_lines  || [ -n "$bug_lines" ]; do
	bug_genus=$(echo "${bug_lines}" | cut -d'	' -f1)
	bug_species=$(echo "${bug_lines}" | cut -d'	' -f2)
	bug_info=$(echo "${bug_lines}" | cut -d'	' -f3-)
	bug_size=$(echo "${bug_lines}" | cut -d'	' -f6)
	bug_name="${bug_genus:0:1}.${bug_species}"
	#echo "Should be adding ${bug_size} for ${bug_name}"
	mmb_bugs["${bug_name}"]="${bug_size}"
done < ${local_DBs}/MMB_Bugs.txt

# Set output folder as directory of input list
output_folder=$(dirname ${1})
> "${output_folder}/Seqlog_output.txt"

# Goes through each item on the list and pulls all relevant info
while IFS= read -r var || [ -n "$var" ]; do
	# Current (12/17/18) order of expected run output
	#  kraken - QC - estimated coverage - #contigs - cumulative length assembly - BUSCO - ANI

	project=$(echo "${var}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
	sample_name=$(echo "${var}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
	OUTDATADIR="${processed}/${project}/${sample_name}"

	#echo "P:${project}:     S:${sample_name}:"
	#echo "O:${OUTDATADIR}:"

	# Creates default values in case they are not filled in later
	g_s_assembled="Unidentified"
	genus_post="not_assigned"
	species_post="not_assigned"
	# Pulls species and genus_post information from kraken out of assembly
	if [[ -s "${OUTDATADIR}/kraken/postAssembly/${sample_name}_kraken_summary_assembled_BP.txt" ]]; then
		while IFS= read -r line  || [ -n "$line" ]; do
			first=${line::1}
			if [ "${first}" = "S" ]
			then
				species_post=$(echo "${line}" | awk -F ' ' '{print $4}')
			elif [ "${first}" = "G" ]
			then
				genus_post=$(echo "${line}" | awk -F ' ' '{print $4}')
			fi
		done < "${OUTDATADIR}/kraken/postAssembly/${sample_name}_kraken_summary_assembled_BP.txt"
		g_s_assembled="${genus_post} ${species_post}"
		#echo "${g_s_assembly}"
	elif [[ ! -f "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]]; then
		#echo "Cant find ${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta"
		g_s_assembled="Failed_Assembly"
	fi

	g_s_reads="Unidentified"
	genus_reads="not_assigned"
	species_reads="not_assigned"
	# Pulls species and genus_16s information from kraken out of reads
	if [[ -s "${OUTDATADIR}/kraken/preAssembly/${sample_name}_kraken_summary_paired.txt" ]]; then
		while IFS= read -r line  || [ -n "$line" ]; do
			first=${line::1}
			if [ "${first}" = "S" ]
			then
				species_reads=$(echo "${line}" | awk -F ' ' '{print $4}')
			elif [ "${first}" = "G" ]
			then
				genus_reads=$(echo "${line}" | awk -F ' ' '{print $4}')
			fi
		done < "${OUTDATADIR}/kraken/preAssembly/${sample_name}_kraken_summary_paired.txt"
		g_s_reads="${genus_reads} ${species_reads}"
		#echo "${g_s}"
	fi
	# Pulls species and genus_16s information from 16s
	g_s_16s="Unidentified"
	genus_16s="not_assigned"
	species_16s="not_assigned"
	if [[ -s "${OUTDATADIR}/16s/${sample_name}_16s_blast_id.txt" ]]; then
		info=$(tail -n 1 "${OUTDATADIR}/16s/${sample_name}_16s_blast_id.txt")
		type=$(echo "${info}" | cut -d' ' -f1)
		genus_16s=$(echo "${info}" | cut -d'	' -f3 | cut -d' ' -f1)
		species_16s=$(echo "${info}" | cut -d'	' -f3 | cut -d' ' -f2)
#		echo "g-${genus_16s};s-${species}"
		if [[ "${genus_16s}" = "No_16s_sequences_found" ]] && [[ "${genus_16s}" = "No_16s_sequences_found" ]]; then
			g_s_16s="${genus_16s}"
		else
			g_s_16s="${genus_16s} ${species_16s}"
		fi
		#		echo "g_s_16-${g_s_16s}"
	elif [[ ! -f "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]]; then
		g_s_16s="Failed_Assembly"
	fi
	# Pulls QC count info from counts file (Order is as follows Q20_Total_[bp]	Q30_Total_[bp]	Q20_R1_[bp]	Q20_R2_[bp]	Q20_R1_[%]	Q20_R2_[%]	Q30_R1_[bp]	Q30_R2_[bp]
	# Q30_R1_[%]	Q30_R2_[%]	Total_Sequenced_[bp]	Total_Sequenced_[reads]
	read_qc_info="N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A	N/A"
	# If the counts file exists take the header line (the only one) and copy all but the first entry (which is the sample name) and store in an array
	if [[ -s "${OUTDATADIR}/preQCcounts/${sample_name}_counts.txt" ]]; then
		line=$(tail -n 1 "${OUTDATADIR}/preQCcounts/${sample_name}_counts.txt")
		IFS='	' read -r -a qcs <<< "${line}"
		read_qc_info=${qcs[@]:1}
		#echo "${read_qc_info}"
	fi

	source_call=$(head -n1 "${OUTDATADIR}/${sample_name}.tax").
	tax_source="UNK"
	while IFS= read -r line  || [ -n "$line" ]; do
		# Grab first letter of line (indicating taxonomic level)
		first=${line:0:1}
		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "s" ]; then
			dec_species=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "G" ]; then
			dec_genus=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "(" ]; then
			tax_source=$(echo "${line}" | cut -d')' -f1 | cut -d'(' -f2)
		fi
	done < "${OUTDATADIR}/${sample_name}.tax"



	# Pulls contig info from toms qc analysis file
	contig_info="0(0)\\t0\tNot_in_DB"
		if [[ -s "${OUTDATADIR}/Assembly_Stats/${sample_name}_report.tsv" ]]; then
		counter=0
		while IFS= read -r line  || [ -n "$line" ]; do
			if [ ${counter} -eq 0 ]
			then
				num_contigs=$(sed -n '14p' "${OUTDATADIR}/Assembly_Stats/${sample_name}_report.tsv"| sed -r 's/[\t]+/ /g' | cut -d' ' -f3 )
			elif [ ${counter} -eq 1 ]
			then
				assembly_length=$(sed -n '16p' "${OUTDATADIR}/Assembly_Stats/${sample_name}_report.tsv" | sed -r 's/[\t]+/ /g' | cut -d' ' -f3)
				#Check Assembly ratio against expected size to see if it is missing a large portion or if there is contamination/double genome
				dec_genus_initial="${dec_genus:0:1}"
				if [[ "${dec_genus_initial}" = "[" ]] || [[ "${dec_genus_initial}" = "(" ]]; then
					dec_genus_initial="${dec_genus:1:1}"
				fi
				assembly_ID="${dec_genus_initial}.${dec_species}"
				echo "About to check mmb_bugs[${assembly_ID}]"
				if [[ ! -z "${mmb_bugs[${assembly_ID}]}" ]]; then
					assembly_ratio=$(awk -v p="${assembly_length}" -v q="${mmb_bugs[${assembly_ID}]}" -v r="${tax_source}" -v s="${assembly_ID}" 'BEGIN{printf("%.2f(%s-%s)",p/q, r, s)}')
				else
					assembly_ratio="Not_in_DB (${tax_source}-${assembly_ID})"
				fi
			elif [ ${counter} -eq 3 ]
			then
				N50=$(sed -n '18p' "${OUTDATADIR}/Assembly_Stats/${sample_name}_report.tsv"  | sed -r 's/[\t]+/ /g'| cut -d' ' -f2)
			fi
			counter=$((counter+1))
		done < ${OUTDATADIR}/Assembly_Stats/${sample_name}_report.tsv
		contig_info=$(echo -e "${num_contigs}\\t${assembly_length}\\t${assembly_ratio}")
		#with N50 size
		#contig_info=$(echo -e "${num_contigs}(${n50_length})\\t${assembly_length}")
	fi

	# Pulls busco info from summary file
	busco_info="No BUSCO performed"
	if [[ -s "${OUTDATADIR}/BUSCO/short_summary_${sample_name}.txt" ]]; then
		while IFS= read -r line  || [ -n "$line" ]; do
			if [[ ${line} == *"Complete BUSCOs (C)"* ]]
			then
				#echo "C-"${line}
				found_buscos=$(echo "${line}" | awk -F ' ' '{print $1}')
			elif [[ ${line} == *"Total BUSCO groups searched"* ]];
			then
				#echo "T-"${line}
				total_buscos=$(echo "${line}" | awk -F ' ' '{print $1}')
			elif [[ "${line}" == *"The lineage dataset is:"* ]];
			then
				#echo "L-"${line}
				db=$(echo "${line}" | awk -F ' ' '{print $6}')
			fi
		done < ${OUTDATADIR}/BUSCO/short_summary_${sample_name}.txt
		busco_info="${found_buscos}/${total_buscos}(${db})"
	fi
	# Pulls ANI info from best_ANI_hits file
	ani_info="No ANI performed"
	# Count the number of matching format files for the current sample
	file_count=$(find "${OUTDATADIR}/ANI/" -name *"${sample_name}"*"_vs_"*".txt" | wc -l)
	# Rename files in old formating convention
	if [[ -s "${OUTDATADIR}/ANI/best_hits_ordered.txt" ]]; then
		mv "${OUTDATADIR}/ANI/best_hits_ordered.txt" "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${genus}).txt"
	fi
	# If 1 and only 1 file exists pull the first line as the best hit information
	# echo "test-${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${dec_genus}).txt"
	if [[ -s "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_All.txt" ]]; then
		ani_info=$(head -n 1 "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_All).txt")
	elif [[ -s "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${dec_genus}).txt" ]]; then
		ani_info=$(head -n 1 "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${dec_genus}).txt")
	# Report that more than one file exists
	else
		for file in "${OUTDATADIR}/ANI/"*
		do
			if [[ "${file}" == *"best_ANI_hits_ordered(${sample_name}_vs_"* ]]; then
				filename=${file}
				#echo "${OUTDATADIR}"
				#echo "${file}"
				#echo "${dec_genus^}"
				ani_info=$(head -n 1 "${file}")
				break
			fi
		done
	fi
	# Extract q30 reads from qcCounts to calculate average coverage as q30_reads/assembly_length
	q30_reads=$(echo "${read_qc_info}" | awk -F ' ' '{print $2}')
	# Change later to AWK as this wont work on ASPEN, but consolidate won't likely be run on cluster
	if [[ ${assembly_length} -gt 0 ]]; then
		avg_coverage=$(bc <<<"scale=2 ; ${q30_reads} / ${assembly_length}")
	else
		avg_coverage="N/A"
	fi
	# Replace all spaces in qc_info as tabs
	read_qc_info=$(echo "${read_qc_info}" | tr '[:blank:]' \\t)

	# Get the date to show when the log was made
	NOW=$(date +"%m/%d/%Y")

	# Add all pertinent info to the output file in the correct formatting to add to MMB_Seq log
	echo -e "${sample_name}\\t${NOW}\\t${g_s_reads}\\t${g_s_assembled}\\t${g_s_16s}\\t${read_qc_info}\\t${avg_coverage}\\t${contig_info}\\t${busco_info}\\t${ani_info}\\r" >> "${output_folder}/Seqlog_output.txt"
 done < ${1}

#Script exited gracefully
exit 0
