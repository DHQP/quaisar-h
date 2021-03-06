#!/bin/sh -l

#$ -o make_Seqlog_from_log.out
#$ -e make_Seqlog_from_log.err
#$ -N make_Seqlog_from_log
#$ -cwd
#$ -q short.q

#
# Description: Creates a tsv file that matches the order of samples on the seq_log (instead of matching the list created in QuAISAR)
#
# Usage: ./make_Seqlog_from_list.sh -p run_ID [-y year] [-c path_to_config_file]
#
# Output location: /deafult_config.sh_output_location/run_ID
#
# Modules required: Python3/3.5.4
#
# v1.0.9 (08/21/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml Python3/3.5.4


#  Function to print out help blurb
show_help () {
	echo "Usage is ./make_Seqlog_from_log.sh -p run_ID [-y year] [-c path_to_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?p:c:y:" option; do
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
		y)
			echo "Option -y triggered, argument = ${OPTARG}"
			year=${OPTARG};;
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

# Show help info for when no options are given
if [[ "${options_found}" -eq 0 ]]; then
	echo "No options found"
	show_help
	exit
fi

# Checks for proper argumentation
if [[ -z "${project}" ]]; then
	echo "List empty or non-existent, exiting"
	exit 1
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



if [[ -z "${year}" ]] || [[ "${year}" = "2019_2020" ]] || [[ "${year}" = "2019" ]] || [[ "${year}" = "2020" ]]; then
	year=2019_2020
elif [[ "${year}" = "2017" ]]; then
	year=2017
elif [[ "${year}" = "2018" ]]; then
	year=2018
else
	echo "Year does not fall between 2017 and 2020, retry, exiting"
	exit 43
fi

# Creates a dictionary of commonly found bugs to use when looking up sizes and assembly ratios later
declare -A mmb_bugs
while IFS= read -r bug_lines || [ -n "$bug_lines" ]; do
	bug_genus=$(echo "${bug_lines}" | cut -d'	' -f1)
	bug_species=$(echo "${bug_lines}" | cut -d'	' -f2)
	bug_info=$(echo "${bug_lines}" | cut -d'	' -f3-)
	bug_size=$(echo "${bug_lines}" | cut -d'	' -f6)
	bug_name="${bug_genus:0:1}.${bug_species}"
	#echo "Should be adding ${bug_size} for ${bug_name}"
	mmb_bugs["${bug_name}"]="${bug_size}"
done < ${local_DBs}/MMB_Bugs.txt



owd=$(pwd)

cd "${processed}/${project}/"
ls  *run_summary* | sort -n -t _  -k 10,10r -k 8,8r -k 9,9r > "sorted_summaries.txt"
cd ${owd}

rm "${processed}/${project}/sorted_summaries.txt"
rm "${processed}/${project}/${project}_list_ordered.txt"

# Order samples (according to logsheet) in folder if not already done so
python3 ${shareScript}/order_samples_3.py -i "${local_DBs}/Seqlog_copies/${year}_MMBSeq_Log.xlsx" -r ${project} -s "Miseq Isolate Log" -o "${processed}/${project}/${project}_list_ordered.txt" -l "${processed}/${project}/${project}_list.txt"
if [[ ! -s "${processed}/${project}/${project}_list_ordered.txt" ]]; then
	echo "Isolates were not able to be sorted, something wrong with MiSeq Log entries or list file, or....?"
	exit
else
	echo "sorted file contains entries"
fi

if [[ -f "${processed}/${project}/seqlog_output.txt" ]]; then
	rm -r "${processed}/${project}/seqlog_output.txt"
fi

> "${processed}/${project}/Seqlog_output.txt"

# Goes through each item on the list and pulls all relevant info
while IFS= read -r var || [ -n "$var" ]; do
	# Current (8/16/17) order of expected run output
	#  kraken - QC - estimated coverage - #contigs - cumulative length assmbly - BUSCO - ANI
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
			if [ "${first}" = "s" ]; then
				species_post=$(echo "${line}" | awk -F ' ' '{print $4}')
			elif [ "${first}" = "G" ]; then
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
	# Pulls species and genus_16s information from kraken out of assembly
	if [[ -s "${OUTDATADIR}/kraken/preAssembly/${sample_name}_kraken_summary_paired.txt" ]]; then
		while IFS= read -r line  || [ -n "$line" ]; do
			first=${line::1}
			if [ "${first}" = "s" ]; then
				species_reads=$(echo "${line}" | awk -F ' ' '{print $4}')
			elif [ "${first}" = "G" ]; then
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
			if [ ${counter} -eq 0 ]; then
				num_contigs=$(sed -n '14p' "${OUTDATADIR}/Assembly_Stats/${sample_name}_report.tsv"| sed -r 's/[\t]+/ /g' | cut -d' ' -f3 )
			elif [ ${counter} -eq 1 ]; then
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
			elif [ ${counter} -eq 3 ]; then
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
			elif [[ ${line} == *"Total BUSCO groups searched"* ]]; then
				#echo "T-"${line}
				total_buscos=$(echo "${line}" | awk -F ' ' '{print $1}')
			elif [[ "${line}" == *"The lineage dataset is:"* ]]; then
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
	ani_dec_genus=${dec_genus}
	if [[ "${ani_dec_genus}" == "Clostridioides" ]] || [[ "${ani_dec_genus}" == "Hungateiclostridium" ]]; then
		ani_dec_genus="Clostridium"
	fi
	if [[ -s "${OUTDATADIR}/ANI/best_hits_ordered.txt" ]]; then
		mv "${OUTDATADIR}/ANI/best_hits_ordered.txt" "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${ani_dec_genus}).txt"
	fi
	# If 1 and only 1 file exists pull the first line as the best hit information
	# echo "test-${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${dec_genus}).txt"
	if [[ -s "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_REFSEQ_${REFSEQ_date}).txt" ]]; then
		ani_info=$(head -n 1 "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_REFSEQ_${REFSEQ_date}).txt")
	elif [[ -s "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_All.txt" ]]; then
		ani_info=$(head -n 1 "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_All).txt")
	elif [[ -s "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${ani_dec_genus}).txt" ]]; then
		ani_info=$(head -n 1 "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${ani_dec_genus}).txt")
	# Report that more than one file exists
	else
		for file in "${OUTDATADIR}/ANI/"*
		do
			if [[ "${file}" == *"best_ANI_hits_ordered(${sample_name}_vs_"* ]]; then
				filename=${file}
				echo "${OUTDATADIR}"
				echo "${file}"
				echo "${ani_dec_genus}"
				echo "${dec_genus^}"
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
	echo -e "${sample_name}\\t${NOW}\\t${g_s_reads}\\t${g_s_assembled}\\t${g_s_16s}\\t${read_qc_info}\\t${avg_coverage}\\t${contig_info}\\t${busco_info}\\t${ani_info}\\r" >> "${processed}/${project}/Seqlog_output.txt"
done < ${processed}/${project}/${project}_list_ordered.txt

ml -Python3/3.5.4

#Script exited gracefully
exit 0
