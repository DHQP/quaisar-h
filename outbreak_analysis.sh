#!/bin/sh -l

#$ -o outbreak_analysis.out
#$ -e outbreak_analysis.err
#$ -N outbreak_analysis
#$ -cwd
#$ -q short.q

#
# Description: Pulls out MLST, AR genes, and plasmid repicons and creates a mashtree for the listed samples and consolidates them into one sheet when run from an alternate or old database
#
# Usage ./outbreak_analysis.sh -l path_to_list -g gapped/ungapped (analysis ran) -s identity (80/95/98/99/100) -t analysis_type(MATRIX|SNV|BOTH) -n analysis_identifier(e.g. outbreak identifier) -k clobberness[keep|clobber] [-c path_to_config_file] [-d path_to_alt_DB]
#
# Output location: Parameter
#
# Modules required: Python3/3.5.2, mashtree/0.29
#		***Must be submitted as a job (or run on the cluster) if there are isolates that need to have csstar, GAMA or srst2 updated
#
# v1.0.4 (08/24/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml Python3/3.5.2 mashtree/0.29


#  Function to print out help blurb
show_help () {
	echo "./outbreak_analysis.sh -l path_to_list -g gapped/ungapped (analysis ran) -s identity (80/95/98/99/100) -t analysis_type (MATRIX|SNV|BOTH)-n analysis_identifier(e.g. outbreak identifier) -k clobberness[keep|clobber] [-c path_to_config_file] [-d path_to_alt_DB]"
}

# Parse command line options
options_found=0
while getopts ":h?l:c:g:n:s:k:d:t:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		c)
			echo "Option -c triggered, argument = ${OPTARG}"
			config=${OPTARG};;
		g)
			echo "Option -g triggered, argument = ${OPTARG}"
			gapping=${OPTARG};;
		n)
			echo "Option -n triggered, argument = ${OPTARG}"
			analysis_name=${OPTARG};;
		s)
			echo "Option -s triggered, argument = ${OPTARG}"
			sim=${OPTARG};;
		k)
			echo "Option -k triggered, argument = ${OPTARG}"
			clobberness=${OPTARG};;
		d)
			echo "Option -d triggered, argument = ${OPTARG}"
			alt_db=${OPTARG};;
		l)
			echo "Option -l triggered, argument = ${OPTARG}"
			list=${OPTARG};;
		t)
			echo "Option -t triggered, argument = ${OPTARG}"
			analysis_requested=${OPTARG^^};;
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
if [[ -z "${list}" ]] || [[ ! -f "${list}" ]]; then
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

# Number regex to test max concurrent submission parametr
number='^[0-9]+$'

# Checks for proper argumentation
if [[ -z "${clobberness}" ]]; then
	echo "Clobberness is empty...keeping"
	clobberness="keep"
elif [[ "${clobberness}" != "keep" ]] && [[ "${clobberness}" != "clobber" ]]; then
	echo "clobberness not input correctly, must be keep or clobber...keeping"
	clobberness="keep"
fi

# Checks that the gapping is set to one of the csstar presets
if [[ "${gapping}" != "gapped" ]] && [[ "${gapping}" != "ungapped" ]]; then
	echo "gapping does not equal gapped or ungapped...exiting"
	exit 1
fi

# Checks that value given for % Identity is one of the presets for csstar
if [[ "${sim}" != 80 ]] && [[ "${sim}" != 95 ]] && [[ "${sim}" != 98 ]] && [[ "${sim}" != 99 ]] && [[ "${sim}" != 100 ]]; then
	echo "Identity is not one of the presets for csstar and therefore will fail, exiting..."
	exit 1
fi

database_path=${ResGANNCBI_srst2}
database_and_version=${ResGANNCBI_srst2_filename}

if [[ ! -z "${alt_db}" ]]; then
	if [[ ! -f "${alt_db}" ]]; then
		echo " No or empty alternate database location supplied to run_c-sstar_altDB.sh, exiting"
		exit 39
	else
		database_path="${alt_DB}"
		database_basename=$(basename -- "${alt_db}")
		database_basename2=$(echo ${database_basename##*/} | cut -d'.' -f1)
		database_and_version=${database_basename2//_srst2/}
	fi
fi

if [[ -f "${shareScript}/outbreak_analysis.out" ]]; then
	truncate -s 0 "${shareScript}/outbreak_analysis.out"
fi

if [[ -f "${shareScript}/outbreak_analysis.err" ]]; then
	truncate -s 0 "${shareScript}/outbreak_analysis.err"
fi

# Creates the output directory if it does not exist
output_directory=${Phyl_OA}/${analysis_name}
if [[ ! -d ${output_directory} ]]; then
	mkdir -p ${output_directory}
fi

# Clean list of any extra spaces and formatting
"${shareScript}/clean_list.sh" -l "${list}"
cp "${list}" "${output_directory}/${analysis_name}_samples.txt"
rm "${list}.original"
list_file="${output_directory}/${analysis_name}_samples.txt"

if [[ "${analysis_requested}" == "MATRIX" ]] || [[ "${analysis_requested}" == "BOTH" ]]; then
	# Remove any pre-existing files from previous runs
	if [[ -f ${output_directory}/${analysis_name}-mlst_summary.txt ]]; then
		rm ${output_directory}/${analysis_name}-mlst_summary.txt
	fi
	if [[ -f ${output_directory}/${analysis_name}-csstar_summary.txt ]]; then
		rm ${output_directory}/${analysis_name}-csstar_summary.txt
	fi
	if [[ -f ${output_directory}/${analysis_name}-plasmid_summary.txt ]]; then
		rm ${output_directory}/${analysis_name}-plasmid_summary.txt
	fi
	if [[ -f ${output_directory}/${analysis_name}_AR_plasmid_report.tsv ]]; then
		rm ${output_directory}/${analysis_name}_AR_plasmid_report.tsv
	fi
	if [[ -f ${output_directory}/${analysis_name}-sample_summary.txt ]]; then
		rm ${output_directory}/${analysis_name}-sample_summary.txt
	fi
	if [[ -f ${output_directory}/${analysis_name}-GAMA_summary.txt ]]; then
		rm ${output_directory}/${analysis_name}-GAMA_summary.txt
	fi
	if [[ -f ${output_directory}/${analysis_name}-GAMA_rejects.txt ]]; then
		rm ${output_directory}/${analysis_name}-GAMA_rejects.txt
	fi
	if [[ -f ${output_directory}/${analysis_name}-srst2.txt ]]; then
		rm ${output_directory}/${analysis_name}-srst2.txt
	fi
	if [[ -f ${output_directory}/${analysis_name}-srst2_rejects.txt ]]; then
		rm ${output_directory}/${analysis_name}-srst2_rejects.txt
	fi



	# Creates a dictionary to match genes to AR conferred when parsing srst files
	declare -A groups
	echo ""
	echo "Creating AR lookup list from ${local_DBs}/star/group_defs.txt"
	counter=0
	while IFS= read -r line || [ -n "$line" ]; do
		line=${line,,}
		gene=$(echo "${line}" | cut -d ':' -f1)
		first=${gene:0:1}
		if [ "$first" == "#" ]; then
			continue
		fi
		confers=$(echo "${line}" | cut -d ':' -f2)
		groups[${gene}]="${confers}"
		#echo "${counter}:${gene}:${confers}"
		counter=$(( counter + 1))
	done < "${local_DBs}/star/group_defs.txt"

	# Set defaults for checking if all isolates have been compared to the newest ResGANNCBI DB file . If any isolates are not up-to-date, they will be submitted with the appropriate abl_mass_qsub.
	run_csstar="false"
	run_srst2="false"
	run_GAMA="false"
	> "${output_directory}/${analysis_name}-csstar_todo.txt"
	> "${output_directory}/${analysis_name}-srst2_todo.txt"
	> "${output_directory}/${analysis_name}-GAMA_todo.txt"

	# Check that each isolate has been compared to the newest ResGANNCBI DB file
	if [[ "${clobberness}" == "keep" ]]; then
		echo -e "\nMaking sure all isolates use the latest AR Database - ${database_and_version}\n"
		while IFS= read -r line || [ -n "$line" ]; do
			sample_name=$(echo "${line}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
			project=$(echo "${line}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
			OUTDATADIR="${processed}/${project}/${sample_name}"
			echo "checking for ${OUTDATADIR}/c-sstar/${sample_name}.${database_and_version}.${gapping}_${sim}_sstar_summary.txt"
			if [[ -s "${OUTDATADIR}/c-sstar/${sample_name}.${database_and_version}.${gapping}_${sim}_sstar_summary.txt" ]];
			then
				#echo "${project}/${sample_name} has newest ResGANNCBI for normal csstar already"
				:
			else
				echo "${project}/${sample_name} - ccstar needs to be run against ${database_and_version} at ${sim}"
				echo "${project}/${sample_name}" >> "${output_directory}/${analysis_name}-csstar_todo.txt"
				run_csstar="true"
			fi
			echo "checking for ${OUTDATADIR}/srst2/${sample_name}__genes__${database_and_version}_srst2__results.txt"
			if [[ -s ${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq ]] && [[ -s ${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq ]] || [[ -s ${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq.gz ]] && [[ -s ${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fastq.gz ]]; then
				#echo "FASTQs exist"
				if [[ -f "${OUTDATADIR}/srst2/${sample_name}__fullgenes__${database_and_version}_srst2__results.txt" ]] || [[ -f "${OUTDATADIR}/srst2/${sample_name}__genes__${database_and_version}_srst2__results.txt" ]]; then
						#echo "${project}/${sample_name} has newest ResGANNCBI for srst2 already"
						:
					else
						echo "${project}/${sample_name} - SRST2 needs to be run against ${database_and_version}"
						echo "${project}/${sample_name}" >> "${output_directory}/${analysis_name}-srst2_todo.txt"
						run_srst2="true"
				fi
			fi
			echo "checking for ${OUTDATADIR}/GAMA/${sample_name}_${database_and_version}.GAMA"
			if [[ -s "${OUTDATADIR}/GAMA/${sample_name}_${database_and_version}.GAMA" ]];
			then
				mv "${OUTDATADIR}/GAMA/${sample_name}_${database_and_version}.GAMA" "${OUTDATADIR}/GAMA/${sample_name}.${database_and_version}.GAMA"
			fi
			if [[ -s "${OUTDATADIR}/GAMA/${sample_name}.${database_and_version}.GAMA" ]];
			then
				#echo "${project}/${sample_name} has newest ResGANNCBI for normal csstar already"
				:
			else
				echo "${project}/${sample_name} - GAMA needs to be run against ${database_and_version}"
				echo "${project}/${sample_name}" >> "${output_directory}/${analysis_name}-GAMA_todo.txt"
				run_GAMA="true"
			fi
		done < ${list_file}
	else
		run_csstar="true"
		run_srst2="true"
		run_GAMA="true"
		rm "${output_directory}/${analysis_name}-csstar_todo.txt"
		rm "${output_directory}/${analysis_name}-srst2_todo.txt"
		rm "${output_directory}/${analysis_name}-GAMA_todo.txt"
		echo "Copying ${list_file} to ${output_directory}/${analysis_name}_*_todo.txt"
		cp ${list_file} "${output_directory}/${analysis_name}-csstar_todo.txt"
		cp ${list_file} "${output_directory}/${analysis_name}-srst2_todo.txt"
		cp ${list_file} "${output_directory}/${analysis_name}-GAMA_todo.txt"
	fi

	# Creating mashtree of all isolates in list
	echo "Creating mashtree of all samples"
	${shareScript}/mashtree_of_list.sh -i "${list_file}" -d "${output_directory}/mashtree" -o "${analysis_name}" -c "${config}"
	cp "${output_directory}/mashtree/${analysis_name}.dnd" "${output_directory}/${analysis_name}_MASH.newick"
	sed -i "s/_scaffolds_trimmed//g" "${output_directory}/${analysis_name}_MASH.newick"
	rm -r ${output_directory}/mashtree

	# Submits the list of isolates that need the newest ResGANNCBI file for csstar
	if [[ "${run_csstar}" = "true" ]]; then
		echo "Submitting list for csstar qsub analysis"
		qsub -sync y ${shareScript}/abl_mass_qsub_csstar.sh -l "${output_directory}/${analysis_name}-csstar_todo.txt" -m 25 -o "${mass_qsub_folder}" -k "${clobberness}" -s "${sim}" -c "${config}" -d "${database_path}"
	fi
	# Submits the list of isolates that need the newest ResGANNCBI file for srst2
	if [[ "${run_srst2}" = "true" ]]; then
		echo "Submitting list for srst2 qsub analysis"
		qsub -sync y ${shareScript}/abl_mass_qsub_srst2.sh -l "${output_directory}/${analysis_name}-srst2_todo.txt" -m 25 -o "${mass_qsub_folder}" -k "${clobberness}" -c "${config}" -d "${database_path}"
	fi
	if [[ "${run_GAMA}" = "true" ]]; then
		echo "Submitting list for GAMA qsub analysis"
		qsub -sync y ${shareScript}/abl_mass_qsub_GAMA.sh -l "${output_directory}/${analysis_name}-GAMA_todo.txt" -m 25 -o "${mass_qsub_folder}" -k "${clobberness}" -c "${config}" -d "${database_path}"
	fi

	echo $(date)
	sleep 10

	# Loop through and extracts and formats AR genes found in all isolates, as well as the primary MLST type and plasmid replicons. Each are output to separate files. Any AR genes that do not meet the length or % identity are copied to the rejects file.
	while IFS= read -r line; do
		sample_name=$(echo "${line}" | awk -F/ '{ print $2}' | tr -d '[:space:]')
		project=$(echo "${line}" | awk -F/ '{ print $1}' | tr -d '[:space:]')
		OUTDATADIR="${processed}/${project}/${sample_name}"

		if [[ -f "${OUTDATADIR}/c-sstar/${sample_name}.${database_and_version}.${gapping}_${sim}_sstar_summary.txt" ]]; then
			csstar_list=""
			ARDB_full="${OUTDATADIR}/c-sstar/${sample_name}.${database_and_version}.${gapping}_${sim}_sstar_summary.txt"
			#echo "${ARDB_full}"
			# Extracts all AR genes from normal csstar output file and creates a lits of all genes that pass the filtering steps
			while IFS= read -r line; do
				# exit if no genes were found for the sample
				if [[ -z "${line}" ]] || [[ "${line}" == *"No anti-microbial genes were found"* ]]; then
					break
				fi
				IFS='	' read -r -a ar_line <<< "$line"
				percent_ID="${ar_line[6]}"
				percent_length="${ar_line[9]}"
				conferred=$(echo "${ar_line[1]}" | rev | cut -d'_' -f2- | rev)
				contig_number=$(echo "${ar_line[5]}" | rev | cut -d'_' -f3 | rev)
				gene="${ar_line[4]}"
				# Ensure that the gene passes % identity and % length threhsolds for reporting
				if [[ ${percent_length} -ge ${project_parser_Percent_length} ]] && [[ ${percent_ID} -ge ${project_parser_Percent_identity} ]] ; then
					if [[ -z "${csstar_list}" ]]; then
						#	echo "First csstar: ${gene}"
						csstar_list="${gene}(${conferred})[${percent_ID}/${percent_length}:#${contig_number}]"
					else
						if [[ ${csstar_list} == *"${gene}"* ]]; then
							#	echo "${gene} already found in ${csstar_list}"
							:
						else
							#	echo "${gene} not found in ${csstar_list}...adding it"
							csstar_list="${csstar_list},${gene}(${conferred})[${percent_ID}/${percent_length}:#${contig_number}]"
						fi
					fi
					# If length is less than predetermined minimum (90% right now) then the gene is added to a rejects list to show it was outside acceptable limits
				else
					echo -e "${project}\t${sample_name}\tfull_assembly\t${line}" >> ${output_directory}/${analysis_name}-csstar_rejects.txt
				fi
			done < ${ARDB_full}
			if [[ -z "${csstar_list}" ]]; then
				#echo "EMPTY-${project}	${sample_name}	No AR genes discovered"
				echo "${project}	${sample_name}	No AR genes discovered" >> ${output_directory}/${analysis_name}-csstar_summary.txt
				csstar_list="No AR genes discovered"
			else
				#echo "OCCUPADO-${project}	${sample_name}	:${csstar_list}:"
				echo "${project}	${sample_name}	${csstar_list}" >> ${output_directory}/${analysis_name}-csstar_summary.txt
			fi
		else
			echo "${project}	${sample_name}	NO CURRENT FILE" >> ${output_directory}/${analysis_name}-csstar_summary.txt
			csstar_list="NO CURRENT FILE"
		fi


		# Extracts all AR genes from normal csstar output file and creates a lits of all genes that pass the filtering steps
	 	#echo "^^^^^^^^^^^^^^^^^^^ ${OUTDATADIR}/GAMA/${sample_name}.${database_and_version}.GAMA"
		if [[ -f "${OUTDATADIR}/GAMA/${sample_name}.${database_and_version}.GAMA" ]]; then
		    GAMA_list=""
			GARDB_full="${OUTDATADIR}/GAMA/${sample_name}.${database_and_version}.GAMA"
			while IFS= read -r line; do
				# exit if no genes were found for the sample
				if [[ -z "${line}" ]]; then
					break
				elif [[ "${line}" == *"DB	Resistance	Gene_Family	Gene	Contig	Start"* ]]; then
					continue
				fi
				#echo "!!!!!!!!!!!!            ${line}"
				IFS='	' read -r -a ar_line <<< "$line"
				percent_BP_ID=$(echo "${ar_line[11]}" | awk '{ printf "%d", ($1*100) }' )
				percent_codon_ID=$(echo "${ar_line[12]}" | awk '{ printf "%d", ($1*100) }' )
				percent_length=$(echo "${ar_line[13]}" | awk '{ printf "%d", ($1*100) }' )
				conferred=$(echo "${ar_line[1]}" | rev | cut -d'_' -f2- | rev)
				contig_number=$(echo "${ar_line[4]}" | rev | cut -d'_' -f3 | rev)
				gene="${ar_line[3]}"
				# Ensure that the gene passes % identity and % length threhsolds for reporting
				if [[ ${percent_length} -ge ${project_parser_Percent_length} ]] && [[ ${percent_codon_ID} -ge ${project_parser_Percent_identity} ]]; then
					if [[ -z "${GAMA_list}" ]]; then
					#	echo "First GAMA: ${gene}"
						GAMA_list="${gene,,}(${conferred})[${percent_BP_ID}/${percent_codon_ID}/${percent_length}:#${contig_number}]"
					else
						if [[ ${GAMA_list} == *"${gene}"* ]]; then
						#	echo "${gene} already found in ${GAMA_list}"
							:
						else
						#	echo "${gene} not found in ${GAMA_list}...adding it"
							GAMA_list="${GAMA_list},${gene,,}(${conferred})[${percent_BP_ID}/${percent_codon_ID}/${percent_length}:#${contig_number}]"
						fi
					fi
				# If length is less than predetermined minimum (90% right now) then the gene is added to a rejects list to show it was outside acceptable limits
				else
					echo -e "${project}\t${sample_name}\tfull_assembly\t${line}" >> ${output_directory}/${analysis_name}-GAMA_rejects.txt
				fi
			done < ${GARDB_full}
			if [[ -z "${GAMA_list}" ]]; then
				echo "${project}	${sample_name}	No AR genes discovered" >> ${output_directory}/${analysis_name}-GAMA_summary.txt
				GAMA_list="No AR genes discovered"
			else
				echo "${project}	${sample_name}	${GAMA_list}" >> ${output_directory}/${analysis_name}-GAMA_summary.txt
			fi
		else
			echo "${project}	${sample_name}	NO CURRENT FILE" >> ${output_directory}/${analysis_name}-GAMA_summary.txt
			GAMA_list="NO CURRENT FILE"
		fi

		# Adding in srst2 output in a similar fashion as to how the csstar genes are output to the file.
		if [[ -s "${OUTDATADIR}/srst2/${sample_name}__fullgenes__${database_and_version}_srst2__results.txt" ]]; then
			srst2_results=""
			while IFS= read -r line || [ -n "$line" ]; do
			#	echo "Start"
				gene=$(echo "${line}" | cut -d'	' -f3)
				#ODD WAY to do this right now, must look into later, but
				confers=$(echo "${line}" | cut -d'	' -f14 | cut -d';' -f3)
			#	echo "${gene}-${confers}"
				if [[ "${confers}" = "annotation" ]]; then
					continue
				fi
				if [[ -z "${confers}" ]]; then
					if [[ -n ${gene} ]]; then
						if [[ "${gene,,}" == "agly_flqn" ]]; then
							confers="aminoglycoside_and_fluoroquinolone_resistance"
						elif [[ "${gene,,}" == "tetracenomycinc" ]]; then
							confers="tetracenomycinC_resistance"
						else
							confers=${groups[${gene:0:3}]}
						fi
					fi
				fi
				confers=${confers//_resistance/}
				allele=$(echo "${line}" | cut -d'	' -f4 | rev | cut -d'_' -f2- | rev)
				if [[ "${allele}" = "Zn-dependent" ]]; then
					allele="${allele}_hydrolase"
				fi
				coverage=$(echo "${line}" | cut -d'	' -f5)
				depth=$(echo "${line}" | cut -d'	' -f6)
				diffs=$(echo "${line}" | cut -d'	' -f7)
				if [[ ${diffs} == *"trunc"* ]]; then
					allele="TRUNC-${allele}"
				fi
				uncertainty=$(echo "${line}" | cut -d'	' -f8)
				divergence=$(echo "${line}" | cut -d'	' -f9)
				``
				length=$(echo "${line}" | cut -d'	' -f10)
				percent_length=$(echo "$coverage / 1" | bc)
				if [[ "${divergence}" = "0.0" ]]; then
					percent_ID=100
				else
					percent_ID=$(echo "100 - (($divergence + 1) / 1)" | bc)
				fi
			#	echo "${allele}/${coverage}/${depth}/${diffs}/${uncertainty}/${divergence}/${length}/${percent_ID}/${percent_length}"
			# Filter genes based on thresholds for length and percent identity
				if [[ "${percent_ID}" -ge ${project_parser_Percent_identity} ]] && [[ "${percent_length}" -ge ${project_parser_Percent_length} ]]; then
					info_line="${allele}(${confers})[${percent_ID}/${percent_length}]"
					if [[ -z "${srst2_results}" ]]; then
						srst2_results=${info_line,,}
					else
						srst2_results="${srst2_results},${info_line,,}"
					fi
				else
					if [[ ${line} = "Sample	DB	gene"* ]]; then
						:
					else
						echo ${line} >> ${output_directory}/${analysis_name}-srst2_rejects.txt
					fi
				fi
			done < "${OUTDATADIR}/srst2/${sample_name}__fullgenes__${database_and_version}_srst2__results.txt"
			#echo "Test1"
			if [[ -z "${srst2_results}" ]]; then
				#echo "1"
				echo "${project}	${sample_name}	No AR genes discovered" >> ${output_directory}/${analysis_name}-srst2.txt
				srst2_results="No AR genes discovered"
			else
				#echo "2"
				echo "${project}	${sample_name}	${srst2_results}" >> ${output_directory}/${analysis_name}-srst2.txt
			fi
		else
			#echo "3"
			echo "${project}	${sample_name}	NO CURRENT FILE" >> ${output_directory}/${analysis_name}-srst2.txt
			srst2_results="NO CURRENT FILE"
		fi

		# Extracts taxonomic info
		if [[ ! -f "${OUTDATADIR}/${sample_name}.tax" ]]; then
			"${shareScript}/determine_taxID.sh" -n "${sample_name}" -p "${project}" -c "${config}"
		fi
		tax_file="${OUTDATADIR}/${sample_name}.tax"
		sed -i '/^$/d' "${OUTDATADIR}/${sample_name}.tax"
		tax_header=$(head -n1 "${OUTDATADIR}/${sample_name}.tax")
		taxonomy_source_type=$(echo "${tax_header}" | cut -d'(' -f2 | cut -d')' -f1)
		taxonomy_source=$(echo "${tax_header}" | cut -d'-' -f4-)
		#echo "Test-${tax_header};${taxonomy_source_type};${taxonomy_source}"

		#echo "Looking at ${OUTDATADIR}/${sample_name}.tax"
		genus=$(tail -n2 "${OUTDATADIR}/${sample_name}.tax" | head -n1 | cut -d'	' -f2)
		species=$(tail -n1 "${OUTDATADIR}/${sample_name}.tax" | cut -d'	' -f2)
		taxonomy="${genus} ${species}"
		if [[ "${taxonomy_source_type}" = "ANI_REFSEQ_UTD" ]]; then
			confidence_info=$(head -n1 "${taxonomy_source}")
		else
			taxonomy_source_type=$(echo "${taxonomy_source_type}" | cut -d'(' -f2 | cut -d')' -f1)
			confidence_percent=$(echo "${tax_header}" | cut -d'-' -f2)
			confidence_info="NO_ANI...${taxonomy_source_type}=${confidence_percent}"
		fi

		# Quick fix to rename mlst filenames after it was decided that all should be _Pasteur
		#echo -e "\n\n\n\n\n\n\n\n Checking to move ${OUTDATADIR}/MLST/${sample_name}_abaumannii.mlst \n\n\n\n\n\n\n\n"

		if [[ "${taxonomy}" == "Acinetobacter baumannii" ]]; then
			mv "${OUTDATADIR}/MLST/${sample_name}_abaumannii.mlst" "${OUTDATADIR}/MLST/${sample_name}_Oxford.mlst"
			mv "${OUTDATADIR}/MLST/${sample_name}.mlst" "${OUTDATADIR}/MLST/${sample_name}_Pasteur.mlst"
		elif [[ "${taxonomy}" == "Escherichia coli" ]]; then
			mv "${OUTDATADIR}/MLST/${sample_name}.mlst" "${OUTDATADIR}/MLST/${sample_name}_Achtman.mlst"
			mv "${OUTDATADIR}/MLST/${sample_name}_ecoli_2.mlst" "${OUTDATADIR}/MLST/${sample_name}_Pasteur.mlst"
		else
			if [[ -f "${OUTDATADIR}/MLST/${sample_name}.mlst" ]] && [[ ! -f "${OUTDATADIR}/MLST/${sample_name}_Pasteur.mlst" ]]; then
				mv "${OUTDATADIR}/MLST/${sample_name}.mlst" "${OUTDATADIR}/MLST/${sample_name}_Pasteur.mlst"
			fi
		fi

		# Pulls MLST type for sample and adds it to the summary file
		if [[ -f "${OUTDATADIR}/MLST/${sample_name}_Pasteur.mlst" ]]; then
			mlst=$(head -n 1 ${OUTDATADIR}/MLST/${sample_name}_Pasteur.mlst)
			mlst=${mlst//,/\/}
			alleles=$(echo "${mlst}" | cut -d'	' -f4-)
			alleles=${alleles//,/\/}
			echo "${alleles}"
			alleles=${alleles//	/.}
			echo "${alleles}"
			#alleles=${alleles// /.}
			#echo "${alleles}"
			mlst=$(echo "${mlst}" | cut -d'	' -f3)
			if [[ "${mlst}" == "SUB" ]] || [[ "${mlst}" == "AU" ]]; then
				:
			else
				mlst="ST${mlst}"
			fi
		else
			mlst="N/A"
			alleles="N/A"
		fi
		echo -e "${project}\t${sample_name}\t${mlst}\t${alleles}" >> ${output_directory}/${analysis_name}-mlst_summary.txt

		# Pulls Alternate MLST type for sample and adds it to the summary file
		if [[ -f "${OUTDATADIR}/MLST/${sample_name}_Oxford.mlst" ]]; then
			alt_mlst_file="${OUTDATADIR}/MLST/${sample_name}_Oxford.mlst"
		elif [[ -f "${OUTDATADIR}/MLST/${sample_name}_Achtman.mlst" ]]; then
			alt_mlst_file="${OUTDATADIR}/MLST/${sample_name}_Achtman.mlst"
		else
			alt_mlst_file=""
		fi
		if [[ -n "${alt_mlst_file}" ]]; then
			alt_mlst=$(tail -n 1 "${alt_mlst_file}")
			alt_alleles=$(echo "${alt_mlst}" | cut -d'	' -f4-)
			alt_alleles=${alt_alleles//	/.}
			alt_mlst=$(echo "${alt_mlst}" | cut -d'	' -f3)
			alt_mlst=${alt_mlst//,/\/}
			alt_alleles=${alt_alleles//,/\/}
			if [[ "${alt_mlst}" == "SUB" ]] || [[ "${alt_mlst}" == "AU" ]]; then
				:
			else
				alt_mlst="ST${alt_mlst}"
			fi
		else
			alt_mlst="N/A"
			alt_alleles="N/A"
		fi
		echo -e "${project}\t${sample_name}\t${alt_mlst}\t${alt_alleles}" >> ${output_directory}/${analysis_name}-alt_mlst_summary.txt

  	# Print all extracted info to primary file
		echo -e "${project}\t${sample_name}\t${taxonomy}\t${taxonomy_source_type}\t${confidence_info}\t${mlst}\t${alleles}\t${alt_mlst}\t${alt_alleles}\t${csstar_list}\t${srst2_results}\t${GAMA_list}" >> ${output_directory}/${analysis_name}-sample_summary.txt

		# Goes through the plasmid file of the sample and adds all found plasmid replicons to the summary file
		#echo "Starting plasmid extraction"
		if [[ -f ${OUTDATADIR}/plasmidFinder/${sample_name}_results_table_summary.txt ]]; then
			#echo "Found plasmid file"
			:
		fi
		full_contigs=">"
		full_contigs=$(grep -c ${full_contigs} "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta")
		added=0
		while IFS= read -r plasmid || [ -n "$plasmid" ]; do
			line_in=$(echo ${plasmid} | cut -d' ' -f1)
			if [[ "${line_in}" = "No" ]] || [[ "${line_in}" = "Enterococcus,Streptococcus,Staphylococcus" ]] || [[ "${line_in}" = "Enterobacteriaceae" ]] || [[ "${line_in}" = "Plasmid" ]]; then
				# echo "Not using line: $plasmid"
				:
			else
				echo -e "${project}\t${sample_name}\tfull_assembly\t${plasmid}" >> ${output_directory}/${analysis_name}-plasmid_summary.txt
				added=1
			fi
		done < ${OUTDATADIR}/plasmidFinder/${sample_name}_results_table_summary.txt
		if [[ "${added}" -eq 0 ]]; then
			echo -e "${project}\t${sample_name}\tfull_assembly\tNo_Plasmids_Found\t${full_contigs}_contigs-${components}_components" >> ${output_directory}/${analysis_name}-plasmid_summary.txt
		fi
		plas_contigs=">"
		plas_contigs=$(grep -c ${plas_contigs} "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta")
		contig_count=0
		while IFS= read -r contigs || [ -n "$contigs" ]; do
			if [[ "${contigs}" = ">"* ]]; then
				contig_count=$(( contig_count + 1 ))
			fi
		done < ${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta
		while IFS= read -r plasmid || [ -n "$plasmid" ]; do
			line_in=$(echo ${plasmid} | cut -d' ' -f1)
			if [[ "${line_in}" = "No" ]] || [[ "${line_in}" = "Enterococcus,Streptococcus,Staphylococcus" ]] || [[ "${line_in}" = "Enterobacteriaceae" ]] || [[ "${line_in}" = "Plasmid" ]]; then
				:
			else
				echo -e "${project}\t${sample_name}\tplasmid_assembly\t${plasmid}" >> ${output_directory}/${analysis_name}-plasmid_summary.txt
				added=1
			fi
		done < ${OUTDATADIR}/plasmidFinder_on_plasFlow/${sample_name}_results_table_summary.txt

		if [[ "${added}" -eq 0 ]]; then
			echo -e "${project}\t${sample_name}\tplasmid_assembly\tNo_Plasmids_Found\t${plas_contigs}_contigs-${components}_components" >> ${output_directory}/${analysis_name}-plasmid_summary.txt
		fi

	done < ${list_file}

	# Calls script that sorts and formats all isolates info into a matrix for easy viewing
	python3 "${shareScript}/project_parser.py" -s "${output_directory}/${analysis_name}-sample_summary.txt" -p "${output_directory}/${analysis_name}-plasmid_summary.txt" -o "${output_directory}/${analysis_name}_AR_plasmid_report.csv" -d "${database_and_version}"

	declare -a move_list
	move_list=(csstar_todo GAMA_todo srst2_todo alt_mlst_summary csstar_rejects csstar_summary GAMA_rejects GAMA_summary mlst_summary plasmid_summary sample_summary srst2 srst2_rejects)
	if [[ ! -d "${output_directory}/OA_files" ]]; then
		mkdir "${output_directory}/OA_files"
	fi
	for mlist in "${move_list[@]}"; do
	 	mv "${output_directory}/${analysis_name}-${mlist}.txt" "${output_directory}/OA_files/${analysis_name}-${mlist}.txt"
	done
fi

if [[ "${analysis_requested}" == "BOTH" ]] || [[ "${analysis_requested}" == "SNV" ]]; then
	"${shareScript}/SNVPhyl_OA.sh" -l "${list_file}" -a "${analysis_name}" -o "${Phyl_OA}" -c "${config}"
fi



submitter=$(whoami)
global_end_time=$(date "+%m-%d-%Y @ %Hh_%Mm_%Ss")
printf "%s %s" "outbreak_analysis.sh has completed" "${global_end_time}" | mail -s "outbreak analysis complete" "${submitter}@cdc.gov"
