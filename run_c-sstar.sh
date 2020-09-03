#!/bin/sh -l

#$ -o run_c-sstar.out
#$ -e run_c-sstar.err
#$ -N run_c-sstar
#$ -cwd
#$ -q short.q

#
# Description: Finds anti-microbial resistance genes in the resFinder and ARG-ANNOT databases and exports a file containing list of all genes found
#
# Usage: ./run_c-sstar.sh -n sample_name -g run_type(g/u for gapped/ungapped) -s similarity(l/m/h/u/p/o for low(80),medium(95),high(98),ultra-high(99),perfect(100),other(set in config.sh)) -p run_ID [-d path_to_alt_DB] [-c path_to_config_file] [-l]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/csstar_plasFlow/
#
# Modules required: Python3/3.5.2, ncbi-blast+/LATEST
#
# v1.0.3 (09/01/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml ncbi-blast+/LATEST Python3/3.5.2

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_c-sstar.sh -n sample_name -g run_type(g/u for gapped/ungapped) -s similarity(l/m/h/u/p/o for low(80),medium(95),high(98),ultra-high(99),perfect(100),other(set in config.sh)) -p miseq_run_ID [-d path_to_alt_DB] [-c path_to_config_file] [-l]"
}

plasmid="false"

# Parse command line options
options_found=0
while getopts ":h?c:p:n:g:s:d:l" option; do
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
			gapping=${OPTARG};;
		s)
			echo "Option -s triggered, argument = ${OPTARG}"
			sim_letter=${OPTARG};;
		d)
			echo "Option -d triggered, argument = ${OPTARG}"
			alt_db=${OPTARG};;
		l)
			echo "Option -l triggered, argument = ${OPTARG}"
			plasmid="true";;
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

database_path="${ResGANNCBI_srst2}"
database_and_version="${ResGANNCBI_srst2_filename}"

if [[ -z "${project}" ]]; then
	echo "No Project/Run_ID supplied to run_c-sstar_altDB.sh, exiting"
	exit 33
elif [[ -z "${sample_name}" ]]; then
	echo "No sample name supplied to run_c-sstar_altDB.sh, exiting"
	exit 34
elif [[ -z "${gapping}" ]] || [[ "${gapping}" != "g" ]] && [[ "${gapping}" != "u" ]]; then
	echo "No or incorrect gapping provided to run_c-sstar_altDB.sh, exiting"
	exit 37
elif [[ -z "${sim_letter}" ]]; then
	echo "No similarity supplied to run_c-sstar_altDB.sh, exiting"
	exit 38
elif [[ ! -z "${alt_db}" ]]; then
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

# Sets the parent output folder as the sample name folder in the processed samples folder in MMB_Data
OUTDATADIR="${processed}/${project}/${sample_name}"


#Set similarity threshold (set in config.sh) to values given in arguments
if [ "${sim_letter}" == "h" ]; then
	sim=${csstar_high}
elif [ "${sim_letter}" == "l" ]; then
	sim=${csstar_low}
elif [ "${sim_letter}" == "u" ]; then
	sim=${csstar_ultrahigh}
elif [ "${sim_letter}" == "m" ]; then
	sim=${csstar_medium}
elif [ "${sim_letter}" == "p" ]; then
	sim=${csstar_perfect}
elif [ "${sim_letter}" == "o" ]; then
	sim=${csstar_other}
else
	echo "Unknown similarity threshold set (use 'l,m,h,u,or p' for 80,95,98,99,or 100% respectively). Defaulting to 98%"
	sim=${csstar_high}
fi


# Check if there was a request to run it on the plasmid assembly of the sample, change fasta source as necessary
if [[ "${plasmid}" == "true" ]]; then
	if [[ -s "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta" ]]; then
			source_assembly="${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta"
			OUTDATADIR="${OUTDATADIR}/c-sstar_plasFlow"
	else
		echo "Not found: ${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta"
		#echo "No anti-microbial genes were found using c-SSTAR because there were No Plasmids Found" > "${OUTDATADIR}/${database_and_version}_${suffix}/${sample_name}.${database_and_version}.${suffix}_${sim}.sstar"
		exit
	fi
else
	OUTDATADIR="${OUTDATADIR}/c-sstar"
	source_assembly="${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta"
fi


# Creates the output c-sstar folder if it does not exist yet
echo "${OUTDATADIR}"
if [ ! -d "$OUTDATADIR" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR"
	mkdir -p "$OUTDATADIR"
fi

# Set and call proper version of script based upon if gaps are allowed or not
# Calls the ungapped version of csstar
if [ "${gapping}" == "u" ]; then
	suffix="ungapped"
	if [ ! -d "$OUTDATADIR/${database_and_version}" ]; then  #create outdir if absent
		echo "Creating $OUTDATADIR/${database_and_version}"
		mkdir -p "$OUTDATADIR/${database_and_version}_${suffix}"
	fi
	owd=$(pwd)
	cd "${OUTDATADIR}/${database_and_version}_${suffix}"
	echo "Running c-SSTAR on ResGANNCBI DB using"
	echo "python \"${shareScript}/c-SSTAR_ungapped.py\" -g \"${source_assembly}\" -s \"${sim}\" -d \"${database_path}\" > \"${OUTDATADIR}/${database_and_version}_${suffix}/${sample_name}.${database_and_version}.${suffix}_${sim}.sstar\""
	python3 "${shareScript}/c-SSTAR_ungapped.py" -g "${source_assembly}" -s "${sim}" -d "${database_path}" > "${OUTDATADIR}/${database_and_version}_${suffix}/${sample_name}.${database_and_version}.${suffix}_${sim}.sstar"
# Calls the gapped version of csstar
elif [ "${gapping}" == "g" ]; then
	suffix="gapped"
	if [ ! -d "$OUTDATADIR/${database_and_version}_${suffix}" ]; then  #create outdir if absent
		echo "Creating $OUTDATADIR/${database_and_version}_${suffix}"
		mkdir -p "$OUTDATADIR/${database_and_version}_${suffix}"
	fi
	owd=$(pwd)
	cd "${OUTDATADIR}/${database_and_version}_${suffix}"
	echo "Running c-SSTAR on ResGANNCBI DB"
	echo "python \"${shareScript}/c-SSTAR_gapped.py\" -g \"${source_assembly}\" -s \"${sim}\" -d \"${database_path}\" > \"${OUTDATADIR}/${database_and_version}_${suffix}/${sample_name}.${database_and_version}.${suffix}_${sim}.sstar\""
	python3 "${shareScript}/c-SSTAR_gapped.py" -g "${source_assembly}" -s "${sim}" -d "${database_path}" > "${OUTDATADIR}/${database_and_version}_${suffix}/${sample_name}.${database_and_version}.${suffix}_${sim}.sstar"
# Unknown gapping parameter when called (not 'g' or 'u')
else
	echo "Unknown run type set (only use 'g' or 'u' for gapped/ungapped analysis"
	exit 1
fi



###################################### FIND WAY TO CATCH FAILURE !!!!!!!!!! ###############################

# Goes through ResGANNCBI outfile and adds labels as well as resistance conferred to the beginning of the line
# Takes .sstar file in and outputs as .sstar_grouped
while IFS= read -r line || [ -n "$line" ]; do
	line=${line}
	#echo ${line}
	# Extract gene (label1) and allele (label2) from line, also force all characters to be lowercase
	label1=$(echo "${line}" | cut -d '	' -f3 | tr '[:upper:]' '[:lower:]')
	label2=$(echo "${line}" | cut -d '	' -f4 | tr '[:upper:]' '[:lower:]')
	# Determine what flags were thrown for this gene by csstar
	info1=""
	# Truncated allele
	if [[ "${label1}" = *"TRUNC" ]] && [[ "${label1}" != "str" ]]; then
		#echo "Label 1 was truncated"
		label1="${label1:0:${#label1} - 2}"
		info1="${info1}trunc-"
	fi
	# Likely novel allele
	if ( [[ "${label1}" = *"*"* ]] || [[ "${label1}" = *"*" ]] ) && [[ "${label1}" != "str" ]]; then
		#echo "Label 1 is likely novel"
		label1="${label1:0:${#label1} - 1}"
		info1="${info1}novel-"
	fi
	# Incomplete alignment length, Uncertainy exists in one allele
	if ( [[ "${label1}" = *"?"* ]] || [[ "${label1}" = *"?" ]] ) && [[ "${label1}" != "str" ]]; then
		#echo "Label 1 is uncertain due to incomplete alignment"
		label1="${label1:0:${#label1} - 1}"
		info1="${info1}alinc-"
	fi
	# Incomplete alignment length at edge
	if ( [[ "${label1}" = *"$"* ]] || [[ "${label1}" = *"$" ]] ) && [[ "${label1}" != "str" ]]; then
		#echo "Label 1 is uncertain due to incomplete alignment"
		label1="${label1:0:${#label1} - 1}"
		info1="${info1}edge-"
	fi
	# Removes character add-ons of genes and alleles, also lower cases all characters for searching later
	label1=$(echo "${label1,,}" | tr -d '*?$')
	label2=$(echo "${label2,,}" | tr -d '*?$')
	# Extract source database that AR gene match came from
	source=$(echo "${line,,}" | cut -d '	' -f1 | tr -d '[:space:]')
	# Extract the type of resistance that is conferred by the gene
	resistance=$(echo "${line}" | cut -d '	' -f2 | tr -d '[:space:]')
	# Trim contig identifier of spaces
	contig=$(echo "${line}" | cut -d '	' -f5 | tr -d '[:space:]')
	# Extract % from line
	percent=$(echo "${line}" | cut -d '	' -f6 | cut -d'%' -f1 | tr -d '[:space:]')
	# Determine length of query and subject sequences
	len1=$(echo "${line}" | cut -d '	' -f7 | tr -d '[:space:]')
	len2=$(echo "${line}" | cut -d '	' -f8 | tr -d '[:space:]')
	plen=$(echo "${line}" | cut -d '	' -f9 | tr -d '[:space:]')
	if [[ -z "${info1}" ]]; then
		info1="normal"
	else
		info1=${info1::-1}
	fi
	#printf "%-10s %-50s %-15s %-25s %-25s %-40s %-4s %-5d %-5d %-5d\\n" "${source}1" "${resistance}2" "${label1}3" "${info1}4" "${label2}5" "${contig}A" "${percent}B" "${len1}C" "${len2}D" "${plen}E"
	echo "${source}	${resistance}	${label1}	${info1}	${label2}	${contig}	${percent}	${len1}	${len2}	${plen}"
done < "${OUTDATADIR}/${database_and_version}_${suffix}/${sample_name}.${database_and_version}.${suffix}_${sim}.sstar" > "${OUTDATADIR}/${database_and_version}_${suffix}/${sample_name}.${database_and_version}.${suffix}_${sim}.sstar_grouped"
# Writes all AR genes to file based on %ID, %length, and finally length of gene
sort -k7,7nr -k10,10nr -k8,8n "${OUTDATADIR}/${database_and_version}_${suffix}/${sample_name}.${database_and_version}.${suffix}_${sim}.sstar_grouped" > "${OUTDATADIR}/${sample_name}.${database_and_version}.${suffix}_${sim}_sstar_summary.txt"

# Catches an empty or missing file, adding that no AMR genes were found if no file was created
if [ ! -s "${OUTDATADIR}/${sample_name}.${database_and_version}.${suffix}_${sim}_sstar_summary.txt" ]; then
	echo "No anti-microbial genes were found using c-SSTAR with both resFinder and ARG-ANNOT DBs" > "${OUTDATADIR}/${sample_name}.${database_and_version}.${suffix}_${sim}_sstar_summary.txt"
fi

#Returns to original directory
cd "${owd}"

ml -ncbi-blast+/LATEST -Python3/3.5.2

#Script exited gracefully (unless something else inside failed)
exit 0
