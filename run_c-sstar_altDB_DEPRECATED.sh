#!/bin/sh -l

#$ -o run_c-sstar.out
#$ -e run_c-sstar.err
#$ -N run_c-sstar
#$ -cwd
#$ -q short.q

#
# Description: Finds anti-microbial resistance genes in the resFinder and ARG-ANNOT databases and exports a file containing list of all genes found using a non-standard srst2 formatted database (see rules for formatting alternate databases in readme ..... eventually)
#
# Usage: ./run_c-sstar_altDB.sh -n sample_name -g run_type(g/u for gapped/ungapped) -s similarity(l/m/h/u/p/o for low(80),medium(95),high(98),ultra-high(99),perfect(100),other(set in config.sh)) -p miseq_run_ID -d path_to_alt_DB [-c path_to_config_file] [-l]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/csstar(_plasFlow)/
#
# Modules required: Python3/3.5.2, ncbi-blast+/LATEST
#
# v1.0.2 (08)/24/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml ncbi-blast+/LATEST Python3/3.5.2

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_c-sstar_altDB.sh -n sample_name -g run_type(g/u for gapped/ungapped) -s similarity(l/m/h/u/p/o for low(80),medium(95),high(98),ultra-high(99),perfect(100),other(set in config.sh)) -p miseq_run_ID -d path_to_alt_DB [-c path_to_config_file]"
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
			echo "Option -c triggered, argument = ${OPTARG}"
			alt_db=${OPTARG};;
		l)
			echo "Option -l triggered, argument = ${OPTARG}"
			plamid="true"
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
	echo "No Project/Run_ID supplied to run_c-sstar_altDB.sh, exiting"
	exit 33
elif [[ -z "${sample_name}" ]]; then
	echo "No sample name supplied to run_c-sstar_altDB.sh, exiting"
	exit 34
if [[ -z "${gapping}" ]] || [[ "${gapping}" != "g"]] && [[ "${gapping}" != "u"]]; then
	echo "No or incorrect gapping provided to run_c-sstar_altDB.sh, exiting"
	exit 37
elif [[ -z "${sim}" ]]; then
	echo "No similarity supplied to run_c-sstar_altDB.sh, exiting"
	exit 38

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
if [[ "${plasmid}" == "true" ]]; then
	if [[ -s "${OUTDATADIR}/plasmidAssembly/${sample_name}_plasmid_scaffolds_trimmed.fasta" ]]; then
		source_assembly="${OUTDATADIR}/plasmidAssembly/${sample_name}_plasmid_scaffolds_trimmed.fasta"
		OUTDATADIR="${OUTDATADIR}/c-sstar_plasmid"
	else
		"No anti-microbial genes were found using c-SSTAR because there were No Plasmids Found" > "${OUTDATADIR}/c-sstar_plasmid/${sample_name}_plasmid_scaffolds_trimmed.fasta"
		exit
	fi
else
	source_assembly="${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta"
	OUTDATADIR="${OUTDATADIR}/c-sstar"
fi

#alt_database="${5##*/}"
alt_database_path=$(basename -- "${alt_db}")
alt_database=$(echo ${alt_database_path##*/} | cut -d'.' -f1)
alt_database=${alt_database//_srst2/}
#echo ${alt_database}

# Creates the output c-sstar folder if it does not exist yet
#echo "${OUTDATADIR}"
if [ ! -d "$OUTDATADIR" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR"
	mkdir -p "$OUTDATADIR"
fi

# Set and call proper version of script based upon if gaps are allowed or not
if [ "${gapping}" == "u" ]; then
	suffix="ungapped"
	if [ ! -d "$OUTDATADIR/${alt_database}" ]; then  #create outdir if absent
		echo "Creating $OUTDATADIR/${alt_database}"
		mkdir -p "$OUTDATADIR/${alt_database}_${suffix}"
	fi
	owd=$(pwd)
	cd "${OUTDATADIR}/${alt_database}_${suffix}"
	echo "Running c-SSTAR on ResGANNCBI DB"
	python "${shareScript}/c-SSTAR_ungapped.py" -g "${source_assembly}" -s "${sim}" -d "${alt_db}" > "${OUTDATADIR}/${alt_database}_${suffix}/${sample_name}.${alt_database}.${suffix}_${sim}.sstar"
elif [ "${gapping}" == "g" ]; then
	suffix="gapped"
	if [ ! -d "$OUTDATADIR/${alt_database}_${suffix}" ]; then  #create outdir if absent
		echo "Creating $OUTDATADIR/${alt_database}_${suffix}"
		mkdir -p "$OUTDATADIR/${alt_database}_${suffix}"
	fi
	owd=$(pwd)
	cd "${OUTDATADIR}/${alt_database}_${suffix}"
	echo "Running c-SSTAR on ResGANNCBI DB"
	python "${shareScript}/c-SSTAR_gapped.py" -g "${source_assembly}" -s "${sim}" -d "${alt_db}" > "${OUTDATADIR}/${alt_database}_${suffix}/${sample_name}.${alt_database}.${suffix}_${sim}.sstar"
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
	label1=$(echo "${line}" | cut -d '	' -f3 | tr '[:upper:]' '[:lower:]')
	label2=$(echo "${line}" | cut -d '	' -f4 | tr '[:upper:]' '[:lower:]')
	info1=""
	#info2=""
	#echo "R1;${label1}-${label2}"
	if [[ "${label1}" = *"TRUNC" ]] && [[ "${label1}" != "str" ]]; then
		#echo "Label 1 was truncated"
		label1="${label1:0:${#label1} - 2}"
		info1="${info1}trunc-"
	fi
	if ( [[ "${label1}" = *"*"* ]] || [[ "${label1}" = *"*" ]] ) && [[ "${label1}" != "str" ]]; then
		#echo "Label 1 is likely novel"
		label1="${label1:0:${#label1} - 1}"
		info1="${info1}novel-"
	fi
	if ( [[ "${label1}" = *"?"* ]] || [[ "${label1}" = *"?" ]] ) && [[ "${label1}" != "str" ]]; then
		#echo "Label 1 is uncertain due to incomplete alignment"
		label1="${label1:0:${#label1} - 1}"
		info1="${info1}alinc-"
	fi
	if ( [[ "${label1}" = *"$"* ]] || [[ "${label1}" = *"$" ]] ) && [[ "${label1}" != "str" ]]; then
		#echo "Label 1 is uncertain due to incomplete alignment"
		label1="${label1:0:${#label1} - 1}"
		info1="${info1}edge-"
	fi
	label1=$(echo "${label1,,}" | tr -d '*?$')
	label2=$(echo "${label2,,}" | tr -d '*?$')
	source=$(echo "${line,,}" | cut -d '	' -f1 | tr -d '[:space:]')
	resistance=$(echo "${line}" | cut -d '	' -f2 | tr -d '[:space:]')
	contig=$(echo "${line}" | cut -d '	' -f5 | tr -d '[:space:]')
	percent=$(echo "${line}" | cut -d '	' -f6 | cut -d'%' -f1 | tr -d '[:space:]')
	len1=$(echo "${line}" | cut -d '	' -f7 | tr -d '[:space:]')
	len2=$(echo "${line}" | cut -d '	' -f8 | tr -d '[:space:]')
	SNPs=$(echo "${line}" | cut -d '	' -f10 | tr -d '[:space:]')
	plen=$(echo "${line}" | cut -d '	' -f9 | tr -d '[:space:]')
	if [[ -z "${info1}" ]]; then
		info1="normal"
	else
		info1=${info1::-1}
	fi
	#printf "%-10s %-50s %-15s %-25s %-25s %-40s %-4s %-5d %-5d %-5d\\n" "${source}" "${resistance}" "${label1}" "${info1}" "${label2}" "${contig}" "${percent}" "${len1}" "${len2}" "${plen}" "${SNPs}"
	echo "${source}	${resistance}	${label1}	${info1}	${label2}	${contig}	${percent}	${len1}	${len2}	${plen}" "${SNPs}"
done < "${OUTDATADIR}/${alt_database}_${suffix}/${sample_name}.${alt_database}.${suffix}_${sim}.sstar" > "${OUTDATADIR}/${alt_database}_${suffix}/${sample_name}.${alt_database}.${suffix}_${sim}.sstar_grouped"
sort -k7,7nr -k10,10nr -k8,8n "${OUTDATADIR}/${alt_database}_${suffix}/${sample_name}.${alt_database}.${suffix}_${sim}.sstar_grouped" > "${OUTDATADIR}/${sample_name}.${alt_database}.${suffix}_${sim}_sstar_summary.txt"

# Catches an empty or missing file
if [ ! -s "${OUTDATADIR}/${sample_name}.${alt_database}.${suffix}_${sim}_sstar_summary.txt" ]; then
	echo "No anti-microbial genes were found using c-SSTAR with both resFinder and ARG-ANNOT DBs" > "${OUTDATADIR}/${sample_name}.${alt_database}.${suffix}_${sim}_sstar_summary.txt"
fi

#Returns to original directory
cd "${owd}"

ml -ncbi-blast+/LATEST -Python3/3.5.2

#Script exited gracefully (unless something else inside failed)
exit 0
