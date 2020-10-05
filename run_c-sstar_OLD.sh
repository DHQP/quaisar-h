#!/bin/sh -l

#$ -o run_c-sstar.out
#$ -e run_c-sstar.err
#$ -N run_c-sstar
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Finds anti-microbial resistance genes in the resFinder and ARG-ANNOT databases and exports a file containing list of all genes found
#
# Usage: ./run_c-sstar.sh sample_name run_type(g/u for gapped/ungapped) similarity(l/m/h/u/p/o for low(80),medium(95),high(98),ultra-high(99),perfect(100),other(set in config.sh)) run_ID [-p])
# 	-p flag is to run it on the plasFlow assembly, assuming it is in the default config location
#
# Output location: default_config.sh_output_location/run_ID/sample_name/csstar(_plasFlow)/
#
# Modules required: Python/3.5.2, ncbi-blast+/LATEST
#
# v1.0.1 (10/29/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml Python3/3.5.2 ncbi-blast+/LATEST

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_c-sstar.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_c-sstar.sh   sample_name   run-type[g/u](for gapped/ungapped)   similarity[l/m/h/u/p/o](for low/medium/high/ultra-high/perfect as 80/95/98/99/100, other(st in config.sh) run_ID [-p]"
	echo "Output is saved to ${processed}/sample_name/c-sstar"
	exit
elif [ -z "$2" ]; then
	echo "Empty run type supplied to run_c-sstar.sh, exiting"
	exit 1
elif [ -z "$3" ]; then
	echo "Empty similarity supplied to run_c-sstar.sh, exiting"
	exit 1
elif [ -z "$4" ]; then
	echo "Empty project id supplied to run_c-sstar.sh, exiting"
	exit 1
fi

# Sets the parent output folder as the sample name folder in the processed samples folder in MMB_Data
OUTDATADIR="${processed}/${4}/${1}"


#Set similarity threshold (set in config.sh) to values given in arguments
if [ "${3}" == "h" ]; then
	sim=${csstar_high}
elif [ "${3}" == "l" ]; then
	sim=${csstar_low}
elif [ "${3}" == "u" ]; then
	sim=${csstar_ultrahigh}
elif [ "${3}" == "m" ]; then
	sim=${csstar_medium}
elif [ "${3}" == "p" ]; then
	sim=${csstar_perfect}
elif [ "${3}" == "o" ]; then
	sim=${csstar_other}
else
	echo "Unknown similarity threshold set (use 'l,m,h,u,or p' for 80,95,98,99,or 100% respectively). Defaulting to 98%"
	sim=${csstar_high}
fi
# Check if there was a request to run it on the plasmid assembly of the sample, change fasta source as necessary
if [[ "${5}" == "--plasmid" ]] || [[ "${5}" == "-p" ]]; then
	if [[ -s "${OUTDATADIR}/plasmidAssembly/${1}_plasmid_scaffolds_trimmed.fasta" ]]; then
	#if [[ -s "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/assembly.fasta" ]]; then
		source_assembly="${OUTDATADIR}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_trimmed.fasta"
		OUTDATADIR="${OUTDATADIR}/c-sstar_plasFlow"
		#source_assembly="${OUTDATADIR}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/assembly.fasta"
		#OUTDATADIR="${OUTDATADIR}/plasFlow_plasmid"
	else
		if [[ "${2}" = "g" ]]; then
			suffix="gapped"
		elif [[ "${2}" = "u" ]]; then
			suffix="ungapped"
		fi
		"No anti-microbial genes were found using c-SSTAR because there were No Plasmids Found" > "${OUTDATADIR}/${ResGANNCBI_srst2_filename}_${suffix}/${1}.${ResGANNCBI_srst2_filename}.${suffix}_${sim}.sstar"
		exit
	fi
else
	if [[ ! -s "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta" ]]; then
		echo "No Assembly found to run c-sstar with (${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta does not exist)"
		exit
	fi
	source_assembly="${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta"
	OUTDATADIR="${OUTDATADIR}/c-sstar"
fi



# Creates the output c-sstar folder if it does not exist yet
echo "${OUTDATADIR}"
if [ ! -d "$OUTDATADIR" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR"
	mkdir -p "$OUTDATADIR"
fi

# Set and call proper version of script based upon if gaps are allowed or not
# Calls the ungapped version of csstar
if [ "${2}" == "u" ]; then
	suffix="ungapped"
	if [ ! -d "$OUTDATADIR/${ResGANNCBI_srst2_filename}" ]; then  #create outdir if absent
		echo "Creating $OUTDATADIR/${ResGANNCBI_srst2_filename}"
		mkdir -p "$OUTDATADIR/${ResGANNCBI_srst2_filename}_${suffix}"
	fi
	owd=$(pwd)
	cd "${OUTDATADIR}/${ResGANNCBI_srst2_filename}_${suffix}"
	echo "Running c-SSTAR on ResGANNCBI DB using"
	echo "python \"${shareScript}/c-SSTAR_ungapped.py\" -g \"${source_assembly}\" -s \"${sim}\" -d \"${ResGANNCBI_srst2}\" > \"${OUTDATADIR}/${ResGANNCBI_srst2_filename}_${suffix}/${1}.${ResGANNCBI_srst2_filename}.${suffix}_${sim}.sstar\""
	python3 "${shareScript}/c-SSTAR_ungapped.py" -g "${source_assembly}" -s "${sim}" -d "${ResGANNCBI_srst2}" > "${OUTDATADIR}/${ResGANNCBI_srst2_filename}_${suffix}/${1}.${ResGANNCBI_srst2_filename}.${suffix}_${sim}.sstar"
# Calls the gapped version of csstar
elif [ "${2}" == "g" ]; then
	suffix="gapped"
	if [ ! -d "$OUTDATADIR/${ResGANNCBI_srst2_filename}_${suffix}" ]; then  #create outdir if absent
		echo "Creating $OUTDATADIR/${ResGANNCBI_srst2_filename}_${suffix}"
		mkdir -p "$OUTDATADIR/${ResGANNCBI_srst2_filename}_${suffix}"
	fi
	owd=$(pwd)
	cd "${OUTDATADIR}/${ResGANNCBI_srst2_filename}_${suffix}"
	echo "Running c-SSTAR on ResGANNCBI DB"
	echo "python \"${shareScript}/c-SSTAR_gapped.py\" -g \"${source_assembly}\" -s \"${sim}\" -d \"${ResGANNCBI_srst2}\" > \"${OUTDATADIR}/${ResGANNCBI_srst2_filename}_${suffix}/${1}.${ResGANNCBI_srst2_filename}.${suffix}_${sim}.sstar\""
	python3 "${shareScript}/c-SSTAR_gapped.py" -g "${source_assembly}" -s "${sim}" -d "${ResGANNCBI_srst2}" > "${OUTDATADIR}/${ResGANNCBI_srst2_filename}_${suffix}/${1}.${ResGANNCBI_srst2_filename}.${suffix}_${sim}.sstar"
# Unknown gapping parameter when called (not 'g' or 'u')
else
	echo "Unknown run type set (only use 'g' or 'u' for gapped/ungapped analysis"
	exit 1
fi



###################################### FIND WAY TO CATCH FAILURE !!!!!!!!!! ###############################

# Goes through ResGANNCBI outfile and adds labels as well as resistance conferred to the beginning of the line
# Takes .sstar file in and outputs as .sstar_grouped
while IFS= read -r line || [ -n "$line" ]; do

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
	# Catch instances where match length is longer than gene (gaps cause extension)
	#if [[ ${len1} -ge ${len2} ]]; then
	#	plen=100
	# Determine % length match
	#else
	#	plen=$( echo "${len1} ${len2}" | awk '{ printf "%d", ($1*100)/$2 }' )
	#fi
	# Check and display any flags found, otherwise mark it as normal
	if [[ -z "${info1}" ]]; then
		info1="normal"
	else
		info1=${info1::-1}
	fi
	#printf "%-10s %-50s %-15s %-25s %-25s %-40s %-4s %-5d %-5d %-5d\\n" "${source}1" "${resistance}2" "${label1}3" "${info1}4" "${label2}5" "${contig}A" "${percent}B" "${len1}C" "${len2}D" "${plen}E"
	echo "${source}	${resistance}	${label1}	${info1}	${label2}	${contig}	${percent}	${len1}	${len2}	${plen}"
done < "${OUTDATADIR}/${ResGANNCBI_srst2_filename}_${suffix}/${1}.${ResGANNCBI_srst2_filename}.${suffix}_${sim}.sstar" > "${OUTDATADIR}/${ResGANNCBI_srst2_filename}_${suffix}/${1}.${ResGANNCBI_srst2_filename}.${suffix}_${sim}.sstar_grouped"
# Writes all AR genes to file based on %ID, %length, and finally length of gene
sort -k7,7nr -k10,10nr -k8,8n "${OUTDATADIR}/${ResGANNCBI_srst2_filename}_${suffix}/${1}.${ResGANNCBI_srst2_filename}.${suffix}_${sim}.sstar_grouped" > "${OUTDATADIR}/${1}.${ResGANNCBI_srst2_filename}.${suffix}_${sim}_sstar_summary.txt"

# Catches an empty or missing file, adding that no AMR genes were found if no file was created
if [ ! -s "${OUTDATADIR}/${1}.${ResGANNCBI_srst2_filename}.${suffix}_${sim}_sstar_summary.txt" ]; then
	echo "No anti-microbial genes were found using c-SSTAR with both resFinder and ARG-ANNOT DBs" > "${OUTDATADIR}/${1}.${ResGANNCBI_srst2_filename}.${suffix}_${sim}_sstar_summary.txt"
fi

#Returns to original directory
cd "${owd}"

ml -ncbi-blast+/LATEST -Python3/3.5.2

#Script exited gracefully (unless something else inside failed)
exit 0
