#!/bin/sh -l

#$ -o kraken_weigh_contigs.out
#$ -e kraken_weigh_contigs.err
#$ -N kraken_weigh_contigs
#$ -cwd
#$ -q short.q

#
# Description: Grabs the best species match based on %/read hits from the kraken tool run
#
# Usage: ./kraken_weigh_contigs.sh -s sample_name -p run_ID [-v kraken_version ([kraken|kraken2)] [-c path_to_config_file]
#
# Output location: config.sh_output_location/run_ID/sample_name/kraken(2)/
#
# No modules required
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./kraken_weigh_contigs.sh -s sample_name -p run_ID [-v kraken_version ([kraken|kraken2)] [-c path_to_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?v:s:p:c:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		s)
			echo "Option -s triggered, argument = ${OPTARG}"
			sample_name=${OPTARG};;
		v)
			echo "Option -v triggered, argument = ${OPTARG}"
			version=${OPTARG};;
		p)
			echo "Option -p triggered, argument = ${OPTARG}"
			project=${OPTARG};;
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
if [[ -z "${sample_name}" ]]; then
	echo "No sample name supplied, exiting"
	exit 1
elif [[ -z "${project}" ]]; then
	echo "No project ID supplied, exiting"
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

if [[ -z "${version}" ]]; then
	version="kraken"
elif [[ "${version}" != "kraken" ]] && [[ "${version}" != "kraken2" ]]; then
	echo "Invalid kraken version supplied, exiting"
	exit 4
fi

#Sets output folder to the correct path relative to assembly completion
OUTDATADIR="${processed}/${project}/${sample_name}/${version}/postAssembly"
echo "-${OUTDATADIR}-"

#Checks to see if the list file used for calculations exists and exits if it does not
if [[ ! -s "${OUTDATADIR}/${sample_name}_assembled.${version}" ]]; then
	echo "${OUTDATADIR}/${sample_name}_assembled.${version} does not exist"
	exit 1
fi

contig_sizes=()
total_size=0
unclassified=0

sort -t$'\t' -k4,4 -n "${OUTDATADIR}/${sample_name}_assembled.${version}" > "${OUTDATADIR}/${sample_name}_assembled_sorted.${version}"

#Parses the kraken output list line by line
while IFS= read -r line  || [ -n "$line" ]; do
		prefix=$(echo "${line}" | cut -d'	' -f1,2,3,4)
		#echo "${prefix}"
		classified=$(echo "${line}" | cut -d'	' -f1)
		#echo "classified as:${classified}"
		if [[ "${classified}" == "C" ]]; then
			contig_size=$(echo "${line}" | cut -d'	' -f4)
			contig_sizes+=(${contig_size})
			total_size=$(( total_size + contig_size ))
		else
			echo "Contig not classified"
			unclassified=$(( unclassified + 1 ))
		fi
done < "${OUTDATADIR}/${sample_name}_assembled_sorted.${version}"

contig_count=${#contig_sizes[@]}
counter=0
#for contiggy in ${contig_sizes[@]}; do
#	echo "${counter}-${contiggy}"
#	counter=$((counter + 1 ))
#done

smallest=${contig_sizes[0]}
adjusted_contig_count=$(( total_size / smallest ))
echo "Contig count = ${contig_count}"
echo "Total Size = ${total_size}"
echo "unclassified = ${unclassified}"
echo "Smallest contig = ${smallest}"
echo "Adjusted contig count = ${adjusted_contig_count}"

while IFS= read -r line  || [ -n "$line" ]; do
		IFS='	' read -a arr_line <<< "$line"
		#echo "Size:${#arr_line}"
		original_size=${arr_line[3]}
		#echo "3-${arr_line[3]}"
		#echo "OS = ${original_size}, smallest = ${smallest}"
		adjusted_size=$(( original_size / smallest ))
		#echo "Adj_size = ${adjusted_size}"
		arr_line[3]=${adjusted_size}
		echo -e "${arr_line[0]}\t${arr_line[1]}\t${arr_line[2]}\t${arr_line[3]}"

done < "${OUTDATADIR}/${sample_name}_assembled_sorted.${version}"

if [[ ! -s "${OUTDATADIR}/${sample_name}_assembled_weighted.mpa" ]]; then
	echo "${OUTDATADIR}/${sample_name}_assembled_weighted.mpa does not exist, cant do mpa adjustment"
	exit 1
fi

#Script exited gracefully (unless something else inside failed)
exit 0
