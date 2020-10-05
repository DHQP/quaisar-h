#!/bin/sh -l

#$ -o run_kraken.out
#$ -e run_kraken.err
#$ -N run_kraken
#$ -cwd
#$ -q short.q

#
# Description: Runs the kraken classification tool which identifies the most likely taxonomic classification for the sample
# 	Can be run using reads or assemblies
#
# Usage: ./run_kraken.sh -n sample_name -p run_ID -r assembly_relation[pre|post] [-c path_to_config_file]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/kraken/
#
# Modules required: kraken/0.10.5, perl/5.12.3, krona/2.7, Python3/3.5.2
#
# v1.0.2 (09/04/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml kraken/0.10.5 perl/5.12.3 krona/2.7 Python3/3.5.2

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_kraken.sh -n sample_name -p run_ID -r assembly_relation[pre|post] [-c path_to_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?c:p:n:r:" option; do
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
		r)
			echo "Option -r triggered, argument = ${OPTARG}"
			assembly_relation=${OPTARG,,};;
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

# Sets the parent output folder as the sample name folder in the processed samples folder in MMB_Data
OUTDATADIR="${processed}/${project}/${sample_name}"


# Checks for proper argumentation
if [[ -z "${sample_name}" ]]; then
	echo "Empty sample name supplied to run_kraken.sh, exiting"
	exit 1
elif [ -z "${assembly_relation}" ]; then
	echo "Empty assembly relativity supplied to run_kraken.sh. Second argument should be 'pre' or 'post' (no quotes). Exiting"
	exit 1
elif [ -z "${project}" ]; then
	echo "Empty project_id supplied to run_kraken.sh. Fourth argument should be the id of the miseq run. Exiting"
	exit 1
fi

if [[ "${assembly_relation}" != "post" ]] && [[ "${assembly_relation}" != "pre" ]]; then
	echo "Improper assembly relation input, must be pre or post"
	exit 487
elif [[ "${assembly_relation}" == "post" ]]; then
	source_type="assembled"
elif [[ "${assembly_relation}" == "pre" ]]; then
	source_type="paired"
fi

# Creates folder for output from kraken
if [ ! -d "$OUTDATADIR/kraken" ]; then
	echo "Creating $OUTDATADIR/kraken"
	mkdir -p "$OUTDATADIR/kraken/${assembly_relation}Assembly"
elif [ ! -d "$OUTDATADIR/kraken/${assembly_relation}Assembly" ]; then
	echo "Creating $OUTDATADIR/kraken/${assembly_relation}Assembly"
	mkdir -p "$OUTDATADIR/kraken/${assembly_relation}Assembly"
fi

# Prints out version of kraken
kraken --version
# Status view of run
echo "[:] Running kraken.  Output: ${sample_name}.kraken / ${sample_name}.classified"
# Runs kraken in paired reads mode
if [ "${source_type}" = "paired" ]; then
	kraken --paired --db "${kraken_mini_db}" --preload --fastq-input --threads "${procs}" --output "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.kraken" --classified-out "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.classified" "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq"
# Runs kraken in single end mode on the concatenated single read file
#elif [ "${source_type}" = "single" ]; then
#	kraken --db "${kraken_mini_db}" --preload --fastq-input --threads "${procs}" --output "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.kraken" --classified-out "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.classified ${OUTDATADIR}/FASTQs/${sample_name}.single.fastq"
# Runs kraken on the assembly
elif [ "${source_type}" = "assembled" ]; then
	kraken --db "${kraken_mini_db}" --preload --threads "${procs}" --output "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.kraken" --classified-out "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.classified" "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta"
	# Attempting to weigh contigs and produce standard krona and list output using a modified version of Rich's weighting scripts (will also be done on pure contigs later)
	#echo "1"
	python3 ${shareScript}/Kraken_Assembly_Converter_2_Exe.py -i "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.kraken"
#	echo "2"
#	kraken-translate --db "${kraken_mini_db}" "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}_BP.kraken" > "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}_BP.labels"
	# Create an mpa report
	#echo "3"
	kraken-mpa-report --db "${kraken_mini_db}" "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}_BP.kraken" > "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}_weighted.mpa"
	# Convert mpa to krona file
	#echo "4"
	python3 "${shareScript}/Metaphlan2krona.py" -p "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}_weighted.mpa" -k "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}_weighted.krona"
	# Create taxonomy list file from kraken file
	#echo "5"
	kraken-report --db "${kraken_mini_db}" "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}_BP.kraken" > "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}_BP.list"
	# Weigh taxonomy list file
#	echo "6"
#	python3 ${shareScript}/Kraken_Assembly_Summary_Exe.py -k "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}_BP.kraken" -l "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}_BP.labels" -t "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}_BP.list" -o "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}_BP_data.list"
	# Run the krona graph generator from krona output
	#echo "7"
	ktImportText "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}_weighted.krona" -o "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}_weighted_BP_krona.html"
	# Runs the extractor for pulling best taxonomic hit from a kraken run
	echo "8"
	"${shareScript}/best_hit_from_kraken.sh" -i ${OUTDATADIR} -r "${assembly_relation}" -t "${source_type}_BP" -k "kraken"
else
	echo "Argument combination is incorrect"
	exit 1
fi

# Run the metaphlan generator on the kraken output
echo "[:] Generating metaphlan compatible report."
kraken-mpa-report --db "${kraken_mini_db}" "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.kraken" > "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.mpa"
# Run the krona generator on the metaphlan output
echo "[:] Generating krona output for ${sample_name}."
# Convert mpa to krona file
python3 "${shareScript}/Metaphlan2krona.py" -p "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.mpa" -k "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.krona"
# Run the krona graph generator from krona output
ktImportText "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.krona" -o "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.html"

# Creates the taxonomy list file from the kraken output
echo "[:] Creating alternate report for taxonomic extraction"
kraken-report --db "${kraken_mini_db}" "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.kraken" > "${OUTDATADIR}/kraken/${assembly_relation}Assembly/${sample_name}_${source_type}.list"
# Parses the output for the best taxonomic hit
echo "[:] Extracting best taxonomic matches"
# Runs the extractor for pulling best taxonomic hit from a kraken run
"${shareScript}/best_hit_from_kraken.sh" -i ${OUTDATADIR} -r "${assembly_relation}" -t "${source_type}" -k "kraken"

ml -kraken/0.10.5 -perl/5.12.3 -krona/2.7 -Python3/3.5.2

#Script exited gracefully (unless something else inside failed)
exit 0
