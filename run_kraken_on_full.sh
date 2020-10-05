#!/bin/sh -l

#$ -o run_kraken.out
#$ -e run_kraken.err
#$ -N run_kraken
#$ -cwd
#$ -q short.q

#
# Description: Runs a fasta (assembly) through the kraken classification tool which identifies the most likely
# 	taxonomic classification for the sample on the full kraken database
#
# Usage: ./run_kraken_on_full.sh -n sample_name -p run_ID [-c path_to_config_file]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/kraken_full/
#
# Modules required: kraken/0.10.5, perl/5.12.3, krona/2.7, Python3/3.5.2
#
# v1.0.1 (09/08/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml kraken/0.10.5 perl/5.12.3 krona/2.7 Python3/3.5.2

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_kraken_on_full.sh -n sample_name -p run_ID [-c path_to_config_file]"
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

# Sets the parent output folder as the sample name folder in the processed samples folder in MMB_Data
OUTDATADIR="${processed}/${project}/${sample_name}"

# Creates folder for output from kraken
if [ ! -d "$OUTDATADIR/kraken" ]; then
	echo "Creating $OUTDATADIR/kraken"
	mkdir -p "$OUTDATADIR/kraken/postAssembly_full"
elif [ ! -d "$OUTDATADIR/kraken/postAssembly_full" ]; then
	echo "Creating $OUTDATADIR/kraken/postAssembly_full"
	mkdir -p "$OUTDATADIR/kraken/postAssembly_full"
fi

# Prints out version of kraken
kraken --version
# Status view of run
echo "[:] Running kraken.  Output: ${sample_name}.kraken / ${sample_name}.classified"
# Runs chosen kraken db on the assembly
kraken --db "${kraken_full_db}" --preload --threads "${procs}" --output "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full.kraken" --classified-out "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full.classified" "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta"
# Attempting to weigh contigs and produce standard krona and list output using a modified version of Rich's weighting scripts (will also be done on pure contigs later)
python ${shareScript}/Kraken_Assembly_Converter_2_Exe.py -i "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full.kraken"
kraken-translate --db "${kraken_full_db}" "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full_BP.kraken" > "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full_BP.labels"
kraken-mpa-report --db "${kraken_full_db}" "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full_BP.kraken" > "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full_weighted.mpa"
python3 "${shareScript}/Metaphlan2krona.py" -p "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full_weighted.mpa" -k "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full_weighted.krona"
kraken-report --db "${kraken_full_db}" "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full_BP.kraken" > "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full_BP.list"
python ${shareScript}/Kraken_Assembly_Summary_Exe.py -k "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full_BP.kraken" -l "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full_BP.labels" -t "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full_BP.list" -o "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full_BP_data.list"
#. "${shareScript}/module_changers/perl_5221_to_5123.sh"
ktImportText "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full_weighted.krona" -o "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full_weighted_BP_krona.html"
#. "${shareScript}/module_changers/perl_5123_to_5221.sh"
"${shareScript}/best_hit_from_kraken.sh" -i "${OUTDATADIR}" -r post -t full_BP -k kraken

# Run the metaphlan generator on the kraken output
echo "[:] Generating metaphlan compatible report."
kraken-mpa-report --db "${kraken_full_db}" "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full.kraken" > "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full.mpa"

# Run the krona generator on the metaphlan output
echo "[:] Generating krona output for ${sample_name}."
python3 "${shareScript}/Metaphlan2krona.py" -p "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full.mpa" -k "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full.krona"

# Run the krona graph generator from krona output
#. "${shareScript}/module_changers/perl_5221_to_5123.sh"
ktImportText "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full.krona" -o "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full.html"
#. "${shareScript}/module_changers/perl_5123_to_5221.sh"

# Creates the parsible report from the kraken output
echo "[:] Creating alternate report for taxonomic extraction"
kraken-report --db "${kraken_full_db}" "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full.kraken" > "${OUTDATADIR}/kraken/postAssembly_full/${sample_name}_full.list"

# Parses the output for the best taxonomic hit
echo "[:] Extracting best taxonomic matches"

# Runs the extractor for pulling best hit from a kraken run
"${shareScript}/best_hit_from_kraken.sh" -i "${OUTDATADIR}" -r post -t full -k kraken

ml -kraken/0.10.5 -perl/5.12.3 -krona/2.7 -Python3/3.5.2

#Script exited gracefully (unless something else inside failed)
exit 0
