#!/bin/sh -l

#$ -o run_gottcha.out
#$ -e run_gottcha.err
#$ -N run_gottcha
#$ -cwd
#$ -q short.q

#
# Description: Runs the gottcha classification tool (now only at species level) which identifies the most likely taxonomic classification for the sample
#
# Usage: ./run_gottcha.sh -n sample_name -p run_ID [-c path_to_config_file]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/gottcha/
#
# Modules required: gottcha, perl/5.12.3
#
# v1.0.2 (9/04/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml perl/5.12.3  gottcha

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_gottcha.sh -n sample_name -p run_ID [-c path_to_config_file]"
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

if [[ -z "${sample_name}" ]]; then
	echo "Empty sample name supplied to run_gottcha.sh, exiting"
	exit 1
elif [ -z "${project}" ]; then
	echo "Empty project id name supplied to run_gottcha.sh, exiting"
	exit 1
fi

# Sets the output folder of gottcha classifier to the gottcha folder under the sample_name folder in processed samples
OUTDATADIR="${processed}/${project}/${sample_name}"

# Create necessary output directories
echo "Running gottcha Taxonomic Classifier"
if [ ! -d "$OUTDATADIR/gottcha" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR/gottcha"
	mkdir -p "$OUTDATADIR/gottcha"
fi
if [ ! -d "$OUTDATADIR/gottcha/gottcha_S" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR/gottcha/gottcha_S"
	mkdir -p "$OUTDATADIR/gottcha/gottcha_S"
fi

if [[ ! -f "${OUTDATADIR}/trimmed/${sample_name}.paired.fq" ]]; then
	if [[ -f "${OUTDATADIR}/trimmed/${sample_name}.paired.fq.gz" ]]; then
		gunzip -c "${OUTDATADIR}/trimmed/${sample_name}.paired.fq.gz" > "${OUTDATADIR}/trimmed/${sample_name}.paired.fq"
	else
		if [[ -f "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" ]]; then
			:
		elif [[ -f "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]]; then
			gunzip -c "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq.gz" > "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq"
		else
			echo "R1_001.paired.fq does not exist (either zipped or not)"
			exit
		fi
		if [[ -f "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" ]]; then
			:
		elif [[ -f "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]]; then
			gunzip -c "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq.gz" > "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq"
		else
			echo "R2_001.paired.fq does not exist (either zipped or not)"
			exit
		fi
		cat "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq" "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq" > "${OUTDATADIR}/trimmed/${sample_name}.paired.fq"
		if [[ -s "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]]; then
			rm "${OUTDATADIR}/trimmed/${sample_name}_R2_001.paired.fq"
		fi
		if [[ -s "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]]; then
			rm "${OUTDATADIR}/trimmed/${sample_name}_R1_001.paired.fq.gz"
		fi
	fi
fi


### Gottcha Taxonomy Classifier ### in species mode
gottcha.pl --mode all --outdir "${OUTDATADIR}/gottcha/gottcha_S" --input "${OUTDATADIR}/trimmed/${sample_name}.paired.fq" --database "${gottcha_db}"

gzip "${OUTDATADIR}/trimmed/${sample_name}.paired.fq"


# Create the krona graphs from each of the analyses
ml krona

ktImportText "${OUTDATADIR}/gottcha/gottcha_S/${sample_name}_temp/${sample_name}.lineage.tsv" -o "${OUTDATADIR}/gottcha/${sample_name}_species.krona.html"

ml -krona
#Create a best hit from gottcha1 file
"${shareScript}/best_hit_from_gottcha1.sh" -i "${OUTDATADIR}"

ml -perl/5.12.3 -gottcha

#Script exited gracefully (unless something else inside failed)
exit 0
