#!/bin/sh -l

#$ -o run_GAMA.out
#$ -e run_GAMA.err
#$ -N run_GAMA
#$ -cwd
#$ -q short.q

#
# Description: Runs the GAMA AR classification tool
#
# Usage: ./run_GAMA.sh -n sample_name -p run_ID [-l] [-d path_to_alt_DB] [-c path_to_config_file]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/GAMA/
#
# Modules required: blat, Python3/3.5.2
#
# v1.0.4 (9/2/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml blat Python3/3.5.4

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_GAMA.sh -n sample_name -p run_ID [-l] [-d path_to_alt_DB] [-c path_to_config_file]"
}

plasmid="false"

# Parse command line options
options_found=0
while getopts ":h?c:p:n:d:l" option; do
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
elif [[ ! -z "${alt_db}" ]]; then
	if [[ ! -f "${alt_db}" ]]; then
		echo " No or empty alternate database location supplied to run_c-sstar_altDB.sh, exiting"
		exit 39
	else
		database_path="${alt_DB}"
		database_basename=$(basename -- "${alt_db}")
		database_basename2=$(echo ${database_basename##*/} | cut -d'.' -f2)
		database_and_version=${database_basename2//_srst2/}
	fi
fi

if [[ "${plasmid}" == "true" ]]; then
	assembly_source="${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta"
	if [ ! -d "$OUTDATADIR/GAMA_plasFlow" ]; then  #create outdir if absent
		echo "Creating $OUTDATADIR/GAMA_plasFlow"
		mkdir -p "$OUTDATADIR/GAMA_plasFlow"
	fi
	OUTDATADIR="${processed}/${project}/${sample_name}/GAMA_plasFlow"
else
	assembly_source="${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta"
	if [ ! -d "$OUTDATADIR/GAMA" ]; then  #create outdir if absent
		echo "Creating $OUTDATADIR/GAMA"
		mkdir -p "$OUTDATADIR/GAMA"
	fi
	OUTDATADIR="${processed}/${project}/${sample_name}/GAMA"
fi


echo "${database_path} - Using DB - ${database_and_version}"


# Create necessary output directories
echo "Running GAMA Antibiotic Resistance Gene Identifier"

### GAMA AR Classifier ### in species mode
python3 GAMA_ResGANNCBI_SciComp_Exe.py -i "${assembly_source}" -d "${database_path}" -o "${OUTDATADIR}/${sample_name}.${database_and_version}.GAMA"

ml -blat -Python3/3.5.4

#Script exited gracefully (unless something else inside failed)
exit 0
