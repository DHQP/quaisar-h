#!/bin/sh -l

#$ -o run_GAMA.out
#$ -e run_GAMA.err
#$ -N run_GAMA
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Runs the GAMA AR classification tool
#
# Usage: ./run_GAMA.sh sample_name run_ID -c|p [path_to_alt_DB]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/GAMA/
#
# Modules required: blat, Python3/3.5.2
#
# v1.0.3 (1/17/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml blat Python3/3.5.4

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_GAMA.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_GAMA.sh   sample_name    run_ID	c|p	[path_to_alt_DB]"
	echo "Output is saved to ${processed}/miseq_run_ID/sample_name/GAMA/"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty project id name supplied to run_GAMA.sh, exiting"
	exit 1
elif [ -z "$3" ]; then
	echo "Empty assembly source supplied to run_GAMA.sh (must be -c or -p, chromosome or plasmid respectively), exiting"
	exit 1
elif [[ "${3}" != "-c" &&  "${3}" != "-p" ]]; then
	echo "Incorrect assembly source supplied to run_GAMA.sh (must be -c or -p, chromosome or plasmid respectively), exiting"
	exit 1
elif [ ! -z "${4}" ]; then
	ARDB="${4}"
	ARDB_short=$(echo "${ARDB}" | rev | cut -d'/' -f1 | rev | cut -d'.' -f1)
else
	ARDB="${ResGANNCBI_srst2}"
	ARDB_short="${ResGANNCBI_srst2_filename}"
fi

echo "${4} - Using DB - ${ARDB}"

# Sets the output folder of GAMA classifier to the GAMA folder under the sample_name folder in processed samples
OUTDATADIR="${processed}/${2}/${1}"

# Create necessary output directories
echo "Running GAMA Antibiotic Resistance Gene Identifier"

OUTDATA="${OUTDATADIR}"

if [[ "${3}" == "-c" ]]; then
	assembly_source="${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta"
	if [ ! -d "$OUTDATADIR/GAMA" ]; then  #create outdir if absent
		echo "Creating $OUTDATADIR/GAMA"
		mkdir -p "$OUTDATADIR/GAMA"
	fi
	OUTDATADIR="${OUTDATADIR}/GAMA"
elif [[ "${3}" == "-p" ]]; then
	assembly_source="${OUTDATADIR}/plasFlow/Unicycler_assemblies/${1}_uni_assembly/${1}_plasmid_assembly_trimmed.fasta"
	if [ ! -d "$OUTDATADIR/GAMA_plasFlow" ]; then  #create outdir if absent
		echo "Creating $OUTDATADIR/GAMA_plasFlow"
		mkdir -p "$OUTDATADIR/GAMA_plasFlow"
	fi
	OUTDATADIR="${OUTDATADIR}/GAMA_plasFlow"
else
	echo "Unknown Assembly source identifier, exiting"
	exit 5564
fi
### GAMA AR Classifier ### in species mode
python3 GAMA_ResGANNCBI_SciComp_Exe.py -i "${assembly_source}" -d "${ARDB}" -o "${OUTDATADIR}/${1}.${ARDB_short}.GAMA"

ml -blat -Python3/3.5.4

#Script exited gracefully (unless something else inside failed)
exit 0
