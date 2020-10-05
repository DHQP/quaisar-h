#!/bin/sh -l

#$ -o do_busco.out
#$ -e do_busco.err
#$ -N do_busco
#$ -cwd
#$ -q short.q

#
# Description: A script that takes a sample and compares it to a busco database to discover number of similar genes (% conserved proteins) from prokka output
#
# Usage ./do_busco.sh -n sample_name -d DB(to search against) -p run_ID [-c path_config_file]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/BUSCO/
#
# Modules required: busco/3.0.1, Python3/3.5.4 (whatever version used to install it, must have pipebricks)
#
# v1.0.3 (08/18/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml busco/3.0.2 Python3/3.6.1
python3 -V

#  Function to print out help blurb
show_help () {
	echo "Usage is ./do_busco.sh -n sample_name -p project_ID -d database_name [-c path_to_config_file]"
	echo "Output is saved to ${processed}/run_ID/sample_name/ where processed is retrieved from config file, either default or imported"
}

# Parse command line options
options_found=0
while getopts ":h?n:p:d:c:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		n)
			echo "Option -n triggered, argument = ${OPTARG}"
			sample_name=${OPTARG};;
		p)
			echo "Option -p triggered, argument = ${OPTARG}"
			project=${OPTARG};;
		d)
			echo "Option -d triggered, argument = ${OPTARG}"
			DB=${OPTARG};;
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
if [ -z "${DB}" ]; then
	echo "Empty database name supplied, exiting"
	exit 1
elif [ ! -s "${local_DBs}/BUSCO/${DB,}" ]; then
	echo "Exists? - ${local_DBs}/BUSCO/${DB,}"
	echo "The taxon does not exist in the BUSCO database. This will be noted and the curator of the database will be notified. However, since nothing can be done at the moment....exiting"
	# Create a dummy folder to put non-results into (if it doesnt exist)
	if [ ! -d "${processed}/${project}/${sample_name}/BUSCO" ]; then  #create outdir if absent
		echo "${processed}/${project}/${sample_name}/BUSCO"
		mkdir -p "${processed}/${project}/${sample_name}/BUSCO"
	fi
	# Write non-results to a file in busco folder
	echo "No matching DB database found for ${sample_name}" >> "${processed}/${project}/${sample_name}/BUSCO/summary_${1}.txt"
	# Add genus to list to download and to database
	global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
	exit 1
elif [ -z "${project}" ]; then
	echo "Empty project id supplied to do_busco.sh, exiting"
	exit 1
fi

# Sets output folder to parent isolate folder
OUTDATADIR="${processed}/${project}/${sample_name}"
# Sets buscoDB to the folder matching the 2nd argument (dbs vary on taxonomic level)
buscoDB="${local_DBs}/BUSCO/${DB,}"

### BUSCO % Conserved Protein Identity ###
echo "Running BUSCO for Conserved Protein Identity"
#Deletes old busco folder if existent
if [ -d "$OUTDATADIR/BUSCO" ]; then
	echo "Deleting old BUSCO folder at $OUTDATADIR/BUSCO"
	rm -rf "$OUTDATADIR/BUSCO"
fi
# Creates busco output folder if not already existent
echo "Creating $OUTDATADIR/BUSCO"
mkdir -p "$OUTDATADIR/BUSCO"

# Get current directory, as it needs to be changed to control busco output
owd=$(pwd)
cd "${OUTDATADIR}"

export PYTHONPATH=$PYTHONPATH:"/apps/x86_64/busco/busco/build/lib/"
echo "${PATH//:/$'\n'}"

# Run busco on the prokka output using the database provided by command line ($2)
#run_BUSCO.py -i "${OUTDATADIR}/prokka/${sample_name}_PROKKA.faa" -o "${sample_name}" -l "${buscoDB}" -m prot -c "${procs}"

# Temporary fix while SCICOMP fixes normal call
singularity exec -B /scicomp:/scicomp /apps/x86_64/singularity-imgs/busco/3.0.2/busco.img run_BUSCO.py -i "${OUTDATADIR}/prokka/${sample_name}_PROKKA.faa" -o "${sample_name}" -l "${buscoDB}" -m prot -c "${procs}"

# Moves output files to proper location and removes temp files
mv -f "${OUTDATADIR}/run_${sample_name}/"* "${OUTDATADIR}/BUSCO"
rm -r "${OUTDATADIR}/run_${sample_name}"
rm -r "${OUTDATADIR}/tmp"

# returns current directory to original location
cd "${owd}"

# Unloads python 3.6.1 (and loads python 3.5.2 back in)
ml -Python/3.6.1 -busco/3.0.1

#Script exited gracefully (unless something else inside failed)
exit 0
