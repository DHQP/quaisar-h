#!/bin/sh -l

#$ -o do_busco.out
#$ -e do_busco.err
#$ -N do_busco
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: A script that takes a sample and compares it to a busco database to discover number of similar genes (% conserved proteins) from prokka output
#
# Usage ./do_busco.sh sample_name DB(to search against) run_ID
#
# Output location: default_config.sh_output_location/run_ID/sample_name/BUSCO/
#
# Modules required: busco/3.0.1, Python3/3.5.4 (whatever version used to install it, must have pipebricks)
#
# v1.0.2 (10/18/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml busco/3.0.1 Python3/3.6.1

python3 -V

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to $0, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./do busco.sh   sample_name   database_name   run_ID"
	echo "Output is saved to ${processed}/miseq_run_ID/sample_name/busco"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty database name supplied to$0, exiting"
	exit 1
elif [ ! -s "${local_DBs}/BUSCO/${2,}" ]; then
	echo "Exists? - ${local_DBs}/BUSCO/${2,}"
	echo "The taxon does not exist in the BUSCO database. This will be noted and the curator of the database will be notified. However, since nothing can be done at the moment....exiting"
	# Create a dummy folder to put non-results into (if it doesnt exist)
	if [ ! -d "${processed}/${4}/${1}/BUSCO" ]; then  #create outdir if absent
		echo "${processed}/${4}/${1}/BUSCO"
		mkdir -p "${processed}/${4}/${1}/BUSCO"
	fi
	# Write non-results to a file in busco folder
	echo "No matching DB database found for ${2}" >> "${processed}/${3}/${1}/BUSCO/summary_${1}.txt"
	# Add genus to list to download and to database
	global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
	echo "BUSCO: ${2} - Found as ${1} on ${global_time}" >> "${shareScript}/maintenance_To_Do.txt"
	exit 1
elif [ -z "$3" ]; then
	echo "Empty project id supplied to do_busco.sh, exiting"
	exit 1
fi

# Sets output folder to parent isolate folder
OUTDATADIR="${processed}/${3}/${1}"
# Sets buscoDB to the folder matching the 2nd argument (dbs vary on taxonomic level)
buscoDB="${local_DBs}/BUSCO/${2,}"

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
run_BUSCO.py -i "${OUTDATADIR}/prokka/${1}_PROKKA.faa" -o "${1}" -l "${buscoDB}" -m prot -c "${procs}"


# Moves output files to proper location and removes temp files
mv -f "${OUTDATADIR}/run_${1}/"* "${OUTDATADIR}/BUSCO"
rm -r "${OUTDATADIR}/run_${1}"
rm -r "${OUTDATADIR}/tmp"

# returns current directory to original location
cd "${owd}"

# Unloads python 3.6.1 (and loads python 3.5.2 back in)
#. "${mod_changers}/unload_python_3.6.1.sh"

ml -Python/3.6.1 -busco/3.0.1

#Script exited gracefully (unless something else inside failed)
exit 0
