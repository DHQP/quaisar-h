#!/bin/sh -l

#$ -o append_tax.out
#$ -e append_tax.err
#$ -N rbl
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
. ./config.sh
#Import the module file that loads all necessary mods
#. "${mod_changers}/pipeline_mods"

#
# Usage ./append_taxonomy_to_ncbi_assembly_filenames.sh path_to_list
# The list needs to have project/old_name:new_name
#


# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./append_taxonomy_to_ncbi_assembly_filenames.sh path_to_folder_of_assemblies"
		exit 0
fi

# Loop through and act on each sample name in the passed/provided list

for i in ${1}/*.gz; do
	old_name=$(basename ${i} | rev | cut -d'.' -f2- | rev)
	new_name=$(echo ${old_name} | tr -d '[],')
	dir_name=$(dirname ${i})
	gunzip ${i}
	tax_genus=$(head -n1 "${dir_name}/${old_name}" | cut -d' ' -f2 | tr -d '[],')
	tax_species=$(head -n1 "${dir_name}/${old_name}" | cut -d' ' -f3 | tr -d '[],')
	echo "Taxes: ${tax_genus}:${tax_species}"
	mv ${dir_name}/${old_name} ${dir_name}/${tax_genus}_${tax_species}_${new_name}
	gzip ${dir_name}/${tax_genus}_${tax_species}_${new_name}
done < "${1}"
