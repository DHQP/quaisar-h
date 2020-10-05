#!/bin/sh -l

#$ -o append_tax.out
#$ -e append_tax.err
#$ -N rbl
#$ -cwd
#$ -q short.q

#
# Usage ./append_taxonomy_to_ncbi_assembly_filenames.sh path_to_list
# The folder needs to have ncbi source assembly files gzipped
#

#  Function to print out help blurb
show_help () {
	echo "Usage is ./append_taxonomy_to_ncbi_assembly_filenames.sh -i path_to_folder"
	echo "Output is done in folder"
}

# Parse command line options
options_found=0
while getopts ":h?i:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		i)
			echo "Option -i triggered, argument = ${OPTARG}"
			folder=${OPTARG};;
		:)
			echo "Option -${OPTARG} requires as argument";;
		h)
			show_help
			exit 0
			;;
	esac
done

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	show_help
	exit 1
elif [[ ! -d "${folder}" ]]; then
	echo "No folder provided to work on...exiting"
	show_help
	exit 0
fi

# Loop through and act on each sample name in the passed/provided list

for i in ${folder}/*.gz; do
	old_name=$(basename ${i} | rev | cut -d'.' -f2- | rev)
	new_name=$(echo ${old_name} | tr -d '[],')
	dir_name=$(dirname ${i})
	gunzip ${i}
	tax_genus=$(head -n1 "${dir_name}/${old_name}" | cut -d' ' -f2 | tr -d '[],')
	tax_species=$(head -n1 "${dir_name}/${old_name}" | cut -d' ' -f3 | tr -d '[],')
	echo "Taxes: ${tax_genus}:${tax_species}"
	mv ${dir_name}/${old_name} ${dir_name}/${tax_genus}_${tax_species}_${new_name}
	gzip ${dir_name}/${tax_genus}_${tax_species}_${new_name}
done
