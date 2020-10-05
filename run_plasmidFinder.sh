#!/bin/sh -l

#$ -o run_plasmidFinder.out
#$ -e run_plasmidFinder.err
#$ -N run_plasmidFinder
#$ -cwd
#$ -q short.q

#
# Description: Will attempt to find any plasmids in sample
#
# Usage ./run_plasmidFinder.sh -n sample_name -p run_ID -o output_folder [plasmidFinder|plasmidFinder_on_plasFlow] [-f force against all] [-c path_to_config_file]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/plasmidFinder(on_plasFlow)/
#
# Modules required: PlasmidFinder/1.3
#
# v1.0.2 (09/08/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml PlasmidFinder/1.3

#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_plasmidFinder.sh -n sample_name -p run_ID -o output_folder [plasmidFinder|plasmidFinder_on_plasFlow] [-f force against all] [-c path_to_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?c:p:n:o:f" option; do
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
		f)
			echo "Option -f triggered, argument = ${OPTARG}"
			force="true";;
		o)
			echo "Option -o triggered, argument = ${OPTARG}"
			output_dir=${OPTARG,,};;
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
elif [ -z "${output_dir}" ] || [[ "${output_dir}" != "plasmidfinder" ]] && [[ "${output_dir}" != "plasmidfinder_on_plasflow" ]]; then
	echo "Output folder must be plasmidFinder or plasmidFinder_on_plasFlow"
	exit 1
fi

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${processed}/${project}/${sample_name}"

# Create output directory
if [[ ! -d ${OUTDATADIR} ]]; then
	echo "Making ${OUTDATADIR}"
	mkdir ${OUTDATADIR}
fi


# Get proper input file based on output directory (whether it is full assembly or plasmid)
# if [[ "${output_dir}" == "plasmid_on_plasFlow" ]]; then
# 	inpath="plasmidAssembly/${sample_name}_plasmid_scaffolds_trimmed.fasta"
# el
if [[ "${output_dir}" == "plasmidfinder" ]]; then
	output_dir="plasmidFinder"
	inpath="Assembly/${sample_name}_scaffolds_trimmed.fasta"
elif [[ "${output_dir}" == "plasmidfinder_on_plasflow" ]]; then
	output_dir="plasmidFinder_on_plasFlow"
	if [[ -f "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/assembly.fasta" ]]; then
		mv "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/assembly.fasta" "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_assembly.fasta"
		mv "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/assembly.gfa" "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_assembly.gfa"
	fi
	if 	[[ -f "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta" ]]; then
		inpath="plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta"
	else
		if [[ ! -d "${OUTDATADIR}/plasFlow" ]]; then
			echo "plasFlow folder does not exist for ${project}/${sample_name}, it likely is not in the Enterobacteriaceae family"
		else
			echo "No ${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta"
		fi
		exit
	fi
else
	echo "Non standard output location, using full assembly to find plasmids"
	inpath="Assembly/${sample_name}_scaffolds_trimmed.fasta"
fi

#If force flag is set, then run it against all databases
if [[ "${force}" == "true" ]]; then
	echo "Checking against ALL plasmids, but unlikely to find anything"
	plasmidfinder -i ${OUTDATADIR}/${inpath} -o ${OUTDATADIR}/${output_dir} -k ${plasmidFinder_identity} -p enterobacteriaceae
	# Rename all files to include ID
	mv ${OUTDATADIR}/${output_dir}/Hit_in_genome_seq.fsa ${OUTDATADIR}/${output_dir}/${sample_name}_Hit_in_genome_seq_entero.fsa
	mv ${OUTDATADIR}/${output_dir}/Plasmid_seq.fsa ${OUTDATADIR}/${output_dir}/${sample_name}_Plasmid_seq_enetero.fsa
	mv ${OUTDATADIR}/${output_dir}/results.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_entero.txt
	mv ${OUTDATADIR}/${output_dir}/results_tab.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_tab_entero.txt
	mv ${OUTDATADIR}/${output_dir}/results_table.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_table_entero.txt
	plasmidfinder -i ${OUTDATADIR}/${inpath} -o ${OUTDATADIR}/${output_dir} -k ${plasmidFinder_identity} -p gram_positive
	# Rename all files to include ID
	mv ${OUTDATADIR}/${output_dir}/Hit_in_genome_seq.fsa ${OUTDATADIR}/${output_dir}/${sample_name}_Hit_in_genome_seq_gramp.fsa
	mv ${OUTDATADIR}/${output_dir}/Plasmid_seq.fsa ${OUTDATADIR}/${output_dir}/${sample_name}_Plasmid_seq_gramp.fsa
	mv ${OUTDATADIR}/${output_dir}/results.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_gramp.txt
	mv ${OUTDATADIR}/${output_dir}/results_tab.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_tab_gramp.txt
	mv ${OUTDATADIR}/${output_dir}/results_table.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_table_gramp.txt
	cat	${OUTDATADIR}/${output_dir}/${sample_name}_results_table_gramp.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_table_entero.txt > ${OUTDATADIR}/${output_dir}/${sample_name}_results_table_summary.txt
# Else, if the force flag is not set, then TRY to limit search to family (it will still check against all if it does not match the family)
else
	# Checks to see if a tax file is available to extract the family of the sample
	if [[ -f ${OUTDATADIR}/${sample_name}.tax ]]; then
		#Extracts the 6th line from the tax file containing all family information
		family=$(sed -n '6p' < ${OUTDATADIR}/${sample_name}.tax)
		genus=$(sed -n '7p' < ${OUTDATADIR}/${sample_name}.tax)
		#Extracts family name from line
		family=$(echo ${family}  | cut -d' ' -f4)
		genus=$(echo ${genus}  | cut -d' ' -f4)
		echo "${family}-${genus}"
		# If family is enterobacteriaceae, then run against that DB
		if [[ "${family,}" == "enterobacteriaceae" ]]; then
			echo "Checking against Enterobacteriaceae plasmids"
			plasmidfinder -i ${OUTDATADIR}/${inpath} -o ${OUTDATADIR}/${output_dir} -k ${plasmidFinder_identity} -p enterobacteriaceae
			# Rename all files to include ID
			mv ${OUTDATADIR}/${output_dir}/Hit_in_genome_seq.fsa ${OUTDATADIR}/${output_dir}/${sample_name}_Hit_in_genome_seq_entero.fsa
			mv ${OUTDATADIR}/${output_dir}/Plasmid_seq.fsa ${OUTDATADIR}/${output_dir}/${sample_name}_Plasmid_seq_enetero.fsa
			mv ${OUTDATADIR}/${output_dir}/results.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_entero.txt
			mv ${OUTDATADIR}/${output_dir}/results_tab.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_tab_entero.txt
			mv ${OUTDATADIR}/${output_dir}/results_table.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_table_summary.txt
		# If family is staph, strp, or enterococcus, then run against the gram positive database
		elif [[ "${genus,}" == "staphylococcus" ]] || [[ "${3,}" == "streptococcus" ]] || [[ "${3,}" == "enterococcus" ]]; then
			echo "Checking against Staph, Strep, and Enterococcus plasmids"
			plasmidfinder -i ${OUTDATADIR}/${inpath} -o ${OUTDATADIR}/${output_dir} -k ${plasmidFinder_identity} -p gram_positive
			# Rename all files to include ID
			mv ${OUTDATADIR}/${output_dir}/Hit_in_genome_seq.fsa ${OUTDATADIR}/${output_dir}/${sample_name}_Hit_in_genome_seq_gramp.fsa
			mv ${OUTDATADIR}/${output_dir}/Plasmid_seq.fsa ${OUTDATADIR}/${output_dir}/${sample_name}_Plasmid_seq_gramp.fsa
			mv ${OUTDATADIR}/${output_dir}/results.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_gramp.txt
			mv ${OUTDATADIR}/${output_dir}/results_tab.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_tab_gramp.txt
			mv ${OUTDATADIR}/${output_dir}/results_table.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_table_summary.txt
		# Family is not one that has been designated by the creators of plasmidFinder to work well, but still attempting to run against both databases
		else
			echo "Checking against ALL plasmids, but unlikely to find anything"
			plasmidfinder -i ${OUTDATADIR}/${inpath} -o ${OUTDATADIR}/${output_dir} -k ${plasmidFinder_identity} -p enterobacteriaceae
			# Rename all files to include ID
			mv ${OUTDATADIR}/${output_dir}/Hit_in_genome_seq.fsa ${OUTDATADIR}/${output_dir}/${sample_name}_Hit_in_genome_seq_entero.fsa
			mv ${OUTDATADIR}/${output_dir}/Plasmid_seq.fsa ${OUTDATADIR}/${output_dir}/${sample_name}_Plasmid_seq_enetero.fsa
			mv ${OUTDATADIR}/${output_dir}/results.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_entero.txt
			mv ${OUTDATADIR}/${output_dir}/results_tab.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_tab_entero.txt
			mv ${OUTDATADIR}/${output_dir}/results_table.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_table_entero.txt
			plasmidfinder -i ${OUTDATADIR}/${inpath} -o ${OUTDATADIR}/${output_dir} -k ${plasmidFinder_identity} -p gram_positive
			# Rename all files to include ID
			mv ${OUTDATADIR}/${output_dir}/Hit_in_genome_seq.fsa ${OUTDATADIR}/${output_dir}/${sample_name}_Hit_in_genome_seq_gramp.fsa
			mv ${OUTDATADIR}/${output_dir}/Plasmid_seq.fsa ${OUTDATADIR}/${output_dir}/${sample_name}_Plasmid_seq_gramp.fsa
			mv ${OUTDATADIR}/${output_dir}/results.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_gramp.txt
			mv ${OUTDATADIR}/${output_dir}/results_tab.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_tab_gramp.txt
			mv ${OUTDATADIR}/${output_dir}/results_table.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_table_gramp.txt
			cat	${OUTDATADIR}/${output_dir}/${sample_name}_results_table_gramp.txt ${OUTDATADIR}/${output_dir}/${sample_name}_results_table_entero.txt > ${OUTDATADIR}/${output_dir}/${sample_name}_results_table_summary.txt

		fi
	# No assembly file exists and cannot be used to determine family of sample
	else
		echo "Cant guess the genus of the sample, please try again with the force option or check the contents of the .tax file for complete taxonomic classification (${OUTDATADIR}/${sample_name}.tax)"
	fi

	ml -PlasmidFinder/1.3
fi
