#!/bin/sh -l

#$ -o run_Assembly_Quality_Check.out
#$ -e run_Assembly_Quality_Check.err
#$ -N run_Assembly_Quality_Check
#$ -cwd
#$ -q short.q

#
# Description: Checks the Assembly quality  using Toms tool and QUAST and comparing the output of both
# 	Important stats are # of contigs, assembly length, n%0 and G/C content
#
# Usage: ./run_Assembly_Quality_Check.sh   -n sample_name   -p run_ID [-l] [-c path_to_config_file]
# 	Optional l flag is to run it only on the plasmid assembly, assuming it is in the default location of the config file and unicycled
#
# Output location: default_config.sh_output_location/run_ID/sample_name/Assembly_Stats(_plasFlow)
#
# Modules required: quast/4.3, python2/2.7.15(loaded by quast/4.3)
#
# v1.0 (10/3/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml quast/4.3 Python2/2.7.15


#  Function to print out help blurb
show_help () {
	echo "Usage: ./run_Assembly_Quality_Check.sh -n sample_name -p run_ID [-l] [-c path_to_config_file]"
}

do_plasFlow_only="false"

# Parse command line options
options_found=0
while getopts ":h?c:p:n:l" option; do
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
		l)
			echo "Option -l triggered, plasmid mode activated"
			do_plasFlow_only="true";;
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

if [[ -z "${project}" ]]; then
	echo "No Project/Run_ID supplied to retry_ANI_best_hit.sh, exiting"
	exit 33
elif [[ -z "${sample_name}" ]]; then
	echo "No sample name supplied to retry_ANI_best_hit.sh, exiting"
	exit 34
fi


# Sets the parent output folder as the sample name folder in the processed samples folder in MMB_Data
OUTDATADIR="${processed}/${project}/${sample_name}"


echo "Checking Assembly QC with QUAST"
# Run QUAST
# Save current directory and move to output directory because it doesnt know how to redirect output
owd="$(pwd)"
if [[ "${do_plasFlow_only}" != "true" ]]; then
	# Checks for output folder existence and creates creates if not
	if [ ! -d "$OUTDATADIR/Assembly_Stats" ]; then
		echo "Creating $OUTDATADIR/Assembly_Stats"
		mkdir -p "$OUTDATADIR/Assembly_Stats"
	fi
	cd "${OUTDATADIR}/Assembly_Stats"
	# Call QUAST
	python2 "/apps/x86_64/quast/quast-4.3/quast.py" -o "${OUTDATADIR}/Assembly_Stats" "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta"
	mv "${OUTDATADIR}/Assembly_Stats/report.txt" "${OUTDATADIR}/Assembly_Stats/${sample_name}_report.txt"
	mv "${OUTDATADIR}/Assembly_Stats/report.tsv" "${OUTDATADIR}/Assembly_Stats/${sample_name}_report.tsv"
fi
if [[ -s "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta" ]]; then
	if [ ! -d "$OUTDATADIR/Assembly_Stats_plasFlow" ]; then
		echo "Creating $OUTDATADIR/Assembly_Stats_plasFlow"
 		mkdir -p "$OUTDATADIR/Assembly_Stats_plasFlow"
 	fi
 	python2 "/apps/x86_64/quast/quast-4.3/quast.py" -o "${OUTDATADIR}/Assembly_Stats_plasFlow" "${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta"
 	mv "${OUTDATADIR}/Assembly_Stats_plasFlow/report.txt" "${OUTDATADIR}/Assembly_Stats_plasFlow/${sample_name}_report.txt"
 	mv "${OUTDATADIR}/Assembly_Stats_plasFlow/report.tsv" "${OUTDATADIR}/Assembly_Stats_plasFlow/${sample_name}_report.tsv"
 else
	echo "No plasFlow assembly (${OUTDATADIR}/plasFlow/Unicycler_assemblies/${sample_name}_uni_assembly/${sample_name}_plasmid_assembly_trimmed.fasta)"
fi

# Return to original directory
cd "${owd}"

ml -quast/4.3 -Python/2.7.15

#Show that the script seemingly completed successfully
exit 0
