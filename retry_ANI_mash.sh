#!/bin/sh -l

#$ -o run_ANI.out
#$ -e run_ANI.err
#$ -N run_ANI
#$ -cwd
#$ -q short.q

#
# Description: Script to reattempt to make mash tree for closest samples to isolate
#
# Usage: ./retry_ANI_mash.sh -n sample_name -g genus -s species -p run_ID [-a list_samples_to_include(optional)] [-c path_to_config_file] [-a path_to_list_of_additional_samples]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/ANI/
#
# Modules required: mashtree/0.29
#
# v1.0.1 (08/24/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml mashtree/0.29


#  Function to print out help blurb
show_help () {
	echo "./retry_ANI_mash.sh -n sample_name -g genus -s	species -p run_ID [-a list_samples_to_include(optional)] [-c path_to_config_file] [-a path_to_list_of_additional_samples]"
}

# Parse command line options
options_found=0
while getopts ":h?c:p:n:g:s:a:" option; do
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
		g)
			echo "Option -g triggered, argument = ${OPTARG}"
			genus_in=${OPTARG};;
		s)
			echo "Option -s triggered, argument = ${OPTARG}"
			species=${OPTARG};;
		a)
			echo "Option -a triggered, argument = ${OPTARG}"
			others="true"
			additional_samples=${OPTARG};;
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
elif [[ -z "${species}" ]]; then
	echo "No species supplied to retry_ANI_best_hit.sh, exiting"
	exit 35
elif [[ -z "${genus}" ]]; then
	echo "No genus supplied to retry_ANI_best_hit.sh, exiting"
	exit 36
elif [[ "${others}" == "true" ]]; then
	if [[ ! -f "${additional_samples}" ]]; then
		echo "Additional sample file list does not exist...continuing without extra samples"
	fi
fi


# Checks for proper argumentation
if [ ! -s "${local_DBs}/aniDB/${2,}" ]; then
	echo "The genus does not exist in the ANI database. This will be noted and the curator of the database will be notified. However, since nothing can be done at the moment....exiting"
	# Create a dummy folder to put non-results into (if it doesnt exist
	if [ ! -d "${processed}/${4}/${sample_name}/ANI" ]; then  #create outdir if absent
		echo "${processed}/${4}/${sample_name}/ANI"
		mkdir -p "${processed}/${4}/${sample_name}/ANI"
	fi
	# Write non-results to a file in ANI folder
	echo "No matching ANI database found for ${genus_in}(genus)" >> "${processed}/${4}/${sample_name}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${genus_in}).txt"
	# Add genus to list to download and to database
	global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
	echo "ANI: ${genus_in} - Found as ${sample_name} on ${global_time}" >> "${shareScript}/maintenance_To_Do.txt"
	exit 1
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "Started ANI at ${start_time}"

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${processed}/${project}/${sample_name}"

# Checks to see if an ANI folder already exists and creates it if not
if [ ! -d "${OUTDATADIR}/ANI" ]; then
	echo "No ${OUTDATADIR}/ANI directory, exiting"
	exit
fi

# Gets persons name to use as email during entrez request to identify best matching sample
me=$(whoami)


#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
if [[ -s "${OUTDATADIR}/ANI/aniM/ANIm_percentage_identity.tab" ]]; then
	firstline=$(head -n 1 "${OUTDATADIR}/ANI/aniM/ANIm_percentage_identity.tab")
else
	echo "No ${OUTDATADIR}/ANI/aniM/ANIm_percentage_identity.tab file, exiting"
	exit 1
fi

#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line; do
#	echo "!-${line}"
	if [[ ${line:0:7} = "sample_" ]]; then
		sampleline=${line}
#		echo "found it-"$sampleline
		break
	fi
done < "${OUTDATADIR}/ANI/aniM/ANIm_percentage_identity.tab"


#Arrays to read sample names and the %ids for the query sample against those other samples
IFS="	" read -r -a samples <<< "${firstline}"
IFS="	" read -r -a percents <<< "${sampleline}"

#How many samples were compared
n=${#samples[@]}

#Extracts all %id against the query sample (excluding itself) and writes them to file
if [[ ! -d "${OUTDATADIR}/ANI/localANIDB" ]]; then
	mkdir "${OUTDATADIR}/ANI/localANIDB"
	for (( i=0; i<n; i++ ));
	do
		temp_ref=$(find ${local_DBs}/aniDB/${genus_in,,} -maxdepth 1 -type f -name "*${samples[i]}.fna.gz")
		echo "Trying to copy ${temp_ref} --- *${samples[i]}.fna.gz"
		if [[ -f ${temp_ref} ]]; then
			cp "${temp_ref}" "${OUTDATADIR}/ANI/localANIDB"
		else
			echo "Could not find ${temp_ref} (${samples[i]}.fna.gz)"
		fi
	done
	gunzip "${OUTDATADIR}/ANI/localANIDB/"*
	for f in ${OUTDATADIR}/ANI/localANIDB/*; do
		if [[ "${f}" == *".fasta" ]]; then
			mv $f `basename $f .fasta`.fna
		fi
	done
else
	echo "Already/still has its localANIDB folder"
fi

#rename 's/.fna$/.fasta/' ${OUTDATADIR}/ANI/localANIDB/*.fna
for foo in ${OUTDATADIR}/ANI/localANIDB/*.fna; do
	#echo "Moving $foo to `basename $foo .fna`.fasta"
	mv $foo ${OUTDATADIR}/ANI/localANIDB/`basename $foo .fna`.fasta;
done

mashtree --numcpus ${procs} *.fasta --tempdir ${OUTDATADIR}/ANI/temp > ${OUTDATADIR}/ANI/"${genus_in}_and_${sample_name}_mashtree.dnd";

#rename 's/.fasta$/.fna/' ${OUTDATADIR}/ANI/localANIDB/*.fasta
for foo in ${OUTDATADIR}/ANI/localANIDB/*.fasta; do mv $foo ${OUTDATADIR}/ANI/localANIDB/`basename $foo .fasta`.fna; done

rm -r ${OUTDATADIR}/ANI/temp

end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "ENDed ANI at ${end_time}"

#Script exited gracefully (unless something else inside failed)
exit 0
