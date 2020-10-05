#!/bin/sh -l

#$ -o run_ANI.out
#$ -e run_ANI.err
#$ -N run_ANI
#$ -cwd
#$ -q short.q

#
# Description: Will attempt to recalculate the best hit from a previously completed ANI analysis
#
# Usage: ./retry_ANI_best_hit.sh -n sample_name -g genus -s	species -p run_ID [-a list_samples_to_include(optional)] [-c path_to_config_file] [-a path_to_list_of_additional_samples]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/ANI/
#
# Modules required: Python/3.5.2
#
# v1.0.1 (08/24/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml Python3/3.5.2


#  Function to print out help blurb
show_help () {
	echo "./retry_ANI_best_hit.sh -n sample_name -g genus -s	species -p run_ID [-a list_samples_to_include(optional)] [-c path_to_config_file] [-a path_to_list_of_additional_samples]"
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
		echo "Trying to copy ${temp_ref} --- ${samples[i]}.fna.gz"
		if [[ -f ${temp_ref} ]]; then
			cp "${temp_ref}" "${OUTDATADIR}/ANI/localANIDB"
		else
			echo "Could not copy ${temp_ref} (*${samples[i]}.fna.gz)"
		fi
	done
	gunzip "${OUTDATADIR}/ANI/localANIDB/"*
else
	echo "Already/still has its localANIDB folder"
fi

temp_ref=""

#ls -l "${OUTDATADIR}/ANI/localANIDB"

for (( i=0; i<n; i++ ));
do
#	echo ${i}-${samples[i]}
	if [[ ${samples[i]:0:7} = "sample_" ]]; then
#		echo "Skipping ${i}"
		continue
	else
		temp_ref=$(find "${OUTDATADIR}/ANI/localANIDB" -type f -name "*${samples[i]}.fna")
		if [[ -f ${temp_ref} ]]; then
			definition=$(head -1 "${temp_ref}")
			# Prints all matching samples to file (Except the self comparison) by line as percent_match  sample_name  fasta_header
			echo "${percents[i+1]} ${samples[i]} ${definition}" >> "${OUTDATADIR}/ANI/best_hits.txt"
		else
			echo "Could not find ${temp_ref} (*${samples[i]}.fna)"
			echo "${percents[i+1]} ${samples[i]} NO_FILE_FOUND-NO_ACCESSION" >> "${OUTDATADIR}/ANI/best_hits.txt"
		fi
	fi
done

#Sorts the list in the file based on %id (best to worst)
sort -nr -t' ' -k1 -o "${OUTDATADIR}/ANI/best_hits_ordered.txt" "${OUTDATADIR}/ANI/best_hits.txt"
#Extracts the first line of the file (best hit)
best=$(head -n 1 "${OUTDATADIR}/ANI/best_hits_ordered.txt")
#Creates an array from the best hit
IFS=' ' read -r -a def_array <<< "${best}"
#echo -${def_array[@]}+
#Captures the assembly file name that the best hit came from
best_file=${def_array[1]}
#Formats the %id to standard percentage (xx.xx%)
best_percent=$(awk -v per="${def_array[0]}" 'BEGIN{printf "%.2f", per * 100}')
#echo "${best_file}"
# If the best match comes from the additional file, extract the taxonomy from that file
if [[ "${best_file}" = *"_scaffolds_trimmed" ]]; then
	best_outbreak_match=$(echo "${best_file}" | rev | cut -d'_' -f3- | rev)
	while IFS= read -r var || [ -n "$var" ]; do
		sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
		if [[ "${sample_name}" = "${best_outbreak_match}" ]]; then
			project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
			while IFS= read -r pstats_line; do
					tool=$(echo "${pstats_line}" | cut -d':' -f1 | tr -s " ")
					#echo ":${tool}:"
					if [[ "${tool}" = "weighted Classify " ]]; then
						best_organism_guess=$(echo "${pstats_line}" | cut -d':' -f3 | cut -d' ' -f3,4)
						break 2
					fi
				done < ${processed}/${project}/${sample_name}/${sample_name}_pipeline_stats.txt
		fi
	done < ${5}
# if the best hit cmoes from the aniDB then lookup the taxonomy on ncbi
else
	#Extracts the accession number from the definition line
	accession=$(echo "${def_array[2]}" | cut -d' ' -f1  | cut -d'>' -f2)
	#Looks up the NCBI genus species from the accession number
	if [[ "${accession}" == "No_Accession_Number" ]]; then
		best_organism_guess="${def_array[3]} ${def_array[4]}"
	else
		ml Entrez/latest
		attempts=0
		while [[ ${attempts} -lt 5 ]]; do
			#echo "Trying to lookup - ${accession}"
			best_organism_guess=$(python3 "${shareScript}/entrez_get_taxon_from_accession.py" -a "${accession}" -e "${me}")
			if [[ ! -z ${best_organism_guess} ]]; then
				best_organism_guess=$(echo "${best_organism_guess}" | tr -d "[]")
				break
			else
				attempts=$(( attempts + 1 ))
			fi
		done
		ml - Entrez/latest
	fi
fi
# Uncomment this if you want to restrict ID to only genus species, without more resolute definition
#best_organism_guess_arr=($best_organism_guess})
#best_organism_guess="${best_organism_guess_arr[@]:0:2}"

#Creates a line at the top of the file to show the best match in an easily readable format that matches the style on the MMB_Seq log
echo -e "${best_percent}%-${best_organism_guess}(${best_file}.fna)\\n$(cat "${OUTDATADIR}/ANI/best_hits_ordered.txt")" > "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${genus_in}).txt"

"${shareScript}/determine_taxID.sh" -n "${sample_name}" -p "${project}" -c "${config}"

genus=$(tail -n3 ${OUTDATADIR}/${sample_name}.tax | head -n1 | cut -d'	' -f2)
species=$(tail -n2 ${OUTDATADIR}/${sample_name}.tax | head -n1 | cut -d'	' -f2)

echo -e "${genus^} ${species}"

#Removes the transient hit files
if [ -s "${OUTDATADIR}/ANI/best_hits.txt" ]; then
	rm "${OUTDATADIR}/ANI/best_hits.txt"
#	echo "1"
fi
if [ -s "${OUTDATADIR}/ANI/best_hits_ordered.txt" ]; then
	rm "${OUTDATADIR}/ANI/best_hits_ordered.txt"
#	echo "2"
fi

end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "ENDed ANI at ${end_time}"

ml -Python3/3.5.2

#Script exited gracefully (unless something else inside failed)
exit 0
