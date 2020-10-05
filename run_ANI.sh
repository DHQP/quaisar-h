#!/bin/sh -l

#$ -o run_ANI.out
#$ -e run_ANI.err
#$ -N run_ANI
#$ -cwd
#$ -q short.q

#
# Description: Script to calculate the average nucleotide identity of a sample
#
# Usage: ./run_ANI.sh -n sample_name -g genus -s species	-p run_ID [-c path_to_config_file] [-a path_to_list_of_additional_samples]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/ANI/
#
# Modules required: Python3/3.5.2, pyani/0.2.7, mashtree/0.29
#
# V1.0.3 (09/30/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml Python3/3.5.2 pyani/0.2.7 mashtree/0.29


#  Function to print out help blurb
show_help () {
	echo "./run_ANI.sh ./run_ANI.sh -n sample_name -g genus -s species	-p run_ID [-c path_to_config_file] [-a path_to_list_of_additional_samples]"
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
elif [[ -z "${genus_in}" ]]; then
	echo "No genus supplied to retry_ANI_best_hit.sh, exiting"
	exit 36
elif [[ "${others}" == "true" ]]; then
	if [[ ! -f "${additional_samples}" ]]; then
		echo "Additional sample file list does not exist...continuing without extra samples"
	fi
fi



# Checks for proper argumentation
if [ ! -s "${local_DBs}/aniDB/${genus_in,,}" ]; then
	echo "The genus does not exist in the ANI database. This will be noted and the curator of the database will be notified. However, since nothing can be done at the moment....exiting"
	# Create a dummy folder to put non-results into (if it doesnt exist
	if [ ! -d "${processed}/${project}/${sample_name}/ANI" ]; then  #create outdir if absent
		echo "${processed}/${project}/${sample_name}/ANI"
		mkdir -p "${processed}/${project}/${sample_name}/ANI"
	fi
	# Write non-results to a file in ANI folder
	echo "No matching ANI database found for ${genus_in}(genus)" >> "${processed}/${project}/${sample_name}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${genus_in}).txt"
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
if [ ! -d "$OUTDATADIR/ANI" ]; then
	echo "Creating $OUTDATADIR/ANI"
	mkdir -p "$OUTDATADIR/ANI"
fi

# Checks to see if the local DB ANI folder already exists and creates it if not. This is used to store a local copy of all samples in DB to be compared to (as to not disturb the source DB)
if [ ! -d "$OUTDATADIR/ANI/localANIDB" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR/ANI/localANIDB"
	mkdir -p "$OUTDATADIR/ANI/localANIDB"
else
	rm -r "$OUTDATADIR/ANI/localANIDB"
	mkdir -p "$OUTDATADIR/ANI/localANIDB"
fi

# Checks to see if the local DB ANI folder already exists and creates it if not. This is used to store a local copy of all samples in DB to be compared to (as to not disturb the source DB)
if [ ! -d "$OUTDATADIR/ANI/temp" ]; then  #create outdir if absent
	echo "Creating $OUTDATADIR/ANI/temp"
	mkdir -p "$OUTDATADIR/ANI/temp"
else
	rm -r "$OUTDATADIR/ANI/temp"
	mkdir -p "$OUTDATADIR/ANI/temp"
fi

# Gets persons name to use as email during entrez request to identify best matching sample
me=$(whoami)
# Sets the genus as the database that was passed in (The $2 seemed to be getting changed somewhere, so I just set it as a local variable)

#Creates a local copy of the database folder
echo "trying to copy ${local_DBs}/aniDB/${genus_in,}/"
cp "${local_DBs}/aniDB/${genus_in,}/"*".fna.gz" "${OUTDATADIR}/ANI/localANIDB/"
gunzip ${OUTDATADIR}/ANI/localANIDB/*.gz

#Copies the samples assembly contigs to the local ANI db folder
cp "${OUTDATADIR}/Assembly/${sample_name}_scaffolds_trimmed.fasta" "${OUTDATADIR}/ANI/localANIDB/sample_${genus_in}_${species}.fasta"


# Add in all other assemblies to compare using list provided as argument
if [[ "${others}" = "true" ]]; then
	if [[ -f "${5}" ]]; then
		while IFS=read -r var  || [ -n "$var" ]; do
			sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
			project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
			echo "${sample_name}, ${project}"
			if [[ "${project}" == "${project}" ]] && [[ "${sample_name}" == "${sample_name}" ]]; then
				echo "Already in there as ref sample"
			else
				cp "${processed}/${project}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" "${OUTDATADIR}/ANI/localANIDB"
			fi
		done < ${5}
	else
		echo "List file: ${5}, does not exist. Analysis will continue with only database samples"
	fi
else
	echo "Analysis will be completed using only database isolates"
fi

#Renames all files in the localANIDB folder by changing extension from fna to fasta (which pyani needs)
for file in ${OUTDATADIR}/ANI/localANIDB/*.fna;
do
	fasta_name=$(basename "${file}" .fna)".fasta"
	mv "${file}" "${OUTDATADIR}/ANI/localANIDB/${fasta_name}"
done

# Mashtree trimming to reduce run time for ANI
owd=$(pwd)
cd ${OUTDATADIR}/ANI/localANIDB/
mashtree --numcpus ${procs} *.fasta --tempdir ${OUTDATADIR}/ANI/temp > ${OUTDATADIR}/ANI/"${genus_in}_and_${sample_name}_mashtree.dnd";

# Get total number of isolates compared in tree
sample_count=$(find ${OUTDATADIR}/ANI/localANIDB/ -type f | wc -l)
# Must remove sample of interest
sample_count=$(( sample_count - 1 ))
# Check if sample count is greater than the max samples for tree size, if so then reduce tree size to max closest samples balanced around submitted isolate
if [[ ${sample_count} -gt ${max_ani_samples} ]]; then
	sleep 2
	tree=$(head -n 1 "${OUTDATADIR}/ANI/${genus_in}_and_${sample_name}_mashtree.dnd")
	echo $tree
	tree=$(echo "${tree}" | tr -d '()')
	echo $tree
	IFS=',' read -r -a samples <<< "${tree}"
	counter=0
	half_max=$(( (max_ani_samples+1) / 2 ))
	echo "Halfsies = ${half_max}"
	for sample in ${samples[@]};
	do
		counter=$(( counter + 1 ))
		filename=$(echo ${sample} | cut -d':' -f1)
		echo "${filename}"
		filename="${filename}.fasta"
		if [[ "${filename}" == "sample_${genus_in}_${species}.fasta" ]]; then
			match=${counter}
			#echo "Match @ ${counter} and half=${half_max}"
			if [[ ${match} -le ${half_max} ]]; then
				#echo "LE"
				samples_trimmed=${samples[@]:0:$(( max_ani_samples + 1 ))}
			elif [[ ${match} -ge $(( sample_count - half_max)) ]]; then
				#echo "GE"
				samples_trimmed=${samples[@]:$(( max_ani_samples * -1 - 1 ))}
			else
				#echo "MID - $(( match - half_max )) to $(( counter + half_max + 1))"
				samples_trimmed=${samples[@]:$(( match - half_max )):${max_ani_samples}}
			fi
				#echo "${#samples_trimmed[@]}-${samples_trimmed[@]}"
				break
		fi
				#echo ${filename}
	done
	mkdir "${OUTDATADIR}/ANI/localANIDB_trimmed"
	for sample in ${samples_trimmed[@]};
	do
		filename=$(echo ${sample} | cut -d':' -f1)
		filename="${filename}.fasta"
		echo "Moving ${filename}"
		cp ${OUTDATADIR}/ANI/localANIDB/${filename} ${OUTDATADIR}/ANI/localANIDB_trimmed/
	done
	if [[ -d "${OUTDATADIR}/ANI/localANIDB_full" ]]; then
		rm -r "${OUTDATADIR}/ANI/localANIDB_full"
	fi
	mv "${OUTDATADIR}/ANI/localANIDB" "${OUTDATADIR}/ANI/localANIDB_full"
	mv "${OUTDATADIR}/ANI/localANIDB_trimmed" "${OUTDATADIR}/ANI/localANIDB"
# Continue without reducing the tree, as there are not enough samples to require reduction
else
	echo "Sample count below limit, not trimming ANI database"
fi

cd ${owd}
# Resume normal ANI analysis after mashtree reduction

# Checks for a previous copy of the aniM folder, removes it if found
if [ -d "${OUTDATADIR}/ANI/aniM" ]; then
	echo "Removing old ANIm results in ${OUTDATADIR}/ANI/aniM"
	rm -r "${OUTDATADIR}/ANI/aniM"
fi

if [[ -d "${OUTDATADIR}/ANI/localANIDB_full" ]]; then
	rm -r "${OUTDATADIR}/ANI/localANIDB_full"
fi

if [[ -d "${OUTDATADIR}/ANI/temp" ]]; then
	rm -r "${OUTDATADIR}/ANI/temp"
fi

#Calls pyani on local db folder
python -V
#python "${shareScript}/pyani/average_nucleotide_identity.py" -i "${OUTDATADIR}/ANI/localANIDB" -o "${OUTDATADIR}/ANI/aniM" --write_excel
average_nucleotide_identity.py -i "${OUTDATADIR}/ANI/localANIDB" -o "${OUTDATADIR}/ANI/aniM" --write_excel

#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line || [ -n "$line" ]; do
#	echo "!-${line}"
	if [[ ${line:0:7} = "sample_" ]]; then
		sampleline=${line}
#		echo "found it-"$sampleline
		break
	fi
done < "${OUTDATADIR}/ANI/aniM/ANIm_percentage_identity.tab"

#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line || [ -n "$line" ]; do
#	echo "!-${line}"
	if [[ ${line:0:6} = "sample" ]]; then
		sample_coverage_line=${line}
#		echo "found it-"$sample_identity_line
		break
	fi
done < "${OUTDATADIR}/ANI/aniM/ANIm_alignment_coverage.tab"

#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
if [[ -s "${OUTDATADIR}/ANI/aniM/ANIm_percentage_identity.tab" ]]; then
	firstline=$(head -n 1 "${OUTDATADIR}/ANI/aniM/ANIm_percentage_identity.tab")
else
	echo "No "${OUTDATADIR}/ANI/aniM/ANIm_percentage_identity.tab" file, exiting"
	exit 1
fi

#Arrays to read sample names and the %ids for the query sample against those other samples
IFS="	" read -r -a samples <<< "${firstline}"
IFS="	" read -r -a percents <<< "${sampleline}"
IFS="	" read -r -a coverages <<< "${sample_coverage_line}"

#How many samples were compared
n=${#samples[@]}

#Extracts all %id against the query sample (excluding itself) and writes them to file
for (( i=0; i<n; i++ ));
do
#	echo ${i}-${samples[i]}
	if [[ ${samples[i]:0:7} = "sample_" ]];
	then
#		echo "Skipping ${i}"
		continue
	fi
	definition=$(head -1 "${OUTDATADIR}/ANI/localANIDB/${samples[i]}.fasta")
	# Prints all matching samples to file (Except the self comparison) by line as percent_match  sample_name  fasta_header
	echo "${percents[i+1]}	${coverages[i+1]}	${samples[i]}	${definition}" >> "${OUTDATADIR}/ANI/best_hits.txt"
done

#Sorts the list in the file based on %id (best to worst)
sort -nr -t' ' -k1 -o "${OUTDATADIR}/ANI/best_hits_ordered.txt" "${OUTDATADIR}/ANI/best_hits.txt"
#Extracts the first line of the file (best hit)
best=$(head -n 1 "${OUTDATADIR}/ANI/best_hits_ordered.txt")
#Creates an array from the best hit
IFS='	' read -r -a def_array <<< "${best}"
#echo -${def_array[@]}+
#Captures the assembly file name that the best hit came from

echo "A-${best}"
echo -e "B1-${def_array[0]}\nB2-${def_array[1]}\nB3-${def_array[2]}\nB4${def_array[3]}"
best_file=${def_array[2]}
echo "C-${best_file}"

#Formats the %id to standard percentage (xx.xx%)
best_percent=$(awk -v per="${def_array[0]}" 'BEGIN{printf "%.2f", per * 100}')
best_coverage=$(awk -v per="${def_array[1]}" 'BEGIN{printf "%.2f", per * 100}')
#echo "${best_file}"
# If the best match comes from the additional file, extract the taxonomy from that file
# if [[ "${best_file}" = *"_scaffolds_trimmed" ]]; then
# 	best_outbreak_match=$(echo "${best_file}" | rev | cut -d'_' -f3- | rev)
# 	while IFS= read -r var || [ -n "$var" ]; do
# 		sample_name=$(echo "${var}" | cut -d'/' -f2 | tr -d '[:space:]')
# 		if [[ "${sample_name}" = "${best_outbreak_match}" ]]; then
# 			project=$(echo "${var}" | cut -d'/' -f1 | tr -d '[:space:]')
# 			while IFS= read -r pstats_line || [ -n "$pstats_line" ]; do
# 					tool=$(echo "${pstats_line}" | cut -d':' -f1 | tr -s " ")
# 					#echo ":${tool}:"
# 					if [[ "${tool}" = "weighted Classify " ]]; then
# 						best_organism_guess=$(echo "${pstats_line}" | cut -d':' -f3 | cut -d' ' -f3,4)
# 						break 2
# 					fi
# 				done < ${processed}/${project}/${sample_name}/${sample_name}_pipeline_stats.txt
# 		fi
# 	done < ${5}
# # if the best hit cmoes from the aniDB then lookup the taxonomy on ncbi
# else
# 	#Extracts the accession number from the definition line
# 	accession=$(echo "${def_array[2]}" | cut -d' ' -f1  | cut -d'>' -f2)
# 	#Looks up the NCBI genus species from the accession number
# 	if [[ "${accession}" == "No_Accession_Number" ]]; then
# 		best_organism_guess="${def_array[3]} ${def_array[4]}"
# 	else
# 		attempts=0
# 		while [[ ${attempts} -lt 25 ]]; do
# 			best_organism_guess=$(python "${shareScript}/entrez_get_taxon_from_accession.py" -a "${accession}" -e "${me}")
# 			if [[ ! -z ${best_organism_guess} ]]; then
# 				break
# 			else
# 				attempts=$(( attempts + 1 ))
# 			fi
# 		done
# 	fi
# fi

best_genus=$(echo "${best_file}" | cut -d'_' -f1)
best_species=$(echo "${best_file}" | cut -d'_' -f2)
best_organism_guess="${best_genus} ${best_species}"

echo "G-${best_genus}:S-${best_species}"

# Uncomment this if you want to restrict ID to only genus species, without more resolute definition
#best_organism_guess_arr=($best_organism_guess})
#best_organism_guess="${best_organism_guess_arr[@]:0:2}"

#Creates a line at the top of the file to show the best match in an easily readable format that matches the style on the MMB_Seq log
echo -e "${best_percent}%ID-${best_coverage}%COV-${best_organism_guess}(${best_file}.fna)\\n$(cat "${OUTDATADIR}/ANI/best_hits_ordered.txt")" > "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${sample_name}_vs_${genus_in}).txt"

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

#Script exited gracefully (unless something else inside failed)
exit 0
