#!/bin/sh -l

#$ -o run_ANI_sketch.out
#$ -e run_ANI_sketch.err
#$ -N run_ANI_sketch
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp config_template.sh config.sh
fi
. ./config.sh

#
# Description: Script to calculate the average nucleotide identity of a sample
#
# Usage: ./run_ANI.sh sample_name	run_ID
#
# Output location: default_config.sh_output_location/run_ID/sample_name/ANI/
#
# Modules required: Python3/3.5.2, pyani/0.2.7, Mash/2.0
#
# V1.0.3(04/29/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml Python3/3.5.2 pyani/0.2.7 Mash/2.0

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty sample name supplied to run_ANI_REFSEQ.sh, exiting"
	exit 1
# Gives the user a brief usage and help section if requested with the -h option argument
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_ANI_REFSEQ.sh sample_name run_ID"
	echo "Output is saved to in ${processed}/sample_name/ANI"
	exit 0
elif [ -z "${2}" ]; then
	echo "Empty miseq_run_ID name supplied to run_ANI_REFSEQ.sh. Fourth argument should be the run id. Exiting"
	exit 1
fi

start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "Started ANI at ${start_time}"

# Sets the output folder to the sample_name folder in processed samples
OUTDATADIR="${processed}/${2}/${1}"

# Checks to see if an ANI folder already exists and creates it if not
if [ ! -d "${OUTDATADIR}/ANI" ]; then
	echo "Creating ${OUTDATADIR}/ANI"
	mkdir -p "${OUTDATADIR}/ANI"
fi

# Checks to see if the local DB ANI folder already exists and creates it if not. This is used to store a local copy of all samples in DB to be compared to (as to not disturb the source DB)
if [ ! -d "${OUTDATADIR}/ANI/localANIDB_REFSEQ" ]; then  #create outdir if absent
	echo "Creating ${OUTDATADIR}/ANI/localANIDB_REFSEQ"
	mkdir -p "${OUTDATADIR}/ANI/localANIDB_REFSEQ"
else
	rm -r "${OUTDATADIR}/ANI/localANIDB_REFSEQ"
	mkdir -p "${OUTDATADIR}/ANI/localANIDB_REFSEQ"
fi

# Gets persons name to use as email during entrez request to identify best matching sample
me=$(whoami)

#Copies the samples assembly contigs to the local ANI db folder
cp "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta" "${OUTDATADIR}/ANI/localANIDB_REFSEQ/sample.fasta"

echo "RS-${REFSEQ} --- RSF-${REFSEQ_date}"

# Mashtree trimming to reduce run time for ANI
mash dist "${OUTDATADIR}/Assembly/${1}_scaffolds_trimmed.fasta" "${REFSEQ}" > "${OUTDATADIR}/ANI/${1}_${REFSEQ_date}_mash.dists"
sort -k3 -n -o "${OUTDATADIR}/ANI/${1}_${REFSEQ_date}_mash_sorted.dists" "${OUTDATADIR}/ANI/${1}_${REFSEQ_date}_mash.dists"
rm "${OUTDATADIR}/ANI/${1}_${REFSEQ_date}_mash.dists"

cutoff=$(head -n${max_ani_samples} "${OUTDATADIR}/ANI/${1}_${REFSEQ_date}_mash_sorted.dists" | tail -n1 | cut -d'	' -f3)

echo "Cutoff IS: ${cutoff}"

while IFS= read -r var; do
	echo "${var}"
	dist=$(echo ${var} | cut -d' ' -f3)
	kmers=$(echo ${var} | cut -d' ' -f5 | cut -d'/' -f1)
	comparison_files=$(ls -lh "${OUTDATADIR}/ANI/localANIDB_REFSEQ/" | wc -l)
	# Subtract one for size and one for sample.fasta
	comparison_files=$(( comparison_files - 2 ))
	echo "dist-${dist}"
	if (( $(echo "$dist <= $cutoff" | bc -l) )) && [ ${kmers} -gt 0 ] || [[ "${comparison_files}" -lt "${max_ani_samples}" ]]; then
		filename=$(echo ${var} | cut -d' ' -f2 | rev | cut -d'/' -f1 | rev | cut -d'_' -f3- | rev | cut -d'_' -f2,3,4 | rev)
		alpha=${filename:4:3}
		beta=${filename:7:3}
		charlie=${filename:10:3}
		echo "Copying - ${filename}"
		echo "Trying - wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${filename}_genomic.fna.gz -P ${OUTDATADIR}/ANI/localANIDB_REFSEQ"
		wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${filename}_genomic.fna.gz -P ${OUTDATADIR}/ANI/localANIDB_REFSEQ
		#curl --remote-name --remote-time "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/${alpha}/${beta}/${charlie}/${filename}/${filename}_genomic.fna.gz"
	else
		break
	fi
done < ${OUTDATADIR}/ANI/${1}_${REFSEQ_date}_mash_sorted.dists

successful_matches=$(ls -l "${OUTDATADIR}/ANI/localANIDB_REFSEQ" | wc -l)
if [[ ${successful_matches} -gt 2 ]]; then
"${shareScript}/append_taxonomy_to_ncbi_assembly_filenames.sh" "${OUTDATADIR}/ANI/localANIDB_REFSEQ"
gunzip ${OUTDATADIR}/ANI/localANIDB_REFSEQ/*.gz
else
echo "No matches were found against REFSEQ Bacterial Database sketch"
echo "0.00%ID-0.00%COV-NO_MATCHES_FOUND_AGAINST_BACTERIAL_REFSEQ(NONE)" > "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${1}_vs_${REFSEQ_date}).txt"
exit
fi

#Renames all files in the localANIDB_REFSEQ folder by changing extension from fna to fasta (which pyani needs)
for file in ${OUTDATADIR}/ANI/localANIDB_REFSEQ/*.fna;
do
	fasta_name=$(basename "${file}" .fna)".fasta"
	mv "${file}" "${OUTDATADIR}/ANI/localANIDB_REFSEQ/${fasta_name}"
done

# Checks for a previous copy of the aniM folder, removes it if found
if [ -d "${OUTDATADIR}/ANI/aniM_REFSEQ" ]; then
	echo "Removing old ANIm results in ${OUTDATADIR}/ANI/aniM_REFSEQ"
	rm -r "${OUTDATADIR}/ANI/aniM_REFSEQ"
fi

#Calls pyani on local db folder
python -V
#python "${shareScript}/pyani/average_nucleotide_identity.py" -i "${OUTDATADIR}/ANI/localANIDB_REFSEQ" -o "${OUTDATADIR}/ANI/aniM" --write_excel
average_nucleotide_identity.py -i "${OUTDATADIR}/ANI/localANIDB_REFSEQ" -o "${OUTDATADIR}/ANI/aniM_REFSEQ" --write_excel

#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line || [ -n "$line" ]; do
#	echo "!-${line}"
	if [[ ${line:0:6} = "sample" ]]; then
		sample_identity_line=${line}
#		echo "found it-"$sample_identity_line
		break
	fi
done < "${OUTDATADIR}/ANI/aniM_REFSEQ/ANIm_percentage_identity.tab"

#Extracts the query sample info line for percentage identity from the percent identity file
while IFS='' read -r line || [ -n "$line" ]; do
#	echo "!-${line}"
	if [[ ${line:0:6} = "sample" ]]; then
		sample_coverage_line=${line}
#		echo "found it-"$sample_identity_line
		break
	fi
done < "${OUTDATADIR}/ANI/aniM_REFSEQ/ANIm_alignment_coverage.tab"

#Extracts the top line from the %id file to get all the sample names used in analysis (they are tab separated along the top row)
if [[ -s "${OUTDATADIR}/ANI/aniM_REFSEQ/ANIm_percentage_identity.tab" ]]; then
	header_line=$(head -n 1 "${OUTDATADIR}/ANI/aniM_REFSEQ/ANIm_percentage_identity.tab")
else
	echo "No "${OUTDATADIR}/ANI/aniM_REFSEQ/ANIm_percentage_identity.tab" file, exiting"
	exit 1
fi

#Arrays to read sample names and the %ids for the query sample against those other samples
IFS="	" read -r -a samples <<< "${header_line}"
IFS="	" read -r -a percents <<< "${sample_identity_line}"
IFS="	" read -r -a coverages <<< "${sample_coverage_line}"

#How many samples were compared
n=${#samples[@]}

#Extracts all %id against the query sample (excluding itself) and writes them to file
for (( i=0; i<n; i++ ));
do
#	echo ${i}-${samples[i]}
	if [[ ${samples[i]:0:6} = "sample" ]];
	then
#		echo "Skipping ${i}"
		continue
	fi
	definition=$(head -1 "${OUTDATADIR}/ANI/localANIDB_REFSEQ/${samples[i]}.fasta")
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
best_file=${def_array[2]}
#Formats the %id to standard percentage (xx.xx%)
best_percent=$(awk -v per="${def_array[0]}" 'BEGIN{printf "%.2f", per * 100}')
best_coverage=$(awk -v per="${def_array[1]}" 'BEGIN{printf "%.2f", per * 100}')
#echo "${best_file}"

# #Extracts the accession number from the definition line
# accession=$(echo "${def_array[3]}" | cut -d' ' -f1  | cut -d'>' -f2)
# echo "Trying to find taxonomy of ${accession}"
# #Looks up the NCBI genus species from the accession number
# if [[ "${accession}" == "No_Accession_Number" ]]; then
# 	best_organism_guess="${def_array[4]} ${def_array[5]}"
# else
# 	attempts=0
# 	while [[ ${attempts} -lt 25 ]]; do
# 		best_organism_guess=$(python "${shareScript}/entrez_get_taxon_from_accession.py" -a "${accession}" -e "${me}")
# 		if [[ ! -z ${best_organism_guess} ]]; then
# 			break
# 		else
# 			attempts=$(( attempts + 1 ))
# 		fi
# 	done
# fi

# Pulling taxonomy from filename which was looked up. Can possibly be out of date. REFSEQ file will ALWAYS be current though
best_genus=$(echo "${best_file}" | cut -d'_' -f1)
best_species=$(echo "${best_file}" | cut -d'_' -f2)
best_organism_guess="${best_genus} ${best_species}"

# Uncomment this if you want to restrict ID to only genus species, without more resolute definition
#best_organism_guess_arr=($best_organism_guess})
#best_organism_guess="${best_organism_guess_arr[@]:0:2}"

#Creates a line at the top of the file to show the best match in an easily readable format that matches the style on the MMB_Seq log
echo -e "${best_percent}%ID-${best_coverage}%COV-${best_organism_guess}(${best_file}.fna)\\n$(cat "${OUTDATADIR}/ANI/best_hits_ordered.txt")" > "${OUTDATADIR}/ANI/best_ANI_hits_ordered(${1}_vs_REFSEQ_${REFSEQ_date}).txt"

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
