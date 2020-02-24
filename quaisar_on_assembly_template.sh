#!/bin/bash -l

#$ -o qoa_X.out
#$ -e qao_X.err
#$ -N qaoX
#$ -pe smp 12
#$ -cwd
#$ -q short.q

#
# Description: Alternate version of the main QuAISAR-H pipeline that (re)starts from after assembly step, project/isolate_name/Assembly must already have an assembly file (scaffolds.fasta) to work with
# 	This script assumes the sample is located in the default location ($processed) specified within the config file
#
# Usage: ./quaisar_on_assembly_template.sh isolate_name project_name
#
# Output location: default_config.sh_output_location
#
# Modules required: None
#
# v1.0.1 (10/22/2019)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./quaisar_on_assembly.sh  sample_name miseq_run_ID(or_project_name) config_file_to_use optional_alternate_directory"
	echo "Output by default is processed to processed/miseq_run_ID/sample_name"
	exit 0
elif [[ -z "{2}" ]]; then
	echo "No Project/Run_ID supplied to quaisar_template.sh, exiting"
	exit 33
fi

if [[ ! -f "${3}" ]]; then
	echo "no config file to load (${3}), exiting"
	exit 223
else
	echo "${2}/${1} is loading config file ${3}"
	. "${3}"
fi


#Time tracker to gauge time used by each step
totaltime=0
start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Set arguments to filename(sample name) project (miseq run id) and outdatadir(${processed}/project/filename)
filename="${1}"
project="${2}"
OUTDATADIR="${processed}/${project}"
if [[ ! -z "${4}" ]]; then
	OUTDATADIR="${4}/${project}"
fi

### Removing Short Contigs  ###
echo "----- Removing Short Contigs -----"
python3 "${shareScript}/removeShortContigs.py" -i "${OUTDATADIR}/${filename}/Assembly/scaffolds.fasta" -t 500 -s "normal_SPAdes"
mv "${OUTDATADIR}/${filename}/Assembly/scaffolds.fasta.TRIMMED.fasta" "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed.fasta"

# Checks to see that the trimming and renaming worked properly, returns if unsuccessful
if [ ! -s "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed.fasta" ]; then
	echo "Trimmed contigs file does not exist continuing to next sample">&2
	return 1
fi

### ReKraken on Assembly ###
echo "----- Running Kraken on Assembly -----"
# Get start time of kraken on assembly
start=$SECONDS
# Run kraken on assembly
"${shareScript}/run_kraken.sh" "${filename}" post assembled "${project}"
# Get end time of kraken on assembly and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeKrakAss=$((end - start))
echo "Kraken on Assembly - ${timeKrakAss} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeKrakAss))

# Get ID fom 16s
echo "----- Identifying via 16s blast -----"
start=$SECONDS
"${shareScript}/16s_blast.sh" "-n" "${filename}" "-p" "${project}"
end=$SECONDS
time16s=$((end - start))
echo "16S - ${time16s} seconds" >> "${time_summary}"
totaltime=$((totaltime + time16s))

# Get taxonomy from currently available files (Only ANI, has not been run...yet, will change after discussions)
"${shareScript}/determine_taxID.sh" "${filename}" "${project}"
# Capture the anticipated taxonomy of the sample using kraken on assembly output
echo "----- Extracting Taxonomy from Taxon Summary -----"
# Checks to see if the kraken on assembly completed successfully
if [ -s "${OUTDATADIR}/${filename}/${filename}.tax" ]; then
	# Read each line of the kraken summary file and pull out each level  taxonomic unit and store for use later in busco and ANI
	while IFS= read -r line  || [ -n "$line" ]; do
		# Grab first letter of line (indicating taxonomic level)
		first=${line::1}
		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "s" ]
		then
			species=$(echo "${line}" | awk -F '	' '{print $2}')
		elif [ "${first}" = "G" ]
		then
			genus=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "F" ]
		then
			family=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "O" ]
		then
			order=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "C" ]
		then
			class=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "P" ]
		then
			phylum=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "K" ]
		then
			kingdom=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "D" ]
		then
			domain=$(echo "${line}" | awk -F ' ' '{print $2}')
		fi
	done < "${OUTDATADIR}/${filename}/${filename}.tax"
# Print out taxonomy for confirmation/fun
echo "Taxonomy - ${domain} ${kingdom} ${phylum} ${class} ${order} ${family} ${genus} ${species}"
# If no kraken summary file was found
else
	echo "No Taxonomy output available to make best call from, skipped"
fi

### Check quality of Assembly ###
echo "----- Running quality checks on Assembly -----"
# Get start time of QC assembly check
start=$SECONDS
# Run qc assembly check
"${shareScript}/run_Assembly_Quality_Check.sh" "${filename}" "${project}"
# Get end time of qc quality check and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeQCcheck=$((end - start))
echo "QC Check of Assembly - ${timeQCcheck} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeQCcheck))

### Prokka on assembly ###
echo "----- Running Prokka on Assembly -----"
# Get start time for prokka
start=$SECONDS
# Run prokka
"${shareScript}/run_prokka.sh" "${filename}" "${project}"
# Get end time of prokka and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeProk=$((end - start))
echo "Identifying Genes (Prokka) - ${timeProk} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeProk))

# Rename contigs to something helpful (Had to wait until after prokka runs due to the strict naming requirements
mv "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed.fasta" "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed_original.fasta"
python3 "${shareScript}/fasta_headers.py" -i "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed_original.fasta" -o "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed.fasta"

# Rename contigs to something helpful (Had to wait until after prokka runs due to the strict naming requirements
mv "${OUTDATADIR}/${filename}/Assembly/${filename}_contigs_trimmed.fasta" "${OUTDATADIR}/${filename}/Assembly/${filename}_contigs_trimmed_original.fasta"
python3 "${shareScript}/fasta_headers.py" -i "${OUTDATADIR}/${filename}/Assembly/${filename}_contigs_trimmed_original.fasta" -o "${OUTDATADIR}/${filename}/Assembly/${filename}_contigs_trimmed.fasta"

### Average Nucleotide Identity ###
echo "----- Running ANI for Species confirmation -----"
# ANI uses assembly and sample would have exited already if assembly did not complete, so no need to check
# Get start time of ANI
start=$SECONDS
# run ANI
# Temp fix for strange genera until we do vs ALL all the time.
if [[ "${genus}" = "Peptoclostridium" ]] || [[ "${genus}" = "Clostridioides" ]]; then
	genus="Clostridium"
elif [[ "${genus}" = "Shigella" ]]; then
	genus="Escherichia"
fi
"${shareScript}/run_ANI.sh" "${filename}" "${genus}" "${species}" "${project}"
# Get end time of ANI and calculate run time and append to time summary (and sum to total time used
end=$SECONDS
timeANI=$((end - start))
echo "autoANI - ${timeANI} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeANI))

# Get taxonomy from currently available files (Only ANI, has not been run...yet, will change after discussions)
"${shareScript}/determine_taxID.sh" "${filename}" "${project}"
"${OUTDATADIR}/${filename}/${filename}.tax"

### BUSCO on prokka output ###
echo "----- Running BUSCO on Assembly -----"
# Check to see if prokka finished successfully
if [ -s "${OUTDATADIR}/${filename}/prokka/${filename}_PROKKA.gbf" ] || [ -s "${OUTDATADIR}/${filename}/prokka/${filename}_PROKKA.gff" ]; then
	# Get start time of busco
	start=$SECONDS
	# Set default busco database as bacteria in event that we dont have a database match for sample lineage
	buscoDB="bacteria_odb10"
	# Iterate through taxon levels (species to domain) and test if a match occurs to entry in database. If so, compare against it
	busco_found=0
	for tax in $species $genus $family $order $class $phylum $kingdom $domain
	do
		if [ -d "${local_DBs}/BUSCO/${tax,}_odb10" ]
		then
			buscoDB="${tax,}_odb10"
			busco_found=1
			break
		fi
	done
	# Report an unknown sample to the maintenance file to look into
	if [[ "${busco_found}" -eq 0 ]]; then
		global_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
		echo "BUSCO: ${domain} ${kingdom} ${phylum} ${class} ${order} ${family} ${genus} ${species} - Found as ${project}/${filename} on ${global_time}" >> "${shareScript}/maintenance_To_Do.txt"
	fi
	# Show which database entry will be used for comparison
	echo "buscoDB:${buscoDB}"
	# Run busco
	"${shareScript}/do_busco.sh" "${filename}" "${buscoDB}" "${project}"
	# Get end time of busco and calculate run time and append to time summary (and sum to total time used
	end=$SECONDS
	timeBUSCO=$((end - start))
	echo "BUSCO - ${timeBUSCO} seconds" >> "${time_summary}"
	totaltime=$((totaltime + timeBUSCO))
# Prokka did not complete successfully and busco cant run (since it relies on prokka output)
else
	echo "Prokka output not found, not able to process BUSCO"
fi

### c-SSTAR for finding AR Genes ###
echo "----- Running c-SSTAR for AR Gene identification -----"
# c-SSTAR uses assembly and sample would have exited already if assembly did not complete, so no need to check
# Get start time of ccstar
start=$SECONDS

# Run csstar in default mode from config.sh
"${shareScript}/run_c-sstar.sh" "${filename}" "${csstar_gapping}" "${csstar_identity}" "${project}"
"${shareScript}/run_c-sstar_altDB.sh" "${filename}" "${csstar_gapping}" "${csstar_identity}" "${project}" "${local_DBs}/star/ResGANNOT_20180608_srst2.fasta"

### GAMA - finding AR Genes ###
echo "----- Running GAMA for AR Gene identification -----"
${shareScript}/run_GAMA.sh "${filename}" "${project}" -c

# Get end time of csstar and calculate run time and append to time summary (and sum to total time used
end=$SECONDS
timestar=$((end - start))
echo "c-SSTAR - ${timestar} seconds" >> "${time_summary}"
totaltime=$((totaltime + timestar))

# Get MLST profile
echo "----- Running MLST -----"
start=$SECONDS
"${shareScript}/run_MLST.sh" "${filename}" "${project}"
python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${filename}/MLST/${filename}.mlst" -t standard
if [[ "${genus}_${species}" = "Acinetobacter_baumannii" ]]; then
	"${shareScript}/run_MLST.sh" "${filename}" "${project}" "-f" "abaumannii"
	python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${filename}/MLST/${filename}_abaumannii.mlst" -t standard
	mv "${processed}/${project}/${filename}/MLST/${filename}_abaumannii.mlst" "${processed}/${project}/${filename}/MLST/${filename}_Oxford.mlst"
	mv "${processed}/${project}/${filename}/MLST/${filename}.mlst" "${processed}/${project}/${filename}/MLST/${filename}_Pasteur.mlst"
	#Check for "-", unidentified type
	type1=$(tail -n1 ${processed}/${project}/${filename}/MLST/${filename}_abaumannii.mlst | cut -d' ' -f3)
	type2=$(head -n1 ${processed}/${project}/${filename}/MLST/${filename}.mlst | cut -d' ' -f3)
	if [[ "${type1}" = "-" ]]; then
		"${shareScript}/run_srst2_mlst.sh" "${filename}" "${project}" "Acinetobacter" "baumannii#1"
		python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${filename}/MLST/${filename}_srst2_acinetobacter_baumannii-baumannii#1.mlst" -t srst2
	fi
	if [[ "${type2}" = "-" ]]; then
		"${shareScript}/run_srst2_mlst.sh" "${filename}" "${project}" "Acinetobacter" "baumannii#2"
		python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${filename}/MLST/${filename}_srst2_acinetobacter_baumannii-baumannii#2.mlst" -t srst2
	fi
elif [[ "${genus}_${species}" = "Escherichia_coli" ]]; then
	# Verify that ecoli_2 is default and change accordingly
	"${shareScript}/run_MLST.sh" "${filename}" "${project}" "-f" "ecoli_2"
	python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${filename}/MLST/${filename}_ecoli_2.mlst" -t standard
	mv "${processed}/${project}/${filename}/MLST/${filename}_ecoli_2.mlst" "${processed}/${project}/${filename}/MLST/${filename}_Pasteur.mlst"
	mv "${processed}/${project}/${filename}/MLST/${filename}.mlst" "${processed}/${project}/${filename}/MLST/${filename}_Achtman.mlst"
	type2=$(tail -n1 ${processed}/${project}/${filename}/MLST/${filename}_ecoli_2.mlst | cut -d' ' -f3)
	type1=$(head -n1 ${processed}/${project}/${filename}/MLST/${filename}.mlst | cut -d' ' -f3)
	if [[ "${type1}" = "-" ]]; then
		"${shareScript}/run_srst2_mlst.sh" "${filename}" "${project}" "Escherichia" "coli#1"
		python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${filename}/MLST/${filename}_srst2_escherichia_coli-coli#1.mlst" -t srst2
	fi
	if [[ "${type2}" = "-" ]]; then
		"${shareScript}/run_srst2_mlst.sh" "${filename}" "${project}" "Escherichia" "coli#2"
		python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${filename}/MLST/${filename}_srst2_escherichia_coli-coli#2.mlst" -t srst2
	fi
else
	mv "${processed}/${project}/${filename}/MLST/${filename}.mlst" "${processed}/${project}/${filename}/MLST/${filename}_Pasteur.mlst"
fi
end=$SECONDS
timeMLST=$((end - start))
echo "MLST - ${timeMLST} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeMLST))

# Try to find any plasmids
echo "----- Identifying plasmids using plasmidFinder -----"
start=$SECONDS
"${shareScript}/run_plasmidFinder.sh" "${filename}" "${project}" "plasmidFinder"
#"${shareScript}/run_plasmidFinder.sh" "${filename}" "${project}" "plasmid_on_plasFlow"
end=$SECONDS
timeplasfin=$((end - start))
echo "plasmidFinder - ${timeplasfin} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeplasfin))

# Run plasFlow if isolate is from the Enterobacteriaceae family  ##### When should we check if this will be expanded?
if [[ "${family}" == "Enterobacteriaceae" ]]; then
	start=$SECONDS
	${shareScript}/run_plasFlow.sh "${filename}" "${project}"
	${shareScript}/run_c-sstar_plasFlow.sh "${filename}" g o "${project}" -p
	${shareScript}/run_plasmidFinder.sh "${filename}" "${project}" plasmidFinder_on_plasFlow
	${shareScript}/run_GAMA.sh "${filename}" "${project}" -p
	end=$SECONDS
	timeplasflow=$((end - start))
	echo "plasmidFlow - ${timeplasflow} seconds" >> "${time_summary_redo}"
	totaltime=$((totaltime + timeplasflow))
fi

"${shareScript}/validate_piperun.sh" "${filename}" "${project}" > "${processed}/${project}/${filename}/${filename}_pipeline_stats.txt"

status=$(tail -n1 "${processed}/${project}/${filename}/${filename}_pipeline_stats.txt" | cut -d' ' -f5)
if [[ "${status}" != "FAILED" ]]; then
	"${shareScript}/sample_cleaner.sh" "${filename}" "${project}"
fi

# Extra dump cleanse in case anything else failed
	if [ -n "$(find "${shareScript}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
		echo "Found core dump files at end of processing ${filename} and attempting to delete"
		find "${shareScript}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
	fi

global_end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Append total time to bottom of time summary
echo "Total time: ${totaltime} seconds" >> "${time_summary}"
echo "Completed at ${global_end_time}"

# Designate end of this sample #
echo "

				End of sample ${filename}
				completed at ${global_end_time}

"
