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
# Usage: ./quaisar_on_assembly_template.sh -n isolate_name -p project_name [-c path_to_config_file]
#
# Output location: default_config.sh_output_location
#
# Modules required: None
#
# v1.0.2 (08/24/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

#  Function to print out help blurb
show_help () {
	echo "Usage 1: ./quaisar_failed_assembly.sh -n isolate_name -p project_name [-c path_to_config_file]"
}

# Parse command line options
options_found=0
while getopts ":h?c:p:n:" option; do
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
	echo "No Project/Run_ID supplied to quaisar_failed_assembly.sh, exiting"
	exit 33
elif [[ -z "${sample_name}" ]]; then
	echo "No sample name supplied to quaisar_failed_assembly.sh, exiting"
	exit 34
fi

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

# Set arguments to sample_name(sample name) project (miseq run id) and outdatadir(${processed}/project/sample_name)
OUTDATADIR="${processed}/${project}"

### Removing Short Contigs  ###
echo "----- Removing Short Contigs -----"
python3 "${shareScript}/removeShortContigs.py" -i "${OUTDATADIR}/${sample_name}/Assembly/scaffolds.fasta" -t 500 -s "normal_SPAdes" -c "${config}"
mv "${OUTDATADIR}/${sample_name}/Assembly/scaffolds.fasta.TRIMMED.fasta" "${OUTDATADIR}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta"

# Checks to see that the trimming and renaming worked properly, returns if unsuccessful
if [ ! -s "${OUTDATADIR}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" ]; then
	echo "Trimmed contigs file does not exist continuing to next sample">&2
	return 1
fi

### ReKraken on Assembly ###
echo "----- Running Kraken on Assembly -----"
# Get start time of kraken on assembly
start=$SECONDS
# Run kraken on assembly
"${shareScript}/run_kraken.sh" -n "${sample_name}" -r post -p "${project}" -c "${config}"
# Get end time of kraken on assembly and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeKrakAss=$((end - start))
echo "Kraken on Assembly - ${timeKrakAss} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeKrakAss))

# Get ID fom 16s
echo "----- Identifying via 16s blast -----"
start=$SECONDS
"${shareScript}/16s_blast.sh" -n "${sample_name}" -p "${project}" -c "${config}"
end=$SECONDS
time16s=$((end - start))
echo "16S - ${time16s} seconds" >> "${time_summary}"
totaltime=$((totaltime + time16s))

# Get taxonomy from currently available files (Only ANI, has not been run...yet, will change after discussions)
"${shareScript}/determine_taxID.sh" -n "${sample_name}" -p "${project}" -c "${config}"
# Capture the anticipated taxonomy of the sample using kraken on assembly output
echo "----- Extracting Taxonomy from Taxon Summary -----"
# Checks to see if the kraken on assembly completed successfully
if [ -s "${OUTDATADIR}/${sample_name}/${sample_name}.tax" ]; then
	# Read each line of the kraken summary file and pull out each level  taxonomic unit and store for use later in busco and ANI
	while IFS= read -r line  || [ -n "$line" ]; do
		# Grab first letter of line (indicating taxonomic level)
		first=${line::1}
		# Assign taxonomic level value from 4th value in line (1st-classification level,2nd-% by kraken, 3rd-true % of total reads, 4th-identifier)
		if [ "${first}" = "s" ]
		then
			species=$(echo "${line}" | awk -F '	' '{print $2}')
		elif [ "${first}" = "G" ]; then
			genus=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "F" ]; then
			family=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "O" ]; then
			order=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "C" ]; then
			class=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "P" ]; then
			phylum=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "K" ]; then
			kingdom=$(echo "${line}" | awk -F ' ' '{print $2}')
		elif [ "${first}" = "D" ]; then
			domain=$(echo "${line}" | awk -F ' ' '{print $2}')
		fi
	done < "${OUTDATADIR}/${sample_name}/${sample_name}.tax"
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
"${shareScript}/run_Assembly_Quality_Check.sh" -n "${sample_name}" -p "${project}" -c "${config}"
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
"${shareScript}/run_prokka.sh" -n "${sample_name}" -p "${project}" -c "${config}"
# Get end time of prokka and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeProk=$((end - start))
echo "Identifying Genes (Prokka) - ${timeProk} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeProk))

# Rename contigs to something helpful (Had to wait until after prokka runs due to the strict naming requirements
mv "${OUTDATADIR}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta" "${OUTDATADIR}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed_original.fasta"
python3 "${shareScript}/fasta_headers.py" -i "${OUTDATADIR}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed_original.fasta" -o "${OUTDATADIR}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta"

# Rename contigs to something helpful (Had to wait until after prokka runs due to the strict naming requirements
mv "${OUTDATADIR}/${sample_name}/Assembly/${sample_name}_contigs_trimmed.fasta" "${OUTDATADIR}/${sample_name}/Assembly/${sample_name}_contigs_trimmed_original.fasta"
python3 "${shareScript}/fasta_headers.py" -i "${OUTDATADIR}/${sample_name}/Assembly/${sample_name}_contigs_trimmed_original.fasta" -o "${OUTDATADIR}/${sample_name}/Assembly/${sample_name}_contigs_trimmed.fasta"

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

#"${shareScript}/run_ANI.sh" -n "${sample_name}" -g "${genus}" -s "${species}" -p "${project}" -c "${config}"
"${shareScript}/run_ANI.sh" -n "${sample_name}" -p "${project}" -c "${config}"

# Get end time of ANI and calculate run time and append to time summary (and sum to total time used
end=$SECONDS
timeANI=$((end - start))
echo "autoANI - ${timeANI} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeANI))

# Get taxonomy from currently available files (Only ANI, has not been run...yet, will change after discussions)
"${shareScript}/determine_taxID.sh" -n "${sample_name}" -p "${project}" -c "${config}"
"${OUTDATADIR}/${sample_name}/${sample_name}.tax"

### BUSCO on prokka output ###
echo "----- Running BUSCO on Assembly -----"
# Check to see if prokka finished successfully
if [ -s "${OUTDATADIR}/${sample_name}/prokka/${sample_name}_PROKKA.gbf" ] || [ -s "${OUTDATADIR}/${sample_name}/prokka/${sample_name}_PROKKA.gff" ]; then
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
		echo "BUSCO: ${domain} ${kingdom} ${phylum} ${class} ${order} ${family} ${genus} ${species} - Found as ${project}/${sample_name} on ${global_time}" >> "${shareScript}/maintenance_To_Do.txt"
	fi
	# Show which database entry will be used for comparison
	echo "buscoDB:${buscoDB}"
	# Run busco
	"${shareScript}/do_busco.sh" -n "${sample_name}" -d "${buscoDB}" -p "${project}" -c "${config}"
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
"${shareScript}/run_c-sstar.sh" -n "${sample_name}" -g "${csstar_gapping}" -s "${csstar_identity}" -p "${project}"
"${shareScript}/run_c-sstar.sh" -n "${sample_name}" -g "${csstar_gapping}" -s "${csstar_identity}" -p "${project}" -d "${local_DBs}/star/ResGANNOT_20180608_srst2.fasta"

### GAMA - finding AR Genes ###
echo "----- Running GAMA for AR Gene identification -----"
${shareScript}/run_GAMA.sh -n "${sample_name}" -p "${project}"  -c "${config}"

# Get end time of csstar and calculate run time and append to time summary (and sum to total time used
end=$SECONDS
timestar=$((end - start))
echo "c-SSTAR - ${timestar} seconds" >> "${time_summary}"
totaltime=$((totaltime + timestar))

# Get MLST profile
echo "----- Running MLST -----"
start=$SECONDS
"${shareScript}/run_MLST.sh" -n "${sample_name}" -p "${project}" -c "${config}"
python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" -t standard
if [[ "${genus}_${species}" = "Acinetobacter_baumannii" ]]; then
	"${shareScript}/run_MLST.sh" -n "${sample_name}" -p "${project}" "-f" "abaumannii"  -c "${config}"
	python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst" -t standard
	mv "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst" "${processed}/${project}/${sample_name}/MLST/${sample_name}_Oxford.mlst"
	mv "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" "${processed}/${project}/${sample_name}/MLST/${sample_name}_Pasteur.mlst"
	#Check for "-", unidentified type
	type1=$(tail -n1 ${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst | cut -d' ' -f3)
	type2=$(head -n1 ${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst | cut -d' ' -f3)
	if [[ "${type1}" = "-" ]]; then
		"${shareScript}/run_srst2_mlst.sh" -n "${sample_name}" -p "${project}" -g "Acinetobacter" -s "baumannii#1"  -c "${config}"
		python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${sample_name}/MLST/${sample_name}_srst2_acinetobacter_baumannii-baumannii#1.mlst" -t srst2
	fi
	if [[ "${type2}" = "-" ]]; then
		"${shareScript}/run_srst2_mlst.sh" -n "${sample_name}" -p "${project}" -g "Acinetobacter" -s "baumannii#2"  -c "${config}"
		python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${sample_name}/MLST/${sample_name}_srst2_acinetobacter_baumannii-baumannii#2.mlst" -t srst2
	fi
elif [[ "${genus}_${species}" = "Escherichia_coli" ]]; then
	# Verify that ecoli_2 is default and change accordingly
	"${shareScript}/run_MLST.sh" -n "${sample_name}" -p "${project}" "-f" "ecoli_2"  -c "${config}"
	python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst" -t standard
	mv "${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst" "${processed}/${project}/${sample_name}/MLST/${sample_name}_Pasteur.mlst"
	mv "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" "${processed}/${project}/${sample_name}/MLST/${sample_name}_Achtman.mlst"
	type2=$(tail -n1 ${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst | cut -d' ' -f3)
	type1=$(head -n1 ${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst | cut -d' ' -f3)
	if [[ "${type1}" = "-" ]]; then
		"${shareScript}/run_srst2_mlst.sh" -n "${sample_name}" -p "${project}" -g "Escherichia" -s "coli#1"  -c "${config}"
		python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${sample_name}/MLST/${sample_name}_srst2_escherichia_coli-coli#1.mlst" -t srst2
	fi
	if [[ "${type2}" = "-" ]]; then
		"${shareScript}/run_srst2_mlst.sh" -n "${sample_name}" -p "${project}" -g "Escherichia" -s "coli#2" -c "${config}"
		python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${sample_name}/MLST/${sample_name}_srst2_escherichia_coli-coli#2.mlst" -t srst2
	fi
else
	mv "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" "${processed}/${project}/${sample_name}/MLST/${sample_name}_Pasteur.mlst"
fi
end=$SECONDS
timeMLST=$((end - start))
echo "MLST - ${timeMLST} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeMLST))

# Try to find any plasmids
echo "----- Identifying plasmids using plasmidFinder -----"
start=$SECONDS
"${shareScript}/run_plasmidFinder.sh" -n "${sample_name}" -p "${project}" -o "plasmidFinder" -c "${config}"
#"${shareScript}/run_plasmidFinder.sh" "${sample_name}" "${project}" "plasmid_on_plasFlow"
end=$SECONDS
timeplasfin=$((end - start))
echo "plasmidFinder - ${timeplasfin} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeplasfin))

# Run plasFlow if isolate is from the Enterobacteriaceae family  ##### When should we check if this will be expanded?
if [[ "${family}" == "Enterobacteriaceae" ]]; then
	start=$SECONDS
	${shareScript}/run_plasFlow.sh -n "${sample_name}" -p "${project}" -c "${config}"
	${shareScript}/run_Assembly_Quality_Check.sh -n "${sample_name}" -p "${project}" -l -c "${config}"
	${shareScript}/run_c-sstar.sh "${sample_name}" -g g -s o -p "${project}" -l -c "${config}"
	${shareScript}/run_plasmidFinder.sh -n "${sample_name}" -p "${project}" -o plasmidFinder_on_plasFlow -c "${config}"
	${shareScript}/run_GAMA.sh -n "${sample_name}" -p "${project}" -l -c "${config}"
	end=$SECONDS
	timeplasflow=$((end - start))
	echo "plasmidFlow - ${timeplasflow} seconds" >> "${time_summary_redo}"
	totaltime=$((totaltime + timeplasflow))
fi

"${shareScript}/validate_piperun.sh" -n "${sample_name}" -p "${project}" -c "${config}" > "${processed}/${project}/${sample_name}/${sample_name}_pipeline_stats.txt"

status=$(tail -n1 "${processed}/${project}/${sample_name}/${sample_name}_pipeline_stats.txt" | cut -d' ' -f5)
if [[ "${status}" != "FAILED" ]]; then
	"${shareScript}/sample_cleaner.sh" -n "${sample_name}" -p "${project}" -c "${config}"
fi

# Extra dump cleanse in case anything else failed
	if [ -n "$(find "${shareScript}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
		echo "Found core dump files at end of processing ${sample_name} and attempting to delete"
		find "${shareScript}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
	fi

global_end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Append total time to bottom of time summary
echo "Total time: ${totaltime} seconds" >> "${time_summary}"
echo "Completed at ${global_end_time}"

# Designate end of this sample #
echo "

				End of sample ${sample_name}
				completed at ${global_end_time}

"
