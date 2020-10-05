#!/bin/bash -l

#$ -o quaisar_X.out
#$ -e quaisar_X.err
#$ -N quasX
#$ -pe smp 12
#$ -cwd
#$ -q all.q

#
# Description: The full QuAISAR-H pipeline start to end serially, project/isolate_name must already have a populated FASTQs folder to work with
#
# Usage: ./quaisar_template.sh -n isolate_name -p project_name [-c path_to_config_file_to_use]
#
# Output location: default_config.sh_output_location
#
# Modules required: Python3/3.5.2
#
# v1.0.5 (08/24/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

version_type="Quaisar Parallel"
version_num="qp1.0.5"


#  Function to print out help blurb
show_help () {
	echo "Usage 1: ./quaisar_template.sh -n isolate_name -p project_name [-c path_to_config_file_to_use]"
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

ml Python3/3.5.2

#Time tracker to gauge time used by each step
totaltime=0
start_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")

# Set arguments to sample_name(sample name) project (miseq run id) and outdatadir(${processed}/project/sample_name)
OUTDATADIR="${processed}/${project}"

# Remove old run stats as the presence of the file indicates run completion
if [[ -f "${OUTDATADIR}/${sample_name}/${sample_name}_pipeline_stats.txt" ]]; then
	rm "${OUTDATADIR}/${sample_name}/${sample_name}_pipeline_stats.txt"
fi

# Create an empty time_summary file that tracks clock time of tools used
touch "${OUTDATADIR}/${sample_name}/${sample_name}_time_summary.txt"
time_summary=${OUTDATADIR}/${sample_name}/${sample_name}_time_summary.txt

echo "QuAISAR v${version}" > "${time_summary}"
echo "Time summary for ${project}/${sample_name}: Started ${global_time}" >> "${time_summary}"
echo "${project}/${sample_name} started at ${global_time}"

echo "Starting processing of ${project}/${sample_name}"
#Checks if FASTQ folder exists for current sample
if [[ -d "${OUTDATADIR}/${sample_name}/FASTQs" ]]; then
	# Checks if FASTQ folder contains any files then continue
	if [[ "$(ls -A "${OUTDATADIR}/${sample_name}/FASTQs")" ]]; then
		# Checks to see if those files in the folder are unzipped fastqs
		count_unzip=`ls -1 ${OUTDATADIR}/${sample_name}/FASTQs/*.fastq 2>/dev/null | wc -l`
		count_zip=`ls -1 ${OUTDATADIR}/${sample_name}/FASTQs/*.fastq.gz 2>/dev/null | wc -l`
		if [[ ${count_unzip} != 0 ]]; then
		#if [[ -f "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}"*".fastq" ]]; then
			echo "----- FASTQ(s) exist, continuing analysis -----"
			if [[ -f "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq" ]] && [[ ! -f "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz" ]]; then
				gzip < "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq" > "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz"
			fi
			if [[ -f "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq" ]] && [[ ! -f "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz" ]]; then
				gzip < "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq" > "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz"
			fi
		# Checks if they are zipped fastqs (checks for R1 first)
		elif [[ -f "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz" ]]; then
			#echo "R1 zipped exists - unzipping"
			gunzip -c "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq.gz" > "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq"
			# Checks for paired R2 file
			if [[ -f "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz" ]]; then
				#echo "R2 zipped exists - unzipping"
				gunzip -c "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz" > "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq"
			else
				echo "No matching R2 to unzip :("
			fi
		# Checks to see if there is an abandoned R2 zipped fastq
		elif [[ -f "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz" ]]; then
			#echo "R2 zipped  exists - unzipping"
			gunzip -c "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq.gz" > "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq"
			echo "No matching R1 to unzip :("
		fi
	# If the folder is empty then return from function
	else
		echo "FASTQs folder empty - No fastqs available for ${sample_name} (and download was not requested). Either unzip fastqs to ${OUTDATADIR}/FASTQs or run the -d flag to trigger unzipping of gzs"
		return 1
	fi
# If the fastq folder does not exist then return out of function
else
	echo "FASTQs not downloaded and FASTQs folder does not exist for ${sample_name}. No fastqs available (and download was not requested). Unzip fastqs to ${OUTDATADIR}/FASTQs"
	return 1
fi

# Get start time for qc check
start=$SECONDS
### Count the number of Q20, Q30, bases and reads within a pair of FASTQ files
echo "----- Counting read quality -----"
# Checks for and creates the specified output folder for the QC counts
if [ ! -d "${OUTDATADIR}/${sample_name}/preQCcounts" ]; then
	echo "Creating ${OUTDATADIR}/${sample_name}/preQCcounts"
	mkdir -p "${OUTDATADIR}/${sample_name}/preQCcounts"
fi
# Run qc count check on raw reads
echo -e "Q20_Total_[bp]	Q30_Total_[bp]	Q20_R1_[bp]	Q20_R2_[bp]	Q20_R1_[%]	Q20_R2_[%]	Q30_R1_[bp]	Q30_R2_[bp]	Q30_R1_[%]	Q30_R2_[%]	Total_Sequenced_[bp]	Total_Sequenced_[reads]" > "${OUTDATADIR}/${sample_name}/preQCcounts/${sample_name}_counts.txt"
python3 "${shareScript}/Fastq_Quality_Printer.py" -1 "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq" -2 "${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq" >> "${OUTDATADIR}/${sample_name}/preQCcounts/${sample_name}_counts.txt"

	# Get end time of qc count and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeQCcount=$((end - start))
echo "QC count - ${timeQCcount} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeQCcount))

###  Trimming and Quality Control  ###
echo "----- Running BBDUK on reads -----"
# Gets start time for bbduk
start=$SECONDS
# Creates folder for BBDUK output
if [ ! -d "${OUTDATADIR}/${sample_name}/removedAdapters" ]; then
	echo "Creating ${OUTDATADIR}/${sample_name}/removedAdapters"
	mkdir -p "${OUTDATADIR}/${sample_name}/removedAdapters"
# It complains if a folder already exists, so the current one is removed (shouldnt happen anymore as each analysis will move old runs to new folder)
else
	echo "Removing old ${OUTDATADIR}/${sample_name}/removedAdapters"
	rm -r "${OUTDATADIR}/${sample_name}/removedAdapters"
	echo "Recreating ${OUTDATADIR}/${sample_name}/removedAdapters"
	mkdir -p "${OUTDATADIR}/${sample_name}/removedAdapters"
fi

# Run bbduk
#bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R1_001.fastq" in2="${OUTDATADIR}/${sample_name}/FASTQs/${sample_name}_R2_001.fastq" out="${OUTDATADIR}/${sample_name}/removedAdapters/${sample_name}-noPhiX-R1.fsq" out2="${OUTDATADIR}/${sample_name}/removedAdapters/${sample_name}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
${shareScript}/run_BBDUK.sh -n "${sample_name}" -p "${project}" -c "${config}"
# Get end time of bbduk and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeAdapt=$((end - start))
echo "Removing Adapters - ${timeAdapt} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeAdapt))


### Quality and Adapter Trimming using trimmomatic ###
echo "----- Running Trimmomatic on reads -----"
# Get start time of trimmomatic
start=$SECONDS
# Creates folder for trimmomatic output if it does not exist
if [ ! -d "${OUTDATADIR}/${sample_name}/trimmed" ]; then
	mkdir -p "${OUTDATADIR}/${sample_name}/trimmed"
fi

ml trimmomatic/0.35
# Run trimmomatic
trimmomatic "${trim_endtype}" -"${trim_phred}" -threads "${procs}" "${OUTDATADIR}/${sample_name}/removedAdapters/${sample_name}-noPhiX-R1.fsq" "${OUTDATADIR}/${sample_name}/removedAdapters/${sample_name}-noPhiX-R2.fsq" "${OUTDATADIR}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq" "${OUTDATADIR}/${sample_name}/trimmed/${sample_name}_R1_001.unpaired.fq" "${OUTDATADIR}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq" "${OUTDATADIR}/${sample_name}/trimmed/${sample_name}_R2_001.unpaired.fq" ILLUMINACLIP:"${trim_adapter_location}:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome}" SLIDINGWINDOW:"${trim_window_size}:${trim_window_qual}" LEADING:"${trim_leading}" TRAILING:"${trim_trailing}" MINLEN:"${trim_min_length}"
# Get end time of trimmomatic and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeTrim=$((end - start))
echo "Trimming - ${timeTrim} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeTrim))
ml -trimmomatic/0.35

# Check differences after QC and trimming (also for gottcha proper read count for assessing unclassified reads)
# Get start time for qc check on trimmed reads
start=$SECONDS
### Count the number of Q20, Q30, bases and reads within the trimmed pair of FASTQ files
echo "----- Counting read quality of trimmed files-----"
# Checks for and creates the specified output folder for the QC counts
if [ ! -d "${OUTDATADIR}/${sample_name}/preQCcounts" ]; then
	echo "Creating ${OUTDATADIR}/${sample_name}/preQCcounts"
	mkdir -p "${OUTDATADIR}/${sample_name}/preQCcounts"
fi
# Run qc count check on filtered reads
echo -e "Q20_Total_[bp]	Q30_Total_[bp]	Q20_R1_[bp]	Q20_R2_[bp]	Q20_R1_[%]	Q20_R2_[%]	Q30_R1_[bp]	Q30_R2_[bp]	Q30_R1_[%]	Q30_R2_[%]	Total_Sequenced_[bp]	Total_Sequenced_[reads]" > "${OUTDATADIR}/${sample_name}/preQCcounts/${sample_name}_trimmed_counts.txt"
python3 "${shareScript}/Fastq_Quality_Printer.py" -1 "${OUTDATADIR}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq" -2 "${OUTDATADIR}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq" >> "${OUTDATADIR}/${sample_name}/preQCcounts/${sample_name}_trimmed_counts.txt"

# Merge both unpaired fq files into one for GOTTCHA
cat "${OUTDATADIR}/${sample_name}/trimmed/${sample_name}_R1_001.unpaired.fq" "${OUTDATADIR}/${sample_name}/trimmed/${sample_name}_R2_001.unpaired.fq" > "${OUTDATADIR}/${sample_name}/trimmed/${sample_name}.single.fq"


# Get end time of qc count and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeQCcount_trimmed=$((end - start))
echo "QC count trimmed - ${timeQCcount_trimmed} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeQCcount))



######  Run Kraken on cleaned reads  ######
echo "----- Running Kraken on cleaned reads -----"
# Get start time of kraken on reads
start=$SECONDS
# Run kraken
"${shareScript}/run_kraken.sh" -n "${sample_name}" -r pre -p "${project}" -c "${config}"
# Get end time of kraken on reads and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeKrak=$((end - start))
echo "Kraken - ${timeKrak} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeKrak))

##### Run gottcha(v1) on cleaned reads #####
echo "----- Running gottcha on cleaned reads -----"
# Get start time of gottcha
start=$SECONDS
# run gootcha
"${shareScript}/run_gottcha.sh" -n "${sample_name}" -p "${project}" -c "${config}"
# Get end time of qc count and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeGott=$((end - start))
echo "Gottcha - ${timeGott} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeGott))

# Check reads using SRST2
echo "----- Running SRST2 -----"
start=$SECONDS
"${shareScript}/run_srst2AR.sh" -n "${sample_name}" -p "${project}" -c "${config}"
"${shareScript}/run_srst2AR.sh" -n "${sample_name}" -p "${project}" -d "${local_DBs}/star/ResGANNOT_20180608_srst2.fasta" -c "${config}"
end=$SECONDS
timesrst2=$((end - start))
echo "SRST2 - ${timesrst2} seconds" >> "${time_summary}"
totaltime=$((totaltime + timesrst2))

######  Assembling Using SPAdes  ######
echo "----- Assembling Using SPAdes -----"
# Get start time of SPAdes
start=$SECONDS
# script tries 3 times for a completed assembly
for i in 1 2 3 4 5 6 7 8 9
do
	# If assembly exists already and this is the first attempt (then the previous run will be used) [should not happen anymore as old runs are now renamed]
	if [ -s "${OUTDATADIR}/${sample_name}/Assembly/scaffolds.fasta" ]; then
		echo "Previous assembly already exists, using it (delete/rename the assembly folder at ${OUTDATADIR}/ if you'd like to try to reassemble"
	# Run normal mode if no assembly file was found
	else
		# Cant figure out why continue does not work yet
		#if [[ "${1}" -gt 1 ]]; then
		#	"${shareScript}/run_SPAdes.sh" "${sample_name}" "continue" "${project}"
		#else
		#	"${shareScript}/run_SPAdes.sh" "${sample_name}" normal "${project}"
		#fi
		"${shareScript}/run_SPAdes.sh" -n "${sample_name}" -t normal -p "${project}" -c "${config}"
	fi
	# Removes any core dump files (Occured often during testing and tweaking of memory parameter
	if [ -n "$(find "${shareScript}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
		echo "Found core dump files in assembly (assumed to be from SPAdes, but could be anything before that as well) and attempting to delete"
		find "${shareScript}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
	fi
done
# Returns if all 3 assembly attempts fail
if [[ -f "${OUTDATADIR}/${sample_name}/Assembly/scaffolds.fasta" ]] && [[ -s "${OUTDATADIR}/${sample_name}/Assembly/scaffolds.fasta" ]]; then
	echo "Assembly completed and created a non-empty scaffolds file"
else
	echo "Assembly FAILED 3 times, continuing to next sample..." >&2
	return 1
fi

# Get end time of SPAdes and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeSPAdes=$((end - start))
echo "SPAdes - ${timeSPAdes} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeSPAdes))

### Removing Short Contigs  ###
echo "----- Removing Short Contigs -----"
python3 "${shareScript}/removeShortContigs.py" -i "${OUTDATADIR}/${sample_name}/Assembly/scaffolds.fasta" -t 500 -s "normal_SPAdes"
mv "${OUTDATADIR}/${sample_name}/Assembly/scaffolds.fasta.TRIMMED.fasta" "${OUTDATADIR}/${sample_name}/Assembly/${sample_name}_scaffolds_trimmed.fasta"

### Removing Short Contigs  ###
echo "----- Removing Short Contigs -----"
python3 "${shareScript}/removeShortContigs.py" -i "${OUTDATADIR}/${sample_name}/Assembly/contigs.fasta" -t 500 -s "normal_SPAdes"
mv "${OUTDATADIR}/${sample_name}/Assembly/contigs.fasta.TRIMMED.fasta" "${OUTDATADIR}/${sample_name}/Assembly/${sample_name}_contigs_trimmed.fasta"


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
"${shareScript}/run_ANI_REFSEQ.sh" -n "${sample_name}" -p "${project}" -c "${config}"
# Get end time of ANI and calculate run time and append to time summary (and sum to total time used
end=$SECONDS
timeANI=$((end - start))
echo "autoANI - ${timeANI} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeANI))

### Average Nucleotide Identity ###
echo "----- Running ANI for Species confirmation -----"
# ANI uses assembly and sample would have exited already if assembly did not complete, so no need to check
# Get start time of ANI
start=$SECONDS
# run ANI
"${shareScript}/run_ANI.sh" -n "${sample_name}" -g "${genus}" -s "${species}" -p "${project}" -c "${config}"
# Get end time of ANI and calculate run time and append to time summary (and sum to total time used
end=$SECONDS
timeANIREF=$((end - start))
echo "ANIREF - ${timeANI} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeANIREF))

# Get taxonomy from currently available files (Only ANI, has not been run...yet, will change after discussions)
"${shareScript}/determine_taxID.sh" -n "${sample_name}" -p "${project}" -c "${config}"

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
"${shareScript}/run_c-sstar.sh" -n "${sample_name}" -g "${csstar_gapping}" -s "${csstar_identity}" -p "${project}" -c "${config}"
"${shareScript}/run_c-sstar.sh" -n "${sample_name}" -g "${csstar_gapping}" -s "${csstar_identity}" -p "${project}" -c "${config}" -d "${local_DBs}/star/ResGANNOT_20180608_srst2.fasta"

### GAMA - finding AR Genes ###
echo "----- Running GAMA for AR Gene identification -----"
${shareScript}/run_GAMA.sh -n "${sample_name}" -p "${project}" -c "${config}"

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
	"${shareScript}/run_MLST.sh" -n "${sample_name}" -p "${project}" "-f" "abaumannii" -c "${config}"
	python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst" -t standard
	mv "${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst" "${processed}/${project}/${sample_name}/MLST/${sample_name}_Oxford.mlst"
	mv "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" "${processed}/${project}/${sample_name}/MLST/${sample_name}_Pasteur.mlst"
	#Check for "-", unidentified type
	type1=$(tail -n1 ${processed}/${project}/${sample_name}/MLST/${sample_name}_abaumannii.mlst | cut -d' ' -f3)
	type2=$(head -n1 ${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst | cut -d' ' -f3)
	if [[ "${type1}" = "-" ]]; then
		"${shareScript}/run_srst2_mlst.sh" -n "${sample_name}" -p "${project}" -g "Acinetobacter" -s "baumannii#1" -c "${config}"
		python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${sample_name}/MLST/${sample_name}_srst2_acinetobacter_baumannii-baumannii#1.mlst" -t srst2
	fi
	if [[ "${type2}" = "-" ]]; then
		"${shareScript}/run_srst2_mlst.sh" -n "${sample_name}" -p "${project}" -g "Acinetobacter" -s "baumannii#2"  -c "${config}"
		python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${sample_name}/MLST/${sample_name}_srst2_acinetobacter_baumannii-baumannii#2.mlst" -t srst2
	fi
elif [[ "${genus}_${species}" = "Escherichia_coli" ]]; then
	# Verify that ecoli_2 is default and change accordingly
	"${shareScript}/run_MLST.sh" -n "${sample_name}" -p "${project}" "-f" "ecoli_2" -c "${config}"
	python3 "${shareScript}/check_and_fix_MLST.py" -i "${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst" -t standard
	mv "${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst" "${processed}/${project}/${sample_name}/MLST/${sample_name}_Pasteur.mlst"
	mv "${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst" "${processed}/${project}/${sample_name}/MLST/${sample_name}_Achtman.mlst"
	type2=$(tail -n1 ${processed}/${project}/${sample_name}/MLST/${sample_name}_ecoli_2.mlst | cut -d' ' -f3)
	type1=$(head -n1 ${processed}/${project}/${sample_name}/MLST/${sample_name}.mlst | cut -d' ' -f3)
	if [[ "${type1}" = "-" ]]; then
		"${shareScript}/run_srst2_mlst.sh" -n "${sample_name}" -p "${project}" -g "Escherichia" -s "coli#1" -c "${config}"
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
	${shareScript}/run_c-sstar.sh -n "${sample_name}" -g g -s o -p "${project}" -l  -c "${config}"
	${shareScript}/run_plasmidFinder.sh -n "${sample_name}" -p "${project}" -o plasmidFinder_on_plasFlow -c "${config}"
	${shareScript}/run_GAMA.sh -n "${sample_name}" -p "${project}" -l -c "${config}"
	end=$SECONDS
	timeplasflow=$((end - start))
	echo "plasmidFlow - ${timeplasflow} seconds" >> "${time_summary_redo}"
	totaltime=$((totaltime + timeplasflow))
fi

"${shareScript}/validate_piperun.sh" -n "${sample_name}" -p "${project}"  -c "${config}" > "${processed}/${project}/${sample_name}/${sample_name}_pipeline_stats.txt"

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
ml -Python3/3.5.2
