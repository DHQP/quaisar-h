#!/bin/bash -l

#$ -o quaisar_X.out
#$ -e quaisar_X.err
#$ -N quasX
#$ -pe smp 14
#$ -cwd
#$ -q all.q

# Copy config file into working directory to allow changes to made to output directory if necessary
shareScript=$(pwd)
echo "${shareScript}"
if [[ -f "${shareScript}/config_template.sh" ]]; then
	if [[ -f "${shareScript}/config.sh" ]]; then
		rm -r "${shareScript}/config.sh"
	fi
	echo "Trying to copy config_template.sh"
	cp "${shareScript}/config_template.sh" "${shareScript}/config.sh"
fi

#Import the config file with shortcuts and settings
. ${shareScript}/config.sh

#Import the module file that loads all necessary mods
. ${mod_changers}/pipeline_mods





#
# The straight pipeline that runs all the tools that have been designated as necessary (and some others that are typically run also)
#
# Usage ./quaisar.sh.sh isolate_name project_name output_directory_to_put_project/isolate_name
#  project/isolate_name must already have a populated FASTQs folder to work with
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to quaisar.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./quaisar.sh  sample_name miseq_run_id(or_project_name) optional_alternate_directory"
	echo "Populated FASTQs folder needs to be present in ${2}/${1}, wherever it resides"
	echo "Output by default is processed to ${processed}/miseq_run_id/sample_name"
	exit 0
elif [[ -z "{2}" ]]; then
	echo "No Project/Run_ID supplied to quaisar.sh, exiting"
	exit 33
fi

#Time tracker to gauge time used by each step
totaltime=0

# Set arguments to filename(sample name) project (miseq run id) and outdatadir(${processed}/project/filename)
filename="${1}"
project="${2}"
OUTDATADIR="${processed}/${2}"
if [[ ! -z "${3}" ]]; then
	OUTDATADIR="${3}/${2}"
fi

# Remove old run stats as the presence of the file indicates run completion
if [[ -f "${processed}/${proj}/${file}/${file}_pipeline_stats.txt" ]]; then
	rm "${processed}/${proj}/${file}/${file}_pipeline_stats.txt"
fi

# Create an empty time_summary file that tracks clock time of tools used
touch "${OUTDATADIR}/${filename}/${filename}_time_summary.txt"
time_summary=${OUTDATADIR}/${filename}/${filename}_time_summary.txt

echo "Time summary for ${project}/${filename}:" >> "${time_summary}"
echo "${project}/${filename} started at ${global_time}"

echo "Starting processing of ${project}/${filename}"
#Checks if FASTQ folder exists for current sample
if [[ -d "$OUTDATADIR/$filename/FASTQs" ]]; then
	# Checks if FASTQ folder contains any files then continue
	if [[ "$(ls -A "${OUTDATADIR}/${filename}/FASTQs")" ]]; then
		# Checks to see if those files in the folder are unzipped fastqs
		count_unzip=`ls -1 ${OUTDATADIR}/${filename}/FASTQs/*.fastq 2>/dev/null | wc -l`
		count_zip=`ls -1 ${OUTDATADIR}/${filename}/FASTQs/*.fastq.gz 2>/dev/null | wc -l`
		if [[ ${count_unzip} != 0 ]]; then
		#if [[ -f "${OUTDATADIR}/${filename}/FASTQs/${filename}"*".fastq" ]]; then
			echo "----- FASTQ(s) exist, continuing analysis -----"
			if [[ -f "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq" ]] && [[ ! -f "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq.gz" ]]; then
				gzip < "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq" > "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq.gz"
			fi
			if [[ -f "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq" ]] && [[ ! -f "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq.gz" ]]; then
				gzip < "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq" > "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq.gz"
			fi
		# Checks if they are zipped fastqs (checks for R1 first)
		elif [[ -f "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq.gz" ]]; then
			#echo "R1 zipped exists - unzipping"
			gunzip -c "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq.gz" > "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq"
			# Checks for paired R2 file
			if [[ -f "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq.gz" ]]; then
				#echo "R2 zipped exists - unzipping"
				gunzip -c "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq.gz" > "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq"
			else
				echo "No matching R2 to unzip :("
			fi
		# Checks to see if there is an abandoned R2 zipped fastq
		elif [[ -f "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq.gz" ]]; then
			#echo "R2 zipped  exists - unzipping"
			gunzip -c "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq.gz" > "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq"
			echo "No matching R1 to unzip :("
		fi
	# If the folder is empty then return from function
	else
		echo "FASTQs folder empty - No fastqs available for ${filename} (and download was not requested). Either unzip fastqs to $OUTDATADIR/FASTQs or run the -d flag to trigger unzipping of gzs"
		return 1
	fi
# If the fastq folder does not exist then return out of function
else
	echo "FASTQs not downloaded and FASTQs folder does not exist for ${filename}. No fastqs available (and download was not requested). Unzip fastqs to ${OUTDATADIR}/FASTQs"
	return 1
fi

# Get start time for qc check
start=$SECONDS
### Count the number of Q20, Q30, bases and reads within a pair of FASTQ files
echo "----- Counting read quality -----"
# Checks for and creates the specified output folder for the QC counts
if [ ! -d "$OUTDATADIR/$filename/preQCcounts" ]; then
	echo "Creating $OUTDATADIR/$filename/preQCcounts"
	mkdir -p "$OUTDATADIR/$filename/preQCcounts"
fi
# Run qc count check on raw reads
python "${shareScript}/Fastq_Quality_Printer.py" "${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq" "${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq" > "${OUTDATADIR}/$filename/preQCcounts/${filename}_counts.txt"

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
if [ ! -d "$OUTDATADIR/$filename/removedAdapters" ]; then
	echo "Creating $OUTDATADIR/$filename/removedAdapters"
	mkdir -p "$OUTDATADIR/$filename/removedAdapters"
# It complains if a folder already exists, so the current one is removed (shouldnt happen anymore as each analysis will move old runs to new folder)
else
	echo "Removing old $OUTDATADIR/$filename/removedAdapters"
	rm -r "$OUTDATADIR/$filename/removedAdapters"
	echo "Recreating $OUTDATADIR/$filename/removedAdapters"
	mkdir -p "$OUTDATADIR/$filename/removedAdapters"
fi

# Run bbduk
bbduk.sh -"${bbduk_mem}" threads="${procs}" in="${OUTDATADIR}/${filename}/FASTQs/${filename}_R1_001.fastq" in2="${OUTDATADIR}/${filename}/FASTQs/${filename}_R2_001.fastq" out="${OUTDATADIR}/${filename}/removedAdapters/${filename}-noPhiX-R1.fsq" out2="${OUTDATADIR}/${filename}/removedAdapters/${filename}-noPhiX-R2.fsq" ref="${phiX_location}" k="${bbduk_k}" hdist="${bbduk_hdist}"
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
if [ ! -d "$OUTDATADIR/$filename/trimmed" ]; then
	mkdir -p "$OUTDATADIR/$filename/trimmed"
fi
# Run trimmomatic
trimmomatic "${trim_endtype}" -"${trim_phred}" -threads "${procs}" "${OUTDATADIR}/${filename}/removedAdapters/${filename}-noPhiX-R1.fsq" "${OUTDATADIR}/${filename}/removedAdapters/${filename}-noPhiX-R2.fsq" "${OUTDATADIR}/${filename}/trimmed/${filename}_R1_001.paired.fq" "${OUTDATADIR}/${filename}/trimmed/${filename}_R1_001.unpaired.fq" "${OUTDATADIR}/${filename}/trimmed/${filename}_R2_001.paired.fq" "${OUTDATADIR}/${filename}/trimmed/${filename}_R2_001.unpaired.fq" LLUMINACLIP:"${trim_adapter_location}:${trim_seed_mismatch}:${trim_palindrome_clip_threshold}:${trim_simple_clip_threshold}:${trim_min_adapt_length}:${trim_complete_palindrome}" SLIDINGWINDOW:"${trim_window_size}:${trim_window_qual}" LEADING:"${trim_leading}" TRAILING:"${trim_trailing}" MINLEN:"${trim_min_length}"
# Get end time of trimmomatic and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeTrim=$((end - start))
echo "Trimming - ${timeTrim} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeTrim))


# Check differences after QC and trimming (also for gottcha proper read count for assessing unclassified reads)
# Get start time for qc check on trimmed reads
start=$SECONDS
### Count the number of Q20, Q30, bases and reads within the trimmed pair of FASTQ files
echo "----- Counting read quality of trimmed files-----"
# Checks for and creates the specified output folder for the QC counts
if [ ! -d "$OUTDATADIR/$filename/preQCcounts" ]; then
	echo "Creating $OUTDATADIR/$filename/preQCcounts"
	mkdir -p "$OUTDATADIR/$filename/preQCcounts"
fi
# Run qc count check on filtered reads
python "${shareScript}/Fastq_Quality_Printer.py" "${OUTDATADIR}/${filename}/trimmed/${filename}_R1_001.paired.fq" "${OUTDATADIR}/${filename}/trimmed/${filename}_R2_001.paired.fq" > "${OUTDATADIR}/${filename}/preQCcounts/${filename}_trimmed_counts.txt"

# Merge both unpaired fq files into one for GOTTCHA
cat "${OUTDATADIR}/${filename}/trimmed/${filename}_R1_001.paired.fq" "${OUTDATADIR}/${filename}/trimmed/${filename}_R2_001.paired.fq" > "${OUTDATADIR}/${filename}/trimmed/${filename}.paired.fq"
cat "${OUTDATADIR}/${filename}/trimmed/${filename}_R1_001.unpaired.fq" "${OUTDATADIR}/${filename}/trimmed/${filename}_R2_001.unpaired.fq" > "${OUTDATADIR}/${filename}/trimmed/${filename}.single.fq"


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
"${shareScript}/run_kraken.sh" "${filename}" pre paired "${project}"
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
"${shareScript}/run_gottcha.sh" "${filename}" "${project}"
# Get end time of qc count and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeGott=$((end - start))
echo "Gottcha - ${timeGott} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeGott))

# Check reads using SRST2
echo "----- Running SRST2 -----"
start=$SECONDS
"${shareScript}/run_srst2_on_singleDB.sh" "${filename}" "${project}"
"${shareScript}/run_srst2_on_singleDB_alternateDB.sh" "${filename}" "${project}" "${local_DBs}/star/ResGANNOT_20180608_srst2.fasta"
end=$SECONDS
timesrst2=$((end - start))
echo "SRST2 - ${timesrst2} seconds" >> "${time_summary}"
totaltime=$((totaltime + timesrst2))

######  Assembling Using SPAdes  ######
echo "----- Assembling Using SPAdes -----"
# Get start time of SPAdes
start=$SECONDS
# script tries 3 times for a completed assembly
for i in 1 2 3
do
	# If assembly exists already and this is the first attempt (then the previous run will be used) [should not happen anymore as old runs are now renamed]
	if [ -s "${OUTDATADIR}/${filename}/Assembly/scaffolds.fasta" ]; then
		echo "Previous assembly already exists, using it (delete/rename the assembly folder at ${OUTDATADIR}/ if you'd like to try to reassemble"
	# Run normal mode if no assembly file was found
	else
		"${shareScript}/run_SPAdes.sh" "${filename}" normal "${project}"
	fi
	# Removes any core dump files (Occured often during testing and tweaking of memory parameter
	if [ -n "$(find "${shareScript}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
		echo "Found core dump files in assembly (assumed to be from SPAdes, but could be anything before that as well) and attempting to delete"
		find "${shareScript}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
	fi
done
# Returns if all 3 assembly attempts fail
if [[ -f "${OUTDATADIR}/${filename}/Assembly/scaffolds.fasta" ]] && [[ -s "${OUTDATADIR}/${filename}/Assembly/scaffolds.fasta" ]]; then
	echo "Assembly completed and created a non-empty scaffolds file"
else
	echo "Assembly FAILED 3 times, continuing to next sample..." >&2
	return 1
fi
for i in 1 2 3
do
	# If plasmid Assembly exists already and this is the first attempt (then the previous run will be used) [should not happen anymore as old runs are now renamed]
	if [ -f "${OUTDATADIR}/${filename}/plasmidAssembly/scaffolds.fasta" ]; then
		echo "Previous plasmid assembly already exists, using it (delete/rename the assembly folder at ${OUTDATADIR}/ if you'd like to try to reassemble"
	else
		"${shareScript}/run_SPAdes.sh" "${filename}" plasmid "${project}"
		# Removes any core dump files created from SPAdes (occurred fairly regularly during testing/tweaking)
		if [ -n "$(find "${shareScript}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
			echo "Found core dump files in plasmid Assembly (assumed to be from SPAdes, but could be anything before that as well) and attempting to delete"
			find "${shareScript}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
			else
			echo "No core dump files during plasmid assembly attempt number ${i}, plasmidSPAdes finished successfully and found nothing, creating dummy scaffolds file"
			>> "${OUTDATADIR}/${filename}/plasmidAssembly/scaffolds.fasta"
		fi
	fi
done
# Returns if all 3 assembly attempts fail
if [ ! -f "${OUTDATADIR}/${filename}/plasmidAssembly/scaffolds.fasta" ]; then
	echo "plasmid Assembly FAILED 3 times, continuing to next step..." >&2
fi
# Get end time of SPAdes and calculate run time and append to time summary (and sum to total time used)
end=$SECONDS
timeSPAdes=$((end - start))
echo "SPAdes - ${timeSPAdes} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeSPAdes))

### Removing Short Contigs  ###
echo "----- Removing Short Contigs -----"
python "${shareScript}/removeShortContigs.py" "${OUTDATADIR}/${filename}/Assembly/scaffolds.fasta" 500
mv "${OUTDATADIR}/${filename}/Assembly/scaffolds.fasta.TRIMMED.fasta" "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed.fasta"

### Removing Short Contigs  ###
echo "----- Removing Short plasmid Contigs -----"
python "${shareScript}/removeShortContigs.py" "${OUTDATADIR}/${filename}/plasmidAssembly/scaffolds.fasta" 500
mv "${OUTDATADIR}/${filename}/plasmidAssembly/scaffolds.fasta.TRIMMED.fasta" "${OUTDATADIR}/${filename}/plasmidAssembly/${filename}_plasmid_scaffolds_trimmed.fasta"

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
"${shareScript}/16s_blast.sh" "${filename}" "${project}"
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
	while IFS= read -r line;
	do
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
python3 "${shareScript}/fasta_headers.py" "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed_original.fasta" "${OUTDATADIR}/${filename}/Assembly/${filename}_scaffolds_trimmed.fasta"
if [[ -s "${OUTDATADIR}/${filename}/plasmidAssembly/${filename}_plasmid_scaffolds_trimmed.fasta" ]]; then
	mv "${OUTDATADIR}/${filename}/plasmidAssembly/${filename}_plasmid_scaffolds_trimmed.fasta" "${OUTDATADIR}/${filename}/plasmidAssembly/${filename}_plasmid_scaffolds_trimmed_original.fasta"
	python3 "${shareScript}/fasta_headers.py" "${OUTDATADIR}/${filename}/plasmidAssembly/${filename}_plasmid_scaffolds_trimmed_original.fasta" "${OUTDATADIR}/${filename}/plasmidAssembly/${filename}_plasmid_scaffolds_trimmed.fasta"
fi

### Average Nucleotide Identity ###
echo "----- Running ANI for Species confirmation -----"
# ANI uses assembly and sample would have exited already if assembly did not complete, so no need to check
# Get start time of ANI
start=$SECONDS
# run ANI
# Temp fix for strange genera until we do vs ALL all the time.
if [[ "${genus}" = "Peptoclostridium" ]] || [[ "${genus}" = "Clostridioides" ]]; then
	genus="Clostridium"
fi
"${shareScript}/run_ANI.sh" "${filename}" "${genus}" "${species}" "${project}"
#"${shareScript}/run_ANI.sh" "${filename}" "All" "All" "${project}"
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
	buscoDB="bacteria_odb9"
	# Iterate through taxon levels (species to domain) and test if a match occurs to entry in database. If so, compare against it
	busco_found=0
	for tax in $species $genus $family $order $class $phylum $kingdom $domain
	do
		if [ -d "${local_DBs}/BUSCO/${tax,}_odb9" ]
		then
			buscoDB="${tax,}_odb9"
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
"${shareScript}/run_c-sstar_on_single.sh" "${filename}" "${csstar_gapping}" "${csstar_identity}" "${project}"
"${shareScript}/run_c-sstar_on_single_alternate_DB.sh" "${filename}" "${csstar_gapping}" "${csstar_identity}" "${project}" "${local_DBs}/star/ResGANNOT_20180608_srst2.fasta"
# Should the parameters be different when checking on plasmids specifically
"${shareScript}/run_c-sstar_on_single.sh" "${filename}" "${csstar_gapping}" "${csstar_plasmid_identity}" "${project}" "--plasmid"
"${shareScript}/run_c-sstar_on_single_alternate_DB.sh" "${filename}" "${csstar_gapping}" "${csstar_plasmid_identity}" "${project}" "--plasmid" "${local_DBs}/star/ResGANNOT_20180608_srst2.fasta"

# Get end time of csstar and calculate run time and append to time summary (and sum to total time used
end=$SECONDS
timestar=$((end - start))
echo "c-SSTAR - ${timestar} seconds" >> "${time_summary}"
totaltime=$((totaltime + timestar))

# Get MLST profile
echo "----- Running MLST -----"
start=$SECONDS
"${shareScript}/run_MLST.sh" "${filename}" "${project}"
if [[ "${genus}_${species}" = "Acinetobacter_baumannii" ]]; then
	"${shareScript}/run_MLST.sh" "${filename}" "${project}" "-f" "abaumannii"
elif [[ "${genus}_${species}" = "Escherichia_coli" ]]; then
	# Verify that ecoli_2 is default and change accordingly
	"${shareScript}/run_MLST.sh" "${filename}" "${project}" "-f" "ecoli_2"
fi
end=$SECONDS
timeMLST=$((end - start))
echo "MLST - ${timeMLST} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeMLST))

# Try to find any plasmids
echo "----- Identifying plasmids using plasmidFinder -----"
start=$SECONDS
"${shareScript}/run_plasmidFinder.sh" "${filename}" "${project}" "plasmid"
"${shareScript}/run_plasmidFinder.sh" "${filename}" "${project}" "plasmid_on_plasmidAssembly"
end=$SECONDS
timeplasfin=$((end - start))
echo "plasmidFinder - ${timeplasfin} seconds" >> "${time_summary}"
totaltime=$((totaltime + timeplasfin))

# Append total time to bottom of time summary
echo "Total time: ${totaltime} seconds" >> "${time_summary}"

# Designate end of this sample #
global_end_time=$(date "+%m-%d-%Y_at_%Hh_%Mm_%Ss")
echo "

				End of sample ${filename}
				completed at ${global_end_time}

"

# Extra dump cleanse in case anything else failed
	if [ -n "$(find "${shareScript}" -maxdepth 1 -name 'core.*' -print -quit)" ]; then
		echo "Found core dump files at end of processing ${filename} and attempting to delete"
		find "${shareScript}" -maxdepth 1 -name 'core.*' -exec rm -f {} \;
	fi

	"${shareScript}/sample_cleaner.sh" "${file}" "${proj}"
	"${shareScript}/validate_piperun.sh" "${file}" "${proj}" > "${processed}/${proj}/${file}/${file}_pipeline_stats.txt"