#!/bin/sh -l

#$ -o best_hit_from_gottcha1.out
#$ -e best_hit_from_gottcha1.err
#$ -N best_hit_from_gottcha1
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Grabs the best species match based on read hits (not relative abundance) from the gottcha tool run
#
# Usage: ./best_hit_from_gottcha1.sh sample_name run_ID
#
# Output location: default_config.sh_output_location/run_ID/sample_name/gottcha/
#
# Modules required: None
#
# V1.0
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [ -z "$1" ]; then
	echo "Empty sample name supplied to $0, exiting"
	exit 1
elif [ -z "$2" ]; then
	echo "Empty project_id supplied to $0, exiting"
	exit 1
# command line version of usage for script
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./best_hit_from_gottcha1.sh   sample_name   run_ID"
	echo "Output is saved to ${processed}/miseq_run_ID/sample_name/gottcha/"
	exit 0
fi

#Sets output directory to the gottcha folder of sample ID that was passed. Sample ID folder is found under processed samples in MMB_Data folder
OUTDATADIR="${processed}/${2}/${1}/gottcha"

#Sets many of the values of percents/descriptions to default values to guard against null/unset values if a failure occurs during processing
unclass_percent=0
domain_percent=100
phylum_percent=0
class_percent=0
order_percent=0
family_percent=0
genus_percent=0
species_percent=0
# Bacteria is set as domain due to gottcha v1 not classifying any reads higher than Phylum level
domain="bacteria"
phylum="N/A"
class="N/A"
order="N/A"
family="N/A"
genus="N/A"
species="N/A"
sum_reads=0
qc_reads=0
phylum_reads=0
class_reads=0
order_reads=0
family_reads=0
genus_reads=0
species_reads=0

# If the source TSV produced from gottcha does not exist then the script will exit with error code 1
if [[ ! -s "${OUTDATADIR}/gottcha_S/${1}.gottcha.tsv" ]]; then
	echo "GOTTCHA output tsv does not exist, exiting..."
	exit 1
fi

# Process the TSV output of gottcha (species level analysis) line by line
while IFS='	' read -r -a line  || [ -n "$line" ]; do
	# Convert the perfect match to proper format from 1.00 to 100
	if [[ "${line[2]}" = "1.0000" ]]; then
		percent=100
	# Convert all non-perfect matches to the correct matching percent values
	else
		percent="${line[2]:2:2}.${line[2]:4:2}"
	fi
	# Convert a no-match to the correct percent value
	if [[ "${percent}" = "00.00" ]]; then
		percent=0
	fi
	# Takes the first letter of the first column as shorthand for identifying the taxonomic level
	classification="${line[0]::1}"
	# Pulls the taxa information from the second column
	description="${line[1]}"
	# Pulls the reads/hits matching the current taxa from the 8th column
	reads="${line[7]}"
	# Acts to properly categorize each lines info based on the classifcation
	# Does not exist in GOTTCHA but kept for later if ever introduced, Place-holder for all unclassified reads
	if [ "${classification}" = "u" ]; then
		unclass_percent="${percent}"
	# Is not classified in GOTTCHA but kept for later if ever introduced, Place-holder for all domain level info
	elif [ "${classification}" = "d" ] && [ "${domain}" = "N/A" ]; then
		domain="${description}"
		domain_percent="${percent}"
	# Phylum level description, percent, and reads are saved if the current reads are greater than the stored reads for phylum level.
	# In all instances of phylum level classification, reads are totalled to provide the classified reads (sum_reads) value for percent calculations
	elif [ "${classification}" = "p" ]; then
		if [ $phylum = "N/A" ] || [ "${reads}" -gt "${phylum_reads}" ]; then
			phylum="${description}"
			phylum_percent="${percent}"
			sum_reads=$(( sum_reads + reads ))
			phylum_reads="${reads}"
		else
			sum_reads=$(( sum_reads + reads ))
		fi
	# If current line is classified as class level and there are more reads than the currently saved value, the stored values are replaced with the current lines info
	elif [ "${classification}" = "c" ]; then
		if [ "${class}" = "N/A" ] || [ "${reads}" -gt "${class_reads}" ]; then
			class="${description}"
			class_percent="${percent}"
			class_reads="${reads}"
		fi
	# If current line is classified as order level and there are more reads than the currently saved value, the stored values are replaced with the current lines info
	elif [ "${classification}" = "o" ]; then
		if [ "${order}" = "N/A" ] || [ "${reads}" -gt "${order_reads}" ]; then
			order="${description}"
			order_percent="${percent}"
			order_reads="${reads}"
		fi
	# If current line is classified as family level and there are more reads than the currently saved value, the stored values are replaced with the current lines info
	elif [ "${classification}" = "f" ]; then
		if [ "${family}" = "N/A" ] || [ "${reads}" -gt "${family_reads}" ]; then
			family="${description}"
			family_percent="${percent}"
			family_reads="${reads}"
		fi
	# If current line is classified as genus level and there are more reads than the currently saved value, the stored values are replaced with the current lines info
	elif [ "${classification}" = "g" ]; then
		if [ "${genus}" = "N/A" ] || [ "${reads}" -gt "${genus_reads}" ]; then
			genus="${description}"
			genus_percent="${percent}"
			genus_reads="${reads}"
		fi
	# If current line is classified as species level and there are more reads than the currently saved value, the stored values are replaced with the current lines info
	elif [ "${classification}" = "s" ]; then
		if [ "${species}" = "N/A" ] || [ "${reads}" -gt "${species_reads}" ]; then
			gs="($description)"			species=${gs[@]:1:-1}
			species_percent="${percent}"
			species_reads="${reads}"
		fi
	fi
done < "${OUTDATADIR}/gottcha_S/${1}.gottcha.tsv"

# Calculate % of unclassified reads using sum of highest taxon level reads against total reads found in QC counts
# Checks for the existence of the trimmed_counts file for the sample. If found true % classification is performed
if [[ -s "${processed}/${2}/${1}/preQCcounts/${1}_trimmed_counts.txt" ]]; then
	# The total reads is determined by pulling from the trimmed_counts
	qc_reads=$(tail -n 1 "${processed}/${2}/${1}/preQCcounts/${1}_trimmed_counts.txt" | cut -d'	' -f13)
	qc_reads=$((qc_reads/2))
	# Unclassified reads are calculated by subtracting the sum of all phylum reads from the total possible reads
	unclass_reads=$(( qc_reads - sum_reads ))
	# The unclassified percent is found using unclassified reads divided by the total possible reads
	u_percent=$( echo "${unclass_reads} ${qc_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	# Each taxas % applied reads is calculated using the applied reads divided by the total possible reads
	domain_percent_total=$( echo "${u_percent}" | awk '{ printf "%2.2f", 100-$1}' )
	phylum_percent_total=$( echo "${phylum_reads} ${qc_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	class_percent_total=$( echo "${class_reads} ${qc_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	order_percent_total=$( echo "${order_reads} ${qc_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	family_percent_total=$( echo "${family_reads} ${qc_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	genus_percent_total=$( echo "${genus_reads} ${qc_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	species_percent_total=$( echo "${species_reads} ${qc_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
# If the trimmed_counts file is not found all true % applied percentages are lableled as unknown
else
	echo "Cant find ${processed}/${2}/${1}/preQCcounts/${1}_trimmed_counts.txt. Continuing without true % reads used..."
	unclass_percent="UNK"
	u_percent="UNK"
	unclass_reads="UNK"
	domain_percent_total="UNK"
	phylum_percent_total="UNK"
	class_percent_total="UNK"
	order_percent_total="UNK"
	family_percent_total="UNK"
	genus_percent_total="UNK"
	species_percent_total="UNK"
fi

#Prints the final best taxon identification output (using species output from gottcha) to the gottcha directory for the sample
echo -e "U: ${unclass_percent} (${u_percent}%-${unclass_reads}) unclassified\\nD: ${domain_percent} (${domain_percent_total}) ${domain}\\nP: ${phylum_percent} (${phylum_percent_total}) ${phylum}\\nC: ${class_percent} (${class_percent_total}) ${class}\\nO: ${order_percent} (${order_percent_total}) ${order}\\nF: ${family_percent} (${family_percent_total}) ${family}\\nG: ${genus_percent} (${genus_percent_total}) ${genus}\\nS: ${species_percent} (${species_percent_total}) ${species}" > "${OUTDATADIR}/${1}_gottcha_species_summary.txt"

#Script exited gracefully (unless something else inside failed...working on finding all possible failure points!)
exit 0
