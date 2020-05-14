#!/bin/sh -l

#$ -o best_hit_kraken.out
#$ -e best_hit_kraken.err
#$ -N best_hit_kraken
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Grabs the best species match based on %/read hits from the kraken tool run
#
# Usage: ./best_hit_from_kraken.sh sample_name pre/post(relative to assembly) paired|assembled(source_type) run_ID kraken|kraken2(source) [alt_path-for-output]
#
# Output location: default_config.sh_output_location/run_ID/sample_name/kraken/pre|post-Assembly/
#				or				alt_path-for-output/kraken/pre|postAssembly
# Modules required: None
#
# v1.0.1 (05/12/2020)
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
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./best_hit_from_kraken.sh  sample_name  pre|post paired|assembled run_ID	kraken|kraken2 [alt_path-for-output]"
	echo "Output is saved to ${processed}/miseq_run_ID_id/sample_name/kraken/pre|postassembly/sample_name_kraken_summary_paired|assembled"
	echo "				or				alt_path-for-output/kraken/pre|postAssembly"
	exit 0
elif [ -z "$2" ]; then
	echo "Empty assembly relativity supplied to $0, exiting"
	exit 1
elif [ -z "$3" ]; then
	echo "Empty source type supplied to $0, exiting"
	exit 1
elif [ -z "$4" ]; then
	echo "Empty project id supplied to $0, exiting"
	exit 1
elif [ -z "$5" ]; then
	echo "Empty out folder supplied to $0, exiting"
	exit 1
elif [ -z "${6}" ]; then
	echo "Using default path - ${processed}"
	OUTDATADIR="${processed}/${4}/${1}/${5}/${2}Assembly"
elif [ -d "${6}" ]; then
	echo "Using given path - ${3}"
	OUTDATADIR="${6}/${4}/${1}/${5}/${2}Assembly"
else
	OUTDATADIR="${processed}/${4}/${1}/${5}/${2}Assembly"
	echo "Path does not exist - ${6}, using default (${processed})"
	exit 1
fi

# Checks what source flag was set as,indicating that it was from kraken or kraken2
if [[ "${5}" != "kraken" ]] && [[ "${5}" != "kraken2" ]]; then
	echo "Source must be kraken or kraken2, exiting"
	exit 334
else
	source="${5}"
fi

#Sets output folder to the correct path relative to assembly completion
OUTDATADIR="${processed}/${4}/${1}/${5}/${2}Assembly"
echo "-${OUTDATADIR}-"

#Creates the default values for output in case any calculations are interrupted.
unclass_reads=0
class_reads=0
u_percent=0
unclass_percent=0
domain_percent=0
phylum_percent=0
class_percent=0
order_percent=0
family_percent=0
genus_percent=0
species_percent=0
domain_percent_total=0
phylum_percent_total=0
class_percent_total=0
order_percent_total=0
family_percent_total=0
genus_percent_total=0
species_percent_total=0
domain_reads=0
phylum_reads=0
class_reads=0
order_reads=0
family_reads=0
genus_reads=0
species_reads=0
classified_reads=0
domain="N/A"
phylum="N/A"
class="N/A"
order="N/A"
family="N/A"
genus="N/A"
species="N/A"

#Checks to see if the list file used for calculations exists and exits if it does not
if [[ ! -s "${OUTDATADIR}/${1}_${3}.list" ]]; then
	echo "${OUTDATADIR}/${1}_${3}.list does not exist"
	exit 1
fi

#Parses the kraken output list line by line
while IFS= read -r line  || [ -n "$line" ]; do
	# Removes the whitespace at the end of the line
	line=${line,,}
	# Turns the line into an array delimited by whitespace
	arrLine=(${line})
	# First element in array is the percent of reads identified as the current taxa
	percent=${arrLine[0]}
	# 3rd element is the taxon level classification
	classification=${arrLine[3]}
	# 2nd element is the actual read count identified as the current taxa
	reads=${arrLine[1]}
	# All elements beyond 5 are stored as the description of the taxa, genus/species/strain and any other included information
	description="${arrLine[*]:5}"
	# Pulls the first word from the description to find if it is root, the descriptor used to identify any number of reads classified, but higher than domain
	first_desc=$(echo "${description}" | cut -d' ' -f1)
	# Unclassified read info (reads and percent) is captured if the classifier claims to be 'u'
	# Kraken lists entries in order of abundance (highest first) therefore each time a new taxon level is encountered it will be the highest ranking representative

	#echo "${reads} ${domain_reads}"
	if [ "${classification}" = "u" ]; then
		unclass_percent=${percent}
		# Will be overwritten by pre flag, but not in post
		unclass_reads=${reads}
		# Grabs all supra domain level read info (reads and percent)
	elif [ "${classification}" = "-" ] && [ "${first_desc}" = "root" ]; then
		classified_reads="${reads}"
	# Grabs all read info (identifier, reads and percent) for best domain level entry
	elif [ "${classification}" = "d" ] && [ "${reads}" -gt "${domain_reads}" ]; then
		domain=${description^}
		domain_percent=${percent}
		domain_reads=${reads}
	# Grabs all read info (identifier, reads and percent) for best phylum level entry
	elif [ "${classification}" = "p" ] && [ "${reads}" -gt "${phylum_reads}" ]; then
		phylum=${description^}
		phylum_percent=${percent}
		phylum_reads=${reads}
	# Grabs all read info (identifier, reads and percent) for best class level entry
	elif [ "${classification}" = "c" ] && [ "${reads}" -gt "${class_reads}" ]; then
		class=${description^}
		class_percent=${percent}
		class_reads=${reads}
	# Grabs all read info (identifier, reads and percent) for best order level entry
	elif [ "${classification}" = "o" ] && [ "${reads}" -gt "${order_reads}" ]; then
		order=${description^}
		order_percent=${percent}
		order_reads=${reads}
	# Grabs all read info (identifier, reads and percent) for best family level entry
	elif [ "${classification}" = "f" ] && [ "${reads}" -gt "${family_reads}" ]; then
		family=${description^}
		family_percent=${percent}
		family_reads=${reads}
	# Grabs all read info (identifier, reads and percent) for best genus level entry
	elif [ "${classification}" = "g" ] && [ "${reads}" -gt "${genus_reads}" ]; then
		genus=${description^}
		genus_percent=${percent}
		genus_reads=${reads}
	# Grabs all read info (identifier, reads and percent) for best species level entry
	elif [ "${classification}" = "s" ] && [ "${reads}" -gt "${species_reads}" ]; then

		echo "Old: ${species}-${species_reads}"
		gs=(${description^})
		species=${gs[@]:1}
		species_percent=${percent}
		species_reads=${reads}
		echo "New: ${species}-${species_reads}"
	fi
done < "${OUTDATADIR}/${1}_${3}.list"

# Calculate % of unclassified reads using sum of highest taxon level reads against total reads found in QC counts
# Grabs total possible reads from preQC counts if kraken was used on reads (pre assembly)
if [[ "${2}" = "pre" ]]; then
	# Checks for the existence of the preQC counts file to get total possible reads
	if [[ -s "${processed}/${4}/${1}/preQCcounts/${1}_trimmed_counts.txt" ]]; then
		# Pulls the total number of possible reads from the preQC counts file
		file_reads=$(tail -n 1 "${processed}/${4}/${1}/preQCcounts/${1}_trimmed_counts.txt" | cut -d'	' -f13)
		# Calculates the true count of unclassified reads rather than the reported value from kraken
		unclass_reads=$(( file_reads - classified_reads ))
		# Calculates the percent of unclassified reads using the total possible reads
		u_percent=$(echo "${unclass_reads} ${file_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	#The file does not exist and something terrible went wrong early in processing of this sample
	else
		# Assigns the true unclassified read count as unknown
		u_percent="UNK"
	fi
# Grabs total possible bases from contigs in trimmed assembly (post assembly, using weighted kraken output)
elif [[ "${2}" = "post" ]]; then
	total_reads=$(( classified_reads + unclass_reads ))
	# Checks for existence of trimmed contig file
	if [[ -s "${processed}/${4}/${1}/Assembly/${1}_scaffolds_trimmed.fasta" ]]; then
		# Sums unclassified and classified weighted base length from contigs
		file_reads=$(( unclass_reads + classified_reads ))
		# Calculates percent of classified reads as 100*classified reads/contigs
		u_percent=$(echo "${unclass_reads} ${file_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	# trimmed contigs file was not found. SPAdes may have failed (along with kraken) and true % classified will not be known
	else
		u_percent="UNK"
	fi
# Some other impossible relative assembly point was used (not pre or post) and so all values will be set to unknown
else
	unclass_percent="UNK"
fi

#If true unclassified percent is unknown then all true classified values for each taxa will not be calculated
if [ "${u_percent}" = "UNK" ]; then
	domain_percent_total="UNK"
	phylum_percent_total="UNK"
	class_percent_total="UNK"
	order_percent_total="UNK"
	family_percent_total="UNK"
	genus_percent_total="UNK"
	species_percent_total="UNK"
# All true classified values are calculated using 100*taxa_reads/total_reads
else
	domain_percent_total=$( echo "${domain_reads} ${file_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	phylum_percent_total=$( echo "${phylum_reads} ${file_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	class_percent_total=$( echo "${class_reads} ${file_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	order_percent_total=$( echo "${order_reads} ${file_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	family_percent_total=$( echo "${family_reads} ${file_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	genus_percent_total=$( echo "${genus_reads} ${file_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
	species_percent_total=$( echo "${species_reads} ${file_reads}" | awk '{ printf "%2.2f", ($1*100)/$2 }' )
fi

#Print out the best taxa for each level and its corresponding % of reads reported by kraken, % reads of total, taxon description
(echo -e "U: ${unclass_percent} (${u_percent}) unclassified\\nD: ${domain_percent} (${domain_percent_total}) ${domain}\\nP: ${phylum_percent} (${phylum_percent_total}) ${phylum}\\nC: ${class_percent} (${class_percent_total}) ${class}\\nO: ${order_percent} (${order_percent_total}) ${order}\\nF: ${family_percent} (${family_percent_total}) ${family}\\nG: ${genus_percent} (${genus_percent_total}) ${genus}\\nS: ${species_percent} (${species_percent_total}) ${species}") > "${OUTDATADIR}/${1}_kraken_summary_${3}.txt"

#Script exited gracefully (unless something else inside failed)
exit 0
