#!/bin/sh -l

#$ -o run_sum.out
#$ -e run_sum.err
#$ -N run_sum
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Parses through summary file for the run and prints out a one word status for each sample in the run
#
# Usage: ./run_sum.sh run_ID
#
# Output location: standard out
#
# Modules required: None
#
# v1.0.2 (05/05/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
elif [[ -z "${1}" ]]; then
	echo "Empty run_ID supplied to run_sum.sh, exiting"
	exit 1
elif [[ "${1}" = "-h" ]]; then
	echo "Usage is ./run_sum.sh miseq_run_ID -l(optional, if using a list instead of project)"
	echo "Output is saved to ${processed}/miseq_run_ID"
	exit 0
fi

# Finds the newest summary file to use when parsing out info
if [[ -z "${2}" ]]; then
	proj="${1}"
	for summary in ${processed}/${proj}/*.sum; do
		third=$(echo "${summary}" | rev | cut -d'_' -f4 | rev)
		second=$(echo "${summary}" | rev | cut -d'_' -f5 | rev)
		first=$(echo "${summary}" | rev | cut -d'_' -f6 | rev)
		post_info=$(echo "${summary}" | rev | cut -d'_' -f1,2,3 | rev)
		pre_info=$(echo "${summary}" | rev | cut -d'_' -f7- | rev)
		#echo "Pre-${pre_info}"
		#echo "Post-${post_info}"
		#echo "1st-${first}"
		#echo "2nd-${second}"
		#echo "3rd-${third}"
		if [[ "${third}" = "201"* ]]; then
			new_name="${pre_info}_${third}_${first}_${second}_${post_info}"
			mv ${summary} ${new_name}
			echo "Tried moving ${summary} to ${new_name}"
		fi
	done
	sum_name=$(find ${processed}/${proj}/*.sum -maxdepth 1 -type f -printf '%p\n' | sort -rt '\0' -n | head -n 1)
	#sum_name=$(find ${processed}/${proj}/*.sum -maxdepth 1 -type f -printf '%h\0%d\0%p\n' | sort -rt '\0' -n)
	echo "${sum_name}"
	sum_name=$(basename ${sum_name})
	sum_file="${processed}/${1}/${sum_name}"
	if [[ -f "${processed}/${1}/${sum_name}" ]]; then
		echo -e "\n\nSummary found: ${sum_name}\n"
	else
		echo "No summary file detected, please execute run_sum.sh ${1} to create it"
	fi
else
	runsumdate=$(date "+%m_%d_%Y_at_%Hh_%Mm")
	echo "Creating run summary at ${runsumdate}"
	# Status of each individual sample is updated in its own folder and the run_summary file
	sum_name="${1}"
	if [[ -f "${1}" ]]; then
		sum_file="${1}"
		echo -e "\n\nSummary found: ${1}\n"
	else
		echo "Summary file not found: ${1}"
	fi
fi

#Full run overview status file is created for quick view
# The full run log is searched and single word status is returned for each sample in the log as SUCCESSFUL,WARNING, or FAILED
warning_samples=0
failed_samples=0
success_samples=0
warning_flags=""
warnings=0
failure_flags=""
failures=0
genus="unk"
species="unk"
sample_name=""
#echo "TEST-${proj}:::::${sum_name}"
while IFS= read -r var || [ -n "$var" ]; do
	#echo "${var}"
	first=$(echo "${var}" | cut -d' ' -f1)
	if [[ "${sample_name}" = "" ]]; then
		sample_name=$(echo "${var}" | cut -d' ' -f2 | cut -d'/' -f2)
	fi
	if [[ "${first}" = "----------" ]]; then
		sample_name=$(echo "${var}" | cut -d' ' -f2)
		status=$(echo "${var}" | cut -d' ' -f5)
		#echo "${sample_name}: ${status} ${failures} failures${failure_flags} and ${warnings} warnings${warning_flags}"
		printf "%-20s: %-7s : %25s : %s\\n" "${sample_name}" "${status}" "${genus} ${species}" "${failures} failures${failure_flags} and ${warnings} warnings${warning_flags}:${notes}"
		if [[ "${status}" == "WARNING" ]]; then
			outarray+=("${sample_name}: ${status}  ${warnings} warnings(${warning_flags})")
			imperfect_samples+=("${1}	${sample_name}	${genus:0:1}.${species}	${status}	${warning_flags:1}")
			warning_samples=$(( warning_samples + 1 ))
		elif [[ "${status}" == "FAILED" ]]; then
			outarray+=("${sample_name}: ${status} ${failures} failures${failure_flags} and ${warnings} warnings${warning_flags}")
			imperfect_samples+=("${1}	${sample_name}	${genus:0:1}.${species}	${status}	F${failure_flags};W${warning_flags}")
			failed_samples=$(( failed_samples + 1 ))
		else
			outarray+=("${sample_name}: ${status}")
			success_samples=$(( success_samples + 1 ))
		fi
		#echo "${outarray}"
		failure_flags=""
		warnings=0
		warning_flags=""
		failures=0
		notes=""
		continue
	fi

	#tool1=$(echo "${var}" | cut -d':' -f1)


	#echo "${var}"
	#echo "+${tool} --- ${tool_details1}"

	tool=$(echo "${var}" | cut -d':' -f1 | tr -d ' ')
	tool_status=$(echo "${var}" | cut -d':' -f2 | tr -d ' ')
	tool_details=$(echo "${var}" | cut -d':' -f3)
	if [[ "${tool}" = "Taxa" ]]; then
		species=$(echo "${tool_details}" | cut -d' ' -f3)
		genus=$(echo "${tool_details}" | cut -d' ' -f2)
	fi
	#echo "TOOL:${tool}:${tool_status}:${tool_details}"
	if [[ "${tool_status}" == "ALERT" ]]; then
		if [[ "${tool}" == "TIME" ]]; then
			if [[ "${notes}" != "" ]]; then
				notes="${notes},"
			fi
			notes="No time summary found"
		elif [[ "${tool}" == "Rawcoverage" ]]; then
			if [[ "${notes}" != "" ]]; then
				notes="${notes},"
			fi
			notes="${notes}Raw read coverage over 90x"
		elif [[ "${tool}" == "Trimmedcoverage" ]]; then
			if [[ "${notes}" != "" ]]; then
				notes="${notes},"
			fi
			notes="${notes}Trimmed coverage not in acceptable range(40-90)"
		elif [[ "${tool}" == "c-SSTAR" ]] || [[ "${tool}" == "c-sstar" ]]; then
			if [[ "${notes}" != "" ]]; then
				notes="${notes},"
			fi
			if [[ "${tool_details}" = *"NO KNOWN AMR genes"* ]]; then
				if [[ "${tool_details}" = *"(DB NOT up to date!"* ]]; then
					notes="${notes}No AMR found AND OLD ResGANNCBI DB used in c-sstar"
				else
					notes="${notes}No AMR found"
				fi
			else
				notes="${notes}OLD ResGANNCBI DB used in c-sstar"
			fi
		elif [[ "${tool}" == "c-SSTAR_plasFlow" ]] || [[ "${tool}" == "c-sstar_plasFlow" ]]; then
			if [[ "${tool_details}" = *"NO KNOWN AMR genes"* ]]; then
				if [[ "${notes}" != "" ]]; then
					notes="${notes},"
				fi
				if [[ "${tool_details}" = *"(DB NOT up to date!"* ]]; then
					notes="${notes}No plasmid AMR found AND OLD ResGANNCBI DB used in c-sstar plasmid"
				else
					notes="${notes}No plasmid AMR found"
				fi
			else
				notes="${notes}OLD ResGANNCBI DB used in c-sstar plasmid"
			fi
		elif [[ "${tool}" == "GAMA" ]]; then
			if [[ "${notes}" != "" ]]; then
				notes="${notes},"
			fi
			if [[ "${tool_details}" = *"NO KNOWN AMR genes"* ]]; then
				if [[ "${tool_details}" = *"(DB NOT up to date!"* ]]; then
					notes="${notes}No AMR found AND OLD ResGANNCBI DB used in GAMA"
				else
					notes="${notes}No AMR found"
				fi
			else
				notes="${notes}OLD ResGANNCBI DB used in GAMA"
			fi
		elif [[ "${tool}" == "GAMA_plasFlow" ]]; then
			if [[ "${tool_details}" = *"NO KNOWN AMR genes"* ]]; then
				if [[ "${notes}" != "" ]]; then
					notes="${notes},"
				fi
				if [[ "${tool_details}" = *"(DB NOT up to date!"* ]]; then
					notes="${notes}No plasmid AMR found AND OLD ResGANNCBI DB used in GAMA_plasFlow"
				else
					notes="${notes}No plasmid AMR found"
				fi
			else
				notes="${notes}OLD ResGANNCBI DB used in c-sstar plasFlow"
			fi
		elif [[ "${tool}" == "srst2" ]]; then
			if [[ "${notes}" != "" ]]; then
				notes="${notes},"
			fi
			if [[ "${tool_details}" = *"NO KNOWN AMR genes"* ]]; then
				if [[ "${tool_details}" = *"(DB NOT up to date!"* ]]; then
					notes="${notes}No srst2 AMR found AND OLD ResGANNCBI DB used in srst2"
				else
					notes="${notes}No srst2 AMR found"
				fi
			else
				notes="${notes}OLD ResGANNCBI DB used in srst2"
			fi
		elif [[ "${tool}" == "MLST-srst2" ]]; then
			if [[ "${notes}" != "" ]]; then
				notes="${notes},"
			fi
			notes="${notes}unexpected number of mlst-srst2 files"
		elif [[ "${tool}" == "ANI_REFSEQ" ]]; then
			if [[ "${notes}" != "" ]]; then
				notes="${notes},"
			fi
			notes="${notes}ANI_REFSEQ is out of date"
		elif [[ "${tool}" == "ANI_OSII" ]]; then
			if [[ "${notes}" != "" ]]; then
				notes="${notes},"
			fi
			notes="${notes}No coverage information extracted"
		elif [[ "${tool}" == "Trimming-R1" ]]; then
			if [[ -s "${processed}/${1}/${sample_name}/FASTQs/${1}_R1_001.fastq" ]]; then
				if [[ "${notes}" != "" ]]; then
					notes="${notes},"
				fi
				notes="${notes}trimmed R1 reads not zipped"
			else
				if [[ "${notes}" != "" ]]; then
					notes="${notes},"
				fi
				notes="${notes}No R1 bbduk counts for comparison"
			fi
		elif [[ "${tool}" == "Trimming-R2" ]]; then
			if [[ -s "${processed}/${1}/${sample_name}/FASTQs/${1}_R2_001.fastq" ]]; then
				if [[ "${notes}" != "" ]]; then
					notes="${notes},"
				fi
				notes="${notes}trimmed R2 reads not zipped"
			else
				if [[ "${notes}" != "" ]]; then
					notes="${notes},"
				fi
				notes="${notes}No R2 bbduk counts for comparison"
			fi
		fi
	elif [[ "${tool_status}" == "WARNING" ]]; then
		#echo "Found warning"
		if [[ "${tool}" == "FASTQs" ]]; then
			if [[ ! -s "${processed}/${1}/${sample_name}/FASTQs/${1}_R1_001.fastq" ]]; then
				warning_flags="${warning_flags}-Missing_R1_reads_file"
				warnings=$(( warnings + 1 ))
			elif [[ ! -s "${processed}/${1}/${sample_name}/FASTQs/${1}_R2_001.fastq" ]]; then
				warning_flags="${warning_flags}-Missing_R2_reads_file"
				warnings=$(( warnings + 1 ))
			fi
		elif [[ "${tool}" == "QCcounts" ]]; then
			warning_flags="${warning_flags}-Raw_read_count_below_1000000"
			warnings=$(( warning_flags + 1 ))
		elif [[ "${tool}" == "Trimming-R1" ]]; then
			warning_flags="${warning_flags}-No bases in trimmed R1 reads"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "Trimming-R2" ]]; then
			warning_flags="${warning_flags}-No bases in trimmed R2 reads"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "Q30_R1%" ]]; then
			warning_flags="${warning_flags}-Q30_R1_below_90%"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "Q30_R2%" ]]; then
			warning_flags="${warning_flags}-Q30_R2_below_70%"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "QCcountaftertrim" ]]; then
			warning_flags="${warning_flags}-Trimmed_reads_below_500000"
			warnings=$(( warnings + 1 ))
		#elif [[ "${tool}" == "Merging" ]]; then
		#	warning_flags="${warning_flags}-Missing_merged_reads"
		#	warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "PreClassify" ]]; then
			# Determine way to check for multiple species above 30%
			warning_flags="${warning_flags}-High_#_unclassified_reads_kraken(pre)"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "preClassContam." ]]; then
			warning_flags="${warning_flags}-Multiple_species_found_above_contamination_threshold(${contamination_threshold}%)_kraken(pre)"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "Gottcha_S" ]]; then
			if [[ ! -s "${processed}/${1}/${sample_name}/gottcha/gottcha_S/${sample_name}.gottcha_full.tsv" ]]; then
				warning_flags="${warning_flags}-Missing_tsv_file"
				warnings=$(( warnings + 1 ))
			elif [[ ! -s "${processed}/${1}/${sample_name}/gottcha/gottcha_S/${sample_name}_species.krona.html" ]]; then
				warning_flags="${warning_flags}-Missing_HTML_file"
				warnings=$(( warnings + 1 ))
			fi
		elif [[ "${tool}" == "GottchaV1Classifier" ]]; then
			# Determine way to check for multiple species above 30%
			warning_flags="${warning_flags}-High_#_unclassified_reads_gottcha"
			warnings=$(( warnings + 1 ))
		#elif [[ "${tool}" == "gottchaContam." ]]; then
		#	warning_flags="${warning_flags}-Multiple_species_found_above_contamination_threshold(${contamination_threshold}%)_gottcha"
		#	warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "ContigTrim" ]]; then
			warning_flags="${warning_flags}->200_contigs_remain"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "postClassify" ]]; then
			# Determine way to check for multiple species above 30%
			if [[ "${tool_details}" = *"possibly contaminated"* ]]; then
				warning_flags="${warning_flags}-Best_hit_below_50%(kraken_post),check_weighted_kraken"
			elif [[ "${tool_details}" = "unclassified reads"* ]]; then
				warning_flags="${warning_flags}-Unclassified_reads_above ${unclass_flag}(kraken_post)"
			fi
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "postClassContam." ]]; then
			warning_flags="${warning_flags}-Multiple_species_found_above_contamination_threshold(${contamination_threshold}%)_kraken(post)"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "weightedClassify" ]]; then
			# Determine way to check for multiple species above 30%
			warning_flags="${warning_flags}-High_#_unclassified_reads_kraken(weighted)"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "Taxa" ]]; then
			warning_flags="${warning_flags}-NO_species_was_determined"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "Assembly_Ratio" ]]; then
			warning_flags="${warning_flags}-Species_not_found_in_MMB_Bugs_database"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "MLST" ]]; then
			warning_flags="${warning_flags}-NOVEL_SCHEME"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "MLST-srst2" ]]; then
			warning_flags="${warning_flags}-ST_indeterminable"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "16s_best_hit" ]]; then
			warning_flags="${warning_flags}-NO_16s_best_species"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "16s_largest_hit" ]]; then
			warning_flags="${warning_flags}-NO_16s_largest_species"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "raw_read_count" ]]; then
			warning_flags="${warning_flags}-Low Raw read count (<1000000)"
			warnings=$(( warnings + 1 ))
		fi
	elif [[ "${tool_status}" == "FAILED" ]]; then
		#echo "Found failure-${tool}-${tool_details}"
		if [[ "${tool}" == "FASTQs" ]]; then
			failure_flags="${failure_flags}-NO_raw_Reads"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "QCcounts" ]]; then
			failure_flags="${failure_flags}-NO_pre_QC_counts.txt"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "Q30_R1%" ]]; then
			failure_flags="${failure_flags}-NO_pre_QC_counts.txt"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "Q30_R2%" ]]; then
			failure_flags="${failure_flags}-NO_pre_QC_counts.txt"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "Trimming" ]]; then
			failure_flags="${failure_flags}-NO_trimmed_reads"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "QCcountaftertrim" ]]; then
			failure_flags="${failure_flags}-NO_QC_trimmed_counts.txt"
			failures=$(( failures + 1 ))
		#elif [[ "${tool}" == "Merging" ]]; then
		#	failure_flags="${failure_flags}-Missing_merged_reads"
		#	failures=$(( failures + 1 ))
		elif [[ "${tool}" == "krona-kraken-preasmb" ]]; then
			if [[ ! -s "${processed}/${1}/${sample_name}/kraken/preAssembly/${1}_paired.krona" ]] && [[ ! -s "${processed}/${1}/${sample_name}/kraken/preAssembly/${1}_paired.html" ]]; then
				failure_flags="${failure_flags}-NO_HTML_or_krona_file(pre)"
				failures=$(( failures + 1 ))
			elif [[ ! -s "${processed}/${1}/${sample_name}/kraken/preAssembly/${1}_paired.html" ]]; then
				failure_flags="${failure_flags}-NO_HTML_file(pre)"
				failures=$(( failures + 1 ))
			elif [[ ! -s "${processed}/${1}/${sample_name}/kraken/preAssembly/${1}_paired.krona" ]]; then
				failure_flags="${failure_flags}-NO_krona_file(pre)"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-kraken_preclassify_failed"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "PreClassify" ]]; then
			if [[ -s "${processed}/${1}/${sample_name}/kraken/preAssembly/${1}_kraken_summary_paired.txt" ]]; then
				failure_flags="${failure_flags}-NO_classified_reads(pre)"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-NO_kraken_file(pre)"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "preClassContam." ]]; then
				failure_flags="${failure_flags}-NO_species_found_above_${contamination_threshold}%(post)"
				failures=$(( failures + 1 ))
		elif [[ "${tool}" == "Gottcha_S" ]]; then
			failure_flags="${failure_flags}-Missing_intermediate_gottcha_files (TSV AND HTML)"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "GottchaV1Classifier" ]]; then
			if [[ -s "${processed}/${1}/${sample_name}/gottcha/${1}_gottcha_species_summary.txt" ]]; then
				failure_flags="${failure_flags}-NO_classified_reads"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-NO_summary_gottcha"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "Assembly" ]]; then
			failure_flags="${failure_flags}-NO_assembly_file"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "plasmidAssembly" ]]; then
			failure_flags="${failure_flags}-NO_assembly_file"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "ContigTrim" ]]; then
			failure_flags="${failure_flags}-NO_trimmed_assembly_file"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "plasmidscontigTrim" ]]; then
			failure_flags="${failure_flags}-NO_trimmed_plasmid_assembly_file"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "krakenpostassembly" ]]; then
			failure_flags="${failure_flags}-NO_kraken_file(post)"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "krona-kraken-pstasmb" ]]; then
			if [[ ! -s "${processed}/${1}/${sample_name}/kraken/postAssembly/${1}_assembled.krona" ]] && [[ ! -s "${processed}/${1}/${sample_name}/kraken/postAssembly/${1}_assembled.html" ]]; then
				failure_flags="${failure_flags}-NO_HTML_or_krona_file(post)"
				failures=$(( failures + 1 ))
			elif [[ ! -s "${processed}/${1}/${sample_name}/kraken/postAssembly/${1}_assembled.html" ]]; then
				failure_flags="${failure_flags}-NO_HTML_file(post)"
				failures=$(( failures + 1 ))
			elif [[ ! -s "${processed}/${1}/${sample_name}/kraken/postAssembly/${1}_assembled.krona" ]]; then
				failure_flags="${failure_flags}-NO_krona_file(post)"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-kraken_postclassify_failed"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "postClassify" ]]; then
			if [[ -s "${processed}/${1}/${sample_name}/kraken/postAssembly/${1}_kraken_summary_assembled.txt" ]]; then
				failure_flags="${failure_flags}-NO_classified_contigs(post)"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-NO_kraken_file(post)"
				failures=$(( failures + 1 ))
			fi
		#elif [[ "${tool}" == "postClassContam." ]]; then
		#	failure_flags="${failure_flags}-Multiple_species_found_above_contamination_threshold(${contamination_threshold})_kraken(post)"
		#	failures=$(( failures + 1 ))
		elif [[ "${tool}" == "krakenweighted" ]]; then
			failure_flags="${failure_flags}-Top_match_under_50%_likely_contaminated"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "krona-kraken-weight" ]]; then
			if [[ ! -s "${processed}/${1}/${sample_name}/kraken/postAssembly/${1}_assembled_weighted.krona" ]] && [[ ! -s "${processed}/${1}/${sample_name}/kraken/postAssembly/${1}_assembled_weighted_BP_krona.html" ]]; then
				failure_flags="${failure_flags}-NO_HTML_or_krona_file(weighted)"
				failures=$(( failures + 1 ))
			elif [[ ! -s "${processed}/${1}/${sample_name}/kraken/postAssembly/${1}_assembled_weighted_BP_krona.html" ]]; then
				failure_flags="${failure_flags}-NO_HTML_file(weighted)"
				failures=$(( failures + 1 ))
			elif [[ ! -s "${processed}/${1}/${sample_name}/kraken/postAssembly/${1}_assembled_weighted.krona" ]]; then
				failure_flags="${failure_flags}-NO_krona_file(weighted)"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-kraken_weighted_classify_failed"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "weightedClassify" ]]; then
			if [[ -s "${processed}/${1}/${sample_name}/kraken/postAssembly/${1}_kraken_summary_assembled_BP.txt" ]]; then
				domain=$(sed -n '2p' "${OUTDATADIR}/kraken/postAssembly/${1}_kraken_summary_assembled_BP.txt" | cut -d' ' -f2)
				if (( $(echo "${domain} <= 0" | bc -l) )); then
					failure_flags="${failure_flags}-NO_weighted_classified_contigs(post)"
					failures=$(( failures + 1 ))
				else
					failure_flags="${failure_flags}-Top_match_under_50%_likely_contaminated"
					failures=$(( failures + 1 ))
				fi
			else
				failure_flags="${failure_flags}-NO_kraken_file(post)"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "weightedContam." ]]; then
			failure_flags="${failure_flags}-Multiple_species_found_above_contamination_threshold(${contamination_threshold})_kraken(weighted)"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "QUAST" ]]; then
			failure_flags="${failure_flags}-NO_QUAST_report"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "QUAST_plasFlow" ]]; then
			failure_flags="${failure_flags}-NO_QUAST_PlasFlow_report"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "Taxa" ]]; then
			warning_flags="${warning_flags}-NO_taxonomy_was_determined"
			warnings=$(( warnings + 1 ))
		elif [[ "${tool}" == "Assembly_Ratio" ]]; then
			if [[ "${tool_details}" == *"Too large"* ]]; then
				failure_flags="${failure_flags}-Assembly_is_>1.2x"
				failures=$(( failures + 1 ))
			elif [[ "${tool_details}" == *"Too small"* ]]; then
				failure_flags="${failure_flags}-Assembly_is_<.8x"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "Rawcoverage" ]]; then
			failure_flags="${failure_flags}-Raw_coverage_below_40x"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "prokka" ]]; then
			failure_flags="${failure_flags}-NO_PROKKA_GBF_file"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "BUSCO" ]]; then
			#echo "${processed}/${1}/${sample_name}/BUSCO/short_summary_${sample_name}.txt"
			if [[ -s "${processed}/${1}/${sample_name}/BUSCO/short_summary_${sample_name}.txt" ]]; then
				failure_flags="${failure_flags}-BUSCO_<90%"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-NO_BUSCO_summary"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "ANI" ]]; then
			#echo "${tool_details}"
			if [[ "${tool_details}" == *"% identity"* ]]; then
				failure_flags="${failure_flags}-ANI_identity<95%"
				failures=$(( failures + 1 ))
			elif [[ "${tool_details}" == *"% coverage"* ]]; then
				failure_flags="${failure_flags}-ANI_coverage<${ani_coverage_threshold}%"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-NO_ANI_output"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "ANI_REFSEQ" ]]; then
			#echo "${tool_details}"
			if [[ "${tool_details}" == *"% identity"* ]]; then
				failure_flags="${failure_flags}-ANI_identity<95%"
				failures=$(( failures + 1 ))
			elif [[ "${tool_details}" == *"% coverage"* ]]; then
				failure_flags="${failure_flags}-ANI_coverage<${ani_coverage_threshold}%"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-NO_ANI_output"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "ANI_OSII" ]]; then
			#echo "${tool_details}"
			if [[ "${tool_details}" == *"% identity"* ]]; then
				failure_flags="${failure_flags}-ANI_identity<95%"
				failures=$(( failures + 1 ))
			elif [[ "${tool_details}" == *"% coverage"* ]]; then
				failure_flags="${failure_flags}-ANI_coverage<${ani_coverage_threshold}%"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-NO_ANI_output"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "c-SSTAR" ]]; then
			failure_flags="${failure_flags}-NO_c-SSTAR_summary"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "c-SSTAR_plasFlow" ]]; then
			failure_flags="${failure_flags}-NO_c-sstar_output_on_plasFlow_assembly"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "GAMA" ]]; then
			failure_flags="${failure_flags}-NO_GAMA_summary"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "GAMA_plasFlow" ]]; then
			failure_flags="${failure_flags}-NO_GAMA_output_on_plasFlow_assembly"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "srst2" ]]; then
			failure_flags="${failure_flags}-NO_srst2_output"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "MLST" ]]; then
			if [[ -s "${processed}/${1}/${sample_name}/MLST/${1}_Pasteur.mlst" ]]; then
				failure_flags="${failure_flags}-NO_SCHEME-Check_that_pubMLST_has_${genus^}_${species})"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-NO_MLST_output"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "MLST-srst2" ]]; then
				failure_flags="${failure_flags}-ST_indeterminable"
				failures=$(( failures + 1 ))
		elif [[ "${tool}" == "16s_best_hit" ]]; then
			if [[ -s "${processed}/${1}/${sample_name}/16s/${sample_name}_16s_blast_id.txt" ]]; then
				failure_flags="${failure_flags}-No_taxonomy_assigned_in_16s_best"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-NO_16s_blast_id.txt file"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "16s_largest_hit" ]]; then
			if [[ -s "${processed}/${1}/${sample_name}/16s/${sample_name}_16s_blast_id.txt" ]]; then
				failure_flags="${failure_flags}-No_taxonomy_assigned_in_16s_largest"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-NO_16s_blast_id.txt file"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "16s" ]]; then
			failure_flags="${failure_flags}-NO_16s"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "Raw_read_count" ]]; then
			failure_flags="${failure_flags}-NO_Raw_Reads"
			failures=$(( failures + 1 ))
		elif [[ "${tool}" == "plasmid" ]]; then
			if [[ -d "${processed}/${1}/${sample_name}/plasmidFinder/" ]]; then
				failure_flags="${failure_flags}-NO_plasmidFinder_summary_file"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-NO_plasmidFinder_directory_for_full_assembly"
				failures=$(( failures + 1 ))
			fi
		elif [[ "${tool}" == "plasmid-plasmidAsmb" ]]; then
			if [[ -d "${processed}/${1}/${sample_name}/plasmidFinder_on_plasFlow/" ]]; then
				failure_flags="${failure_flags}-NO_plasmidFinder_summary_file"
				failures=$(( failures + 1 ))
			else
				failure_flags="${failure_flags}-NO_plasmidFinder_directory_for_plasmid_assembly"
				failures=$(( failures + 1 ))
			fi
		fi
	fi
done < "${sum_file}"
total_samples=$(( success_samples + failed_samples + warning_samples ))
printf "\n%-20s: %-8s\\n" "Total samples" "${total_samples}"
printf "%-20s: %-8s\\n" "Successful samples" "${success_samples}"
printf "%-20s: %-8s\\n" "Warning samples" "${warning_samples}"
printf "%-20s: %-8s\\n\n" "Failed samples" "${failed_samples}"

#if [[ ${#imperfect_samples[@]} -gt 0 ]]; then
#	for imperfect in "${imperfect_samples[@]}"; do
#		echo -e "${imperfect}"
#	done
#	echo -e "\n\n"
#else
#	echo "All samples completed successfully, with nothing noteworthy"
#fi
