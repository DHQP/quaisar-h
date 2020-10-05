#! /usr/bin/env bash
#!/bin/sh -l


#
# Description: Script to consolidate all configuration type settings for quasar pipeline and any tools contained within
# 	Just needs to be sourced within a script to acquire all variables stored within
#
# Usage: . ./config_template.sh
#
# Output location: No output created
#
# Modules required: None
#
# v1.0.2 (05/08/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

# Get hostname to help determine if certain tools can be run and how to specifically get others to run with the right options
hostname=$(hostname -f)
host=$(echo ${hostname} | cut -d'.' -f1)
#echo ${hostname}
if [[ "${host}" = "scicomp-mue-01" ]];
then
	host="Biolinux"
elif [[ "${host}" =~ ^("login01"|"aspen"|"login.aspen"|"login02"|"login2.aspen") ]];
then
	host="aspen_login"
elif [[ "${host:0:4}" = "node" ]];
then
	host="cluster:${host}"
elif [[ "${host}" = "scbs-mue-prod-01" ]];
then
	host="Biolinux2020"
else
	echo "Hostname (${host}) not recognized, exiting"
	exit 1
fi

############# General Options #############

#Location to store all Quaisar-H pipeline run logs
Quaisar_H_log_directory="/scicomp/groups/OID/NCEZID/DHQP/CEMB/QuAISAR_logs"
#shortcut to processed samples folder
processed="/scicomp/groups/OID/NCEZID/DHQP/CEMB/MiSeqAnalysisFiles"
# Locations of all scripts and necessary accessory files
shareScript="$(pwd)"
# Location to keep all temp files when doing mass qsubmissions
mass_qsub_folder="/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/mass_subs"
if [[ ! -d "/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/mass_subs" ]]; then
	mkdir -p "/scicomp/groups/OID/NCEZID/DHQP/CEMB/Nick_DIR/mass_subs"
fi
# Location of default Outbreak Analyses files
Phyl_OA="/scicomp/groups/OID/NCEZID/DHQP/CEMB/PhylogenyAnalysis"
# Local databases that are necessary for pipeline...ANI, BUSCO, star, adapters, phiX
local_DBs="/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases"
# Scicomp databases that are necessary for pipeline...eventually refseq, kraken, gottcha,
scicomp_DBs="/scicomp/reference"
# Maximum number of quaisar pipelines to be running concurrently
max_quaisars=25

#Instruments and locations of files stored by those instruments
miseq1="/scicomp/instruments/17-4-4248_Illumina-MiSeq-M04765"
miseq2="/scicomp/instruments/17-4-4248_Illumina-MiSeq-M01025"
miseq3="/scicomp/instruments/17-5-5248_Illumina-MiSeq-M02103"
miseq4="/scicomp/instruments/17-4-4123_Illumina-MiSeq-M03961"
pacbio="/scicomp/instruments/23-12-651_PacBio-RSII-RS42135"

# Create a list of instruments to check output of
all_instruments=($miseq1 $miseq2 $miseq3 $miseq4) # $pacbio)

# Number of processors requested by numerous applications within the pipeline
procs=12 # Number of processors

# Phred scoring scale to be used (33 or 64)
phred=33

############# Tool versions for recordings sake ###########
version_BBMap="38.26"
version_Trimmomatic="0.35"
version_Kraken="0.10.5"
version_GOTTCHA="1.0c"
version_SPAdes="3.13.0"
version_QUAST="4.3"
version_PROKKA="1.12"
version_BUSCO="3.0.2"
version_pyANI="0.2.7"
version_CSSTAR="1.1.01"
version_srst2="0.2.0"
version_GAMA="4.7.4"
version_MLST="2.16"
version_plasFlow="1.1"
version_BARRNAP="0.8"


############# Application Specific Options #############

#####BBDUK specific config options #####
#requested memory size block
bbduk_mem=Xmx20g
#Kmer length (1-31). Larger Kmer size results in greater specificity
bbduk_k=31
#hamming distance
bbduk_hdist=1
#location of phiX sequences
phiX_location="${local_DBs}/phiX.fasta"

#####Trimmomatic specific config options #####
#Tells trimmomatic to use single or paired ends
#trim_ends=SE
trim_endtype=PE

#Which scoring scale to use
trim_phred="phred$phred"

#Location of the Adapter FASTA file used for trimming
trim_adapter_location="${local_DBs}/adapters.fasta"

#Seeding mismatches
trim_seed_mismatch=2
#palindrome clip threshold
trim_palindrome_clip_threshold=30
#simple clip threshold
trim_simple_clip_threshold=10
#Minimum adapter length (in palindrome)
trim_min_adapt_length=8
#Keeps read if complete forward and reverse sequences are palindromes
trim_complete_palindrome=TRUE
#Window Size
trim_window_size=20
#Window quality
trim_window_qual=30
#Specifies minimum quality to keep a read
trim_leading=20
#Specifies minimum quality to keep a read
trim_trailing=20
#Specifies minimum length to keep a read
trim_min_length=50

##### SPAdes specific options #####
#Phred quality offset based on scoring used
spades_phred_offset=$phred
#Max memory in Gbs
spades_max_memory=32
#Coverage threshold (positive float, off or auto)
spades_cov_cutoff="auto"

##### ANI specific options #####
#Max number of samples to be kept (not including source sample) when creating the mash tree
max_ani_samples=20
ani_coverage_threshold=70
REFSEQ=$(find ${local_DBs}/ANI/REFSEQ_*.msh -maxdepth 1 -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
REFSEQ_date=$(echo ${REFSEQ} | rev | cut -d'/' -f1 | rev | cut -d'_' -f2 | cut -d'.' -f1,2)

##### c-SSTAR identity options #####
csstar_perfect=100
csstar_ultrahigh=99
csstar_high=98
csstar_medium=95
csstar_low=80
# Change to your liking
csstar_other=40

##### c-SSTAR standard settings #####
argannot_srst2=$(find ${local_DBs}/star/argannot_*_srst2.fasta -maxdepth 1 -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
#echo "ARG Summary found: ${argannot_srst2}"
resFinder_srst2=$(find ${local_DBs}/star/resFinder_*_srst2.fasta -maxdepth 1 -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
#echo "RES Summary found: ${resFinder_srst2}"
resGANNOT_srst2=$(find ${local_DBs}/star/ResGANNOT_*_srst2.fasta -maxdepth 1 -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
resGANNOT_previous_srst2=$(find ${local_DBs}/star/ResGANNOT_*_srst2.fasta -maxdepth 1 -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 2 | tail -n 1)
#echo "ResGANNOT Summary found: ${resGANNOT_srst2}"
ResGANNCBI_srst2=$(find ${local_DBs}/star/ResGANNCBI_*_srst2.fasta -maxdepth 1 -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
resGANNCBI_previous_srst2=$(find ${local_DBs}/star/ResGANNCBI_*_srst2.fasta -maxdepth 1 -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 2 | tail -n 1)
#echo "ResGANNOT Summary found: ${resGANNOT_srst2}"
argannot_srst2_filename=$(echo "${argannot_srst2}" | rev | cut -d'/' -f1 | rev | cut -d'_' -f1,2)
resFinder_srst2_filename=$(echo "${resFinder_srst2}" | rev | cut -d'/' -f1 | rev | cut -d'_' -f1,2)
resGANNOT_srst2_filename=$(echo "${resGANNOT_srst2}" | rev | cut -d'/' -f1 | rev | cut -d'_' -f1,2)
ResGANNCBI_srst2_filename=$(echo "${ResGANNCBI_srst2}" | rev | cut -d'/' -f1 | rev | cut -d'_' -f1,2)
# gapped (g) versus ungapped(u)
csstar_gapping="g"
# Identity % 100(p), 99(u), 98(h), 95(m), 80(low)
csstar_identity="h"
csstar_plasmid_identity="o"

##### kraken unclassified threshold ######
# Will throw a warning flag during run summary if percent of unclassified reads/contigs are above this value
unclass_flag=30
# MiniKraken DB (smaller, but effective option)
kraken_mini_db="${local_DBs}/minikrakenDB/"
#kraken_mini_db="/scicomp/agave/execution/database/public/references/organizations/CDC/NCEZID/kraken/cdc-20171227"
kraken_full_db="${scicomp_DBs}/kraken/1.0.0/kraken_db/"
# MiniKraken DB (smaller, but effective option)
kraken2_mini_db="${local_DBs}/minikraken2DB/"
#kraken_mini_db="/scicomp/agave/execution/database/public/references/organizations/CDC/NCEZID/kraken/cdc-20171227"
kraken2_full_db="${scicomp_DBs}/kraken/2.0.0/kraken_db/"
### MOVE THESE TO SHARE/DBS
# Kraken normal, specially made by Tom with bacteria,archae, and viruses
# kraken_db="/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases/kraken_BAV_17/"
# alternate one with bacteria, fungus, and viruses
# kraken_db="/scicomp/groups/OID/NCEZID/DHQP/CEMB/databases/kraken_BVF_16/"
contamination_threshold=25

##### gottcha #####
# gottcha DB
gottcha_db="${local_DBs}/gottcha/GOTTCHA_BACTERIA_c4937_k24_u30.species"

##### plasmidFinder ######
#percent identity to match
plasmidFinder_identity=95.00
#percent length minimum (although not found in command line version, yet)
plasmidFinder_length=60
#DB
#plasmidFinder_all_DB=$(find ${local_DBs}/plasmidFinder_DB/all_*.fsa -maxdepth 1 -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
#plasmidFinder_entero_DB=$(find ${local_DBs}/plasmidFinder_DB/enterobacteriaceae_*.fsa -maxdepth 1 -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)
#plasmidFinder_gpos_DB=$(find ${local_DBs}/plasmidFinder_DB/gram_positive_*.fsa -maxdepth 1 -type f -printf '%p\n' | sort -k2,2 -rt '_' -n | head -n 1)


########## Settings used by downstream scripts ##########

##### Project Parser (run_csstar_MLST_project_parser.sh) #####
# Cutoff values when consolidating csstar hits from sets of samples (projects)
# Minimum % length required to be included in report, otherwise gets placed in rejects file
project_parser_Percent_length=90
# Minimum % identity required to be included in report, otherwise gets placed in rejects file
project_parser_Percent_identity=98
# Minimum % length required to be included in report when looking at plasmid assembly, typically more leeway is given to plasmid only hits, otherwise gets placed in rejects file
project_parser_plasmid_Percent_identity=40

if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
