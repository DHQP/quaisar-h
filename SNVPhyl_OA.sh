#!/bin/bash -l

#$ -o snv_OA.out
#$ -e snv_OA.err
#$ -N snv_OA
#$ -cwd
#$ -q short.q

#Import the config file with shortcuts and settings
if [[ ! -f "./config.sh" ]]; then
	cp ./config_template.sh ./config.sh
fi
. ./config.sh

#
# Description: Script that sorts a list of samples into similar MLST groups and runs snvPhyl on each group
#
# Usage ./SNVPhyl_OA.sh   list_file   Outbreak_Name Output_folder
#
# Output location: Parameter
#
# Modules required: None
#
# v1.0.1 (1/15/2019)
#
# Created by Rich Stanton (njr5@cdc.gov)
#


# Checks for proper argumentation
if [[ $# -eq 0 ]]; then
	echo "No argument supplied to $0, exiting"
	exit 1
# Shows a brief uasge/help section if -h option used as first argument
elif [[ "$1" = "-h" ]]; then
	echo "Usage is ./snvphyl_OA.sh  list_file  Outbreak_Name Output_folder"
	echo "Will group similar MLST types from the list together and then run snvphyl on each group"
	exit 0
elif [[ ! "${1}" ]]; then
	echo "List ${1} does not exist"
	echo "EXITING..."
	exit 1
else
	echo "Cleaning ${processed}/${2}/${1}"
fi

List=$1
Name=$2
Folder=$3

if [[ ! -d ${Folder} ]]; then
	mkdir ${Folder}
fi

ml Python3/3.5.2

python MLST_compare.py -i $List -o $Folder/$Name

for k in $Folder/*.samples
do
	sample=$(basename $k)
	qsub qSNVPhyl.sh $k $Folder ${sample:0: -8}
done

cp $List $Folder/$List.original

ml -Python3/3.5.2
