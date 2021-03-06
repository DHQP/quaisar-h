#!/bin/sh -l

#$ -o SNVPhyl_X.out
#$ -e SNVPhyl_X.err
#$ -N SNVPhyl_X
#$ -cwd
#$ -q short.q

#
# Description: Runs SNVPhyl on a group of samples to determine relatedness
#
# Usage: ./run_SNVPhyl.sh -l path_to_list_file (First sample on list will be reference) -o output_directory -n analysis_identifier [-c path_to_config]
#
# Output location: parameter
#
# Modules required: snvphyl-galaxy-cli/1.3.0, Python/2.7.13 Mash/2.0
#
# v1.0.2 (09/04/2020)
#
# Created by Nick Vlachos (nvx4@cdc.gov)
#

ml purge
ml snvphyl-galaxy-cli/1.3.0 -Python/2.7.15 Python2/2.7.13 Mash/2.0 Python3/3.5.2

#  Function to print out help blurb
show_help () {
	echo "Usage is ./run_SNVPhyl.sh -l path_to_list_file (First sample on list will be reference) -o output_directory -n analysis_identifier [-c path_to_config]"
	echo "Output is saved to output_directory/analysis_identifier"
}

options_found=0
while getopts ":h?l:n:o:c:" option; do
	options_found=$(( options_found + 1 ))
	case "${option}" in
		\?)
			echo "Invalid option found: ${OPTARG}"
      show_help
      exit 0
      ;;
		l)
			echo "Option -l triggered, argument = ${OPTARG}"
			list=${OPTARG};;
		o)
			echo "Option -o triggered, argument = ${OPTARG}"
			outdir=${OPTARG};;
		n)
			echo "Option -n triggered, argument = ${OPTARG}"
			run_name=${OPTARG};;
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


# Checks for proper argumentation
if [[ -z "${list}" ]] || [[ ! -f ${list} ]] ; then
	echo "Empty group name (${list}) or non-existent sample list file supplied to run_SNVPhyl.sh, exiting"
	exit 1
elif [[ -z "${outdir}" ]]; then
	echo "Empty output directory name, exiting"
	exit 1
elif [[ -z "${run_name}" ]]; then
	echo "Empty analysis identifier, exiting"
	exit 1
fi

# Not being run on cluster=no run
if [[ ${host} != "cluster"* ]]; then
	echo "No scheduling system, can not run SNVPhyl"
	exit 1
fi

# Sets output folder to group_name under Phylogeny_analyses in MMB_Data folder
OUTDATADIR=${outdir}/${run_name}
if [[ ! -d "${OUTDATADIR}/FASTQs" ]]; then
	mkdir -p "${OUTDATADIR}/FASTQs"
fi
if [[ -d "${OUTDATADIR}/output" ]]; then
	rm -r "${OUTDATADIR}/output"
fi

echo $(python3 -V)

${shareScript}/clean_list.sh -l ${list}
cp ${list} ${OUTDATADIR}
centroid_filename=$(basename ${list}).centroid
python3 ${shareScript}/Mash_centroid.py -i ${list} -o ${OUTDATADIR}/${centroid_filename}

ml -Python3/3.5.2 Python2/2.7.13

counter=0
while IFS= read -r var || [ -n "$var" ]; do
	echo "var:$var"
	sample_name=$(echo "${var}" | awk -F"/" '{print $2}' | tr -d '[:space:]')
	# SNVPhyl can simulate reads on assemblies, :asm at the end of the filename is the designation for this action, It is unused in SNVPhyl and just removed
	if [[ ${#sample} -gt 4 ]]; then
		if [[ ${sample: -4} = ":asm" ]]; then
			sample=${sample::-4}
		fi
	fi
	echo "sample:$sample"
	project=$(echo "${var}" | awk -F"/" '{print $1}' | tr -d '[:space:]')
	echo "project:$project"
	if [[ ${counter} -eq 0 ]]; then
		#echo "Setting reference as ${sample_name} from ${project}"
		ref=${sample_name}
		ref_proj=${project}
		counter=$(( counter + 1))
		continue
	fi
	echo "Copying: ${sample_name} from ${project}"
	# Copy over standard FASTQs not compressed
	if [[ -f "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq.gz" ]]; then
		cp "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq.gz" "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fq.gz"
	elif [[ -f "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.fq.gz" ]]; then
		cp "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.fq.gz" "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fq.gz"
	elif [[ -f "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq" ]]; then
		echo "Copying ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq"
		cp "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R1_001.paired.fq" "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fq"
		gzip "${OUTDATADIR}/FASTQs/${sample_name}_R1_001.fq"
	else
		echo "No zipped or unzipped trimmed R1 exists...."
		exit
	fi
	if [[ -f "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq.gz" ]]; then
		cp "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq.gz" "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fq.gz"
	elif [[ -f "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.fq.gz" ]]; then
		cp "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.fq.gz" "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fq.gz"
	elif [[ -f "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq" ]]; then
		echo "Copying ${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq"
		cp "${processed}/${project}/${sample_name}/trimmed/${sample_name}_R2_001.paired.fq" "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fq"
		gzip "${OUTDATADIR}/FASTQs/${sample_name}_R2_001.fq"
	else
		echo "No zipped or unzipped trimmed R2 exists...."
		exit
	fi
	counter=$((counter + 1))
done < ${OUTDATADIR}/${centroid_filename}

submitter=$(whoami)
echo "Reference is ${ref} from ${ref_proj}"
cp "${processed}/${ref_proj}/${ref}/Assembly/${ref}_scaffolds_trimmed.fasta" "${OUTDATADIR}/reference(${ref}-${submitter}).fasta"

owd=$(pwd)
cd ${OUTDATADIR}/

# Must change api-key everytime container gets restarted. This will only work on CDC network
snvphyl.py --galaxy-url http://snvphyl-galaxy.biotech.cdc.gov --galaxy-api-key a4e89f39f940dadf2cdea6cbdd753041 --fastq-dir ./FASTQs --reference-file "./reference(${ref}-${submitter}).fasta" --output-dir ./output --relative-snv-abundance 0.75 --min-coverage 10 --min-mean-mapping 30 --filter-density-threshold 2 --filter-density-window 11 --workflow-id "f2db41e1fa331b3e"

snv_all_est=$(tail -n 1 "${OUTDATADIR}/output/vcf2core.tsv")
snv_est=$(echo "${snv_all_est}" | cut -d '	' -f7)

sed -i "s/reference/${ref}/g" "${OUTDATADIR}/output/snvMatrix.tsv"
sed -i "s/reference/${ref}/g" "${OUTDATADIR}/output/phylogeneticTree.newick"

echo -e "\nReference:\t${ref}\nSNVPhyl core estimate:\t${snv_est}%\n" >> "${OUTDATADIR}/output/snvMatrix.tsv"

cp "${OUTDATADIR}/output/snvMatrix.tsv" "${OUTDATADIR}/${run_name}_snvMatrix.tsv"
cp "${OUTDATADIR}/output/phylogeneticTree.newick" "${OUTDATADIR}/${run_name}_SNVPhyl.newick"

ml -snvphyl-galaxy-cli/1.3.0 -Python2/2.7.13 -Mash/2.0

cd ${owd}

exit 0
