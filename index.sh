#!/bin/bash

# SLURM job options
#SBATCH --job-name=Extraction_pipeline
#SBATCH --mail-user=
#SBATCH --mail-type=ALL
#SBATCH --account=
#SBATCH --partition=
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB 
#SBATCH --time=168:00:00
#SBATCH --output=Logs/%x.%j.out
#SBATCH --error=Logs/%x.%j.err

echo "Hello, World!"

set -a  # Export all variables to child processes

# Function: Print a help message.
usage() {								 
	echo "Usage: $0 [ -R Reference Genome Location ] [ -B BAM files location ] [ -G BED file with genes ] [ -E BED file with exons ] [ -S List of samples ] [ -p padding ] [ -d minimum depth ] [ -c coverage threshold ] [ -i Job ID from previous run] [ -r Continue from step ]" 1>&2 
}

# Function: Exit with error.
exit_abnormal() {						 
	usage
	exit 1
}

# Get input options
OPTIND=1  # Reset option index
padding=10  # Default padding
mindepth=1  # Default minimum depth
coverage_threshold=0.9  # Default coverage threshold
resume_step="New_run"  # Default resume step
OUT_ID=$SLURM_JOB_ID  # Default output ID
i_flag=0  # Flag for custom ID

# Parse command-line options
while getopts ":R:B:G:E:S:p:d:c:i:r:" options; do	   
	case "${options}" in
		R)  # Reference Genome Location
			REF_GNM_DIR=${OPTARG}
			fna_files=("$REF_GNM_DIR"/*.fna)
			if [ ${#fna_files[@]} -eq 0 ]; then  
				echo "Error: -R $REF_GNM_DIR doesn't contain a .fna file"
				exit_abnormal
			elif [ ${#fna_files[@]} -gt 1 ]; then
				echo "Error: -R $REF_GNM_DIR contains more than one .fna file"
				exit_abnormal
			else
				REF_GNM="${fna_files[0]}"
			fi
			;;
		B)  # BAM files location
			BAM_DIR=${OPTARG}
			if [[ -z $(ls "$BAM_DIR"/*.bam 2> /dev/null) ]]; then  
				echo "Error: -B $BAM_DIR doesn't contain a .bam file"
				exit_abnormal
			elif [[ -z $(ls "$BAM_DIR"/*.bam.bai 2> /dev/null) ]]; then  
				echo "Error: -B $BAM_DIR doesn't contain a .bam.bai file"
				exit_abnormal
			else
				BAM_DIR=$(realpath $BAM_DIR)
			fi
			;;	
		G)  # BED file with genes
			GENE_BED=${OPTARG}
			if [ ! -f "$GENE_BED" ] ; then  
				echo "Error: -G $GENE_BED doesn't exist"
				exit_abnormal			
			else
				GENE_BED=$(realpath $GENE_BED)
				EXON_BED=$GENE_BED  # Default exon BED file
			fi			
			;;
		E)  # BED file with exons
			EXON_BED=${OPTARG}
			if [ ! -f "$EXON_BED" ]; then
				EXON_BED=$GENE_BED  # Use gene BED file if exon BED file doesn't exist
				echo "No valid exon file provided. Using $GENE_BED instead"			
			else
				EXON_BED=$(realpath $EXON_BED)
			fi
			;;
		S)  # List of samples
			SAMPLES=${OPTARG}
			num_cols=$(awk '{print NF}' $SAMPLES | sort -nu | tail -n 1)
			if (( num_cols < 2 || num_cols > 3 )); then
				echo "Error: -S ${SAMPLES} must contain two or three columns"
				exit_abnormal			
			else
				SAMPLES=$(realpath $SAMPLES)
			fi
			;;
		p)  # Padding
			padding=${OPTARG}
			if [[ $padding ]]; then
				echo "padding=${padding}"
			fi
			;;
		d)  # Minimum depth
			mindepth=${OPTARG}
			if [[ $mindepth ]]; then
				echo "mindepth=${mindepth}"
			fi
			;;
		c)  # Coverage threshold
			coverage_threshold=${OPTARG}
			if [[ $coverage_threshold ]]; then
				echo "coverage_threshold=${coverage_threshold}"
			fi
			;;
		i)  # Job ID from previous run
			OUT_ID=${OPTARG}
			if [[ ! -d "/group/pathogens/IAWS/Personal/WStephen/$SLURM_JOB_NAME.$OUT_ID" ]]; then
				echo "Error: -i is provided but $SLURM_JOB_NAME.$OUT_ID doesn't exist"
				exit_abnormal
			else
				i_flag=1  # Set custom ID flag
			fi
			;;
		r)  # Continue from step
			resume_step=${OPTARG}
			if [[ ${i_flag} == 0 ]]; then
				echo "Error: -r is provided but -i is not"
				exit_abnormal
			fi
			;;
		*)  # Invalid option
			exit_abnormal
			;;
	esac
done
shift $((OPTIND -1))  # Shift off the options and optional --.

pwd
DEF_DIR=$(pwd)

# Create output directory and logs directory
OUT_DIR=$DEF_DIR/$SLURM_JOB_NAME.$SLURM_JOB_ID
mkdir $OUT_DIR
mkdir $OUT_DIR/Logs
cd $OUT_DIR

# Save OUT_DIR with random variables
pwd
OUT_DIR=$(pwd)

# Set log file
MAS_LOG=$DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_JOB_ID.out

awk -v padding=${padding} 'BEGIN { FS="\t"; OFS=FS } {$2-=padding;$3+=padding}1' $GENE_BED | sed 's/-[0-9]\+/1/g' > $OUT_DIR/genes.bed

# Stop exporting all variables to child processes
set +a

#For cleaner logs
timestamp=$(date +%Y-%m-%d\ %H:%M:%S)
echo "${timestamp}"
MAS_LOG_DIR=$DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_JOB_ID
mkdir $MAS_LOG_DIR


# Go to a temporary directory
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

# Copy previous run to current run
if [ "$resume_step" != "New_run" ]; then
	cp -R $DEF_DIR/$SLURM_JOB_NAME.$OUT_ID/* $DEF_DIR/$SLURM_JOB_NAME.$SLURM_JOB_ID
	rm -r $OUT_DIR/Logs
	mkdir $OUT_DIR/Logs
	if [ "$resume_step" == "variant" ]; then
		rm -r $OUT_DIR/${resume_step}
		rm -r $OUT_DIR/extract
		rm -r $OUT_DIR/align
		rm -r $OUT_DIR/tree
	elif [ "$resume_step" == "extract" ]; then
		rm -r $OUT_DIR/${resume_step}
		rm -r $OUT_DIR/align
		rm -r $OUT_DIR/tree
	elif [ "$resume_step" == "align" ]; then
		rm -r $OUT_DIR/${resume_step}
		rm -r $OUT_DIR/tree
	elif [ "$resume_step" == "tree" ]; then
		rm -r $OUT_DIR/${resume_step}
	fi
fi

#----------------------------------------------------------------------#

# Count number of BAM files
cat $SAMPLES | awk '{ print $1}' | sed "s|^|$BAM_DIR\/|" | sed 's/$/.bam/' > $OUT_DIR/bam_list.txt
bam_count=$(wc -l < $OUT_DIR/bam_list.txt)

# If not resuming from Variant, Extract, Align or Tree, submit Coverage.sh
if [[ "$resume_step" != "variant" && "$resume_step" != "extract" && "$resume_step" != "align" && "$resume_step" != "tree" ]]; then
	echo "Coverage started at $(date +%Y-%m-%d\ %H:%M:%S)!"
	sbatch --export=ALL --array=1-${bam_count} --mem-per-cpu=10GB --output=$DEF_DIR/Logs/%x.%A/%x.%j.out --error=$DEF_DIR/Logs/%x.%A/%x.%j.err $DEF_DIR/coverage.sh

	# Test if Coverage.sh is done
	err_count=$(ls -1 "$OUT_DIR/Logs/coverage"/*.err 2>/dev/null | wc -l)  # Count the number of .err files
	last_count=0
	while [[ ${err_count} -ne ${bam_count} ]]; do # Wait for the number of .err files to match the number of genes
		sleep 60 # Wait for 60 seconds before checking again
		err_count=$(ls -1 "$OUT_DIR/Logs/coverage"/*.err 2>/dev/null | wc -l)
		if [[ ${err_count} -eq ${last_count} ]]; then
			((stall_count++))
		else
			stall_count=0
		fi
		last_count=${err_count}
		if [[ ${stall_count} -ge 60 ]]; then  # If no progress for 60 minutes
			# Check if the current job is the only one running
			if [[ $(squeue -h -t RUNNING | wc -l) -eq 1 ]]; then
				echo "Current job is the only one running. Something went wrong. Resubmitting failed jobs"
				# Obtain files that are finished
				ls $OUT_DIR/coverage > coverage_status.txt

				# Extract base names and sort them
				awk -F'/' '{print $NF}' $OUT_DIR/bam_list.txt | awk -F'.' '{print $1}' | sort > bam_names.txt
				awk -F'/' '{print $NF}' coverage_status.txt | awk -F'.' '{print $1}' | sort > coverage_status.txt

				# Find incomplete
				comm -23 bam_names.txt coverage_status.txt > incompleate_jobs.txt

				# Convert the line numbers to a comma-separated list
				failed_jobs=$(echo $(while read name; do grep -n "^$name$" bam_names.txt | cut -d: -f1; done < incompleate_jobs.txt) | tr ' ' ',')

				# Use the list in the sbatch command
				sbatch --export=ALL --array=${failed_jobs} --mem-per-cpu=10GB --output=$DEF_DIR/Logs/%x.%A/%x.%j.out --error=$DEF_DIR/Logs/%x.%A/%x.%j.err $DEF_DIR/coverage.sh

				stall_count=0  # Reset the stall count
			else
				echo "Jobs still running"
				stall_count=0  # Reset the stall count
			fi
		fi
	done
	sbatch $DEF_DIR/stats.sh -d "$OUT_DIR/Logs/coverage/"
	echo "Coverage completed!"
fi

#----------------------------------------------------------------------#

if [ ! -d "$OUT_DIR/bam_lists" ]; then

	echo "Coverage threshold of $coverage_threshold"

	# Get the list of txt files in the current directory
	files=("$OUT_DIR/coverage/"*.bam_coverage.txt)

	# Get the list of unique gene names from the files
	readarray -t genes < <(awk -F'[_\t]' '{print $4}' "${files[@]}" | sort -u)

	# Initialize the output file
	IFS=$'\t'; echo -e "bam\t${genes[*]}" > coverage_table.txt

	# Calculate the averages
	for file in "${files[@]}"; do
		echo -n "$(basename "${file}" | cut -d'_' -f1)" >> coverage_table.txt
		for gene in ${genes[@]}; do
			avg=$(awk -v gene="${gene}" '$4 ~ gene {sum += $8; count++} END {print (count ? sum / count : 0)}' "${file}")
			echo -ne "\t$avg" >> coverage_table.txt
		done
		echo "" >> coverage_table.txt
	done

	mkdir $OUT_DIR/bam_lists

	# Extract the columns, remove the first row, filter rows based on the threshold, and remove the second row
	for ((i=0; i<${#genes[@]}; i++)); do
		awk -v col=$((i+2)) -v threshold="${coverage_threshold}" 'NR != 1 && $col >= threshold {print $1}' coverage_table.txt > $OUT_DIR/bam_lists/${genes[i]}_bam_list.txt
	done

fi

#----------------------------------------------------------------------#

# Count number of genes
gene_count=$(wc -l < $OUT_DIR/genes.bed)

# If not resuming from Extract, Align or Tree, submit Variant.sh
if [[ "$resume_step" != "extract" && "$resume_step" != "align" && "$resume_step" != "tree" ]]; then
	echo "Variant started at $(date +%Y-%m-%d\ %H:%M:%S)!"
	sbatch --export=ALL --array=1-${gene_count} --output=$DEF_DIR/Logs/%x.%A/%x.%j.out --error=$DEF_DIR/Logs/%x.%A/%x.%j.err $DEF_DIR/variant.sh

	# Test if Variant.sh is done
	err_count=$(ls -1 "$OUT_DIR/Logs/variant"/*.err 2>/dev/null | wc -l)  # Count the number of .err files
	while [[ ${err_count} -ne ${gene_count} ]]; do  # Wait for the number of .err files to match the number of .bam files
		sleep 60  # Wait for 60 seconds before checking again
		err_count=$(ls -1 "$OUT_DIR/Logs/variant"/*.err 2>/dev/null | wc -l)
	done
	sbatch $DEF_DIR/stats.sh -d "$OUT_DIR/Logs/variant/"
	echo "Variant completed!"
fi

#----------------------------------------------------------------------#

# If not resuming from Align or Tree, submit Extract.sh
if [[ "$resume_step" != "align" && "$resume_step" != "tree" ]]; then
	echo "Extract started at $(date +%Y-%m-%d\ %H:%M:%S)!"
	sbatch --export=ALL --array=1-${gene_count} --output=$DEF_DIR/Logs/%x.%A/%x.%j.out --error=$DEF_DIR/Logs/%x.%A/%x.%j.err $DEF_DIR/extract.sh

	# Test if Extract.sh is done
	err_count=$(ls -1 "$OUT_DIR/Logs/extract"/*.err 2>/dev/null | wc -l)  # Count the number of .err files
	while [[ ${err_count} -ne ${gene_count} ]]; do  # Wait for the number of .err files to match the number of .bam files
		sleep 60  # Wait for 60 seconds before checking again
		err_count=$(ls -1 "$OUT_DIR/Logs/extract"/*.err 2>/dev/null | wc -l)
	done
	sbatch $DEF_DIR/stats.sh -d "$OUT_DIR/Logs/extract/"
	echo "Extract completed!"
fi

#----------------------------------------------------------------------#

# If not resuming from Tree, submit Align.sh
if [[ "$resume_step" != "tree" ]]; then
	echo "Align started at $(date +%Y-%m-%d\ %H:%M:%S)!"
	sbatch --export=ALL --array=1-${gene_count} --mem-per-cpu=512MB --output=$DEF_DIR/Logs/%x.%A/%x.%j.out --error=$DEF_DIR/Logs/%x.%A/%x.%j.err $DEF_DIR/align.sh

	# Test if Align.sh is done
	err_count=$(ls -1 "$OUT_DIR/Logs/align"/*.err 2>/dev/null | wc -l)  # Count the number of .err files
	while [[ ${err_count} -ne ${gene_count} ]]; do  # Wait for the number of .err files to match the number of genes
		sleep 60  # Wait for 60 seconds before checking again
		err_count=$(ls -1 "$OUT_DIR/Logs/align"/*.err 2>/dev/null | wc -l)
	done
	sbatch $DEF_DIR/stats.sh -d "$OUT_DIR/Logs/align/"
	echo "Align completed!"
fi

#----------------------------------------------------------------------#

# Submit Tree.sh
echo "Tree started at $(date +%Y-%m-%d\ %H:%M:%S)!"
sbatch --export=ALL --array=1-${gene_count} --mem-per-cpu=1GB --output=$DEF_DIR/Logs/%x.%A/%x.%j.out --error=$DEF_DIR/Logs/%x.%A/%x.%j.err $DEF_DIR/tree.sh

# Test if Align.sh is done
err_count=$(ls -1 "$OUT_DIR/Logs/tree"/*.err 2>/dev/null | wc -l)  # Count the number of .err files
while [[ ${err_count} -ne ${gene_count} ]]; do  # Wait for the number of .err files to match the number of genes
	sleep 60  # Wait for 60 seconds before checking again
	err_count=$(ls -1 "$OUT_DIR/Logs/tree"/*.err 2>/dev/null | wc -l)
done
sbatch $DEF_DIR/stats.sh -d "$OUT_DIR/Logs/tree/"
echo "Align completed!"

echo "Completed!"

# Output useful job stats
/usr/local/bin/showJobStats.scr 

# Copy scripts and logs to output directory
cp $0 $OUT_DIR/Logs/$SLURM_JOB_NAME.sh
cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_JOB_ID.out $OUT_DIR/Logs
mv $OUT_DIR/Logs/$SLURM_JOB_NAME.$SLURM_JOB_ID.out $OUT_DIR/Logs/$SLURM_JOB_NAME.out
cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_JOB_ID.err $OUT_DIR/Logs
mv $OUT_DIR/Logs/$SLURM_JOB_NAME.$SLURM_JOB_ID.err $OUT_DIR/Logs/$SLURM_JOB_NAME.err

find $DEF_DIR/Logs -type d -newermt "${timestamp}" -not -path "$MAS_LOG_DIR" -not -path "$DEF_DIR/Logs" > file_list.txt

while IFS= read -r file
do
    mv "$file" $MAS_LOG_DIR
done < file_list.txt
