#!/bin/bash

# SLURM job options
#SBATCH --job-name=coverage
#SBATCH --account=pathogens
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00

echo "Hello, World!"

# Go to a temporary directory
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

# Load necessary modules	
module purge
module load BEDTools/2.30.0-GCC-11.2.0
module load SAMtools/1.17-GCC-11.2.0

# Get the BAM file for this task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $OUT_DIR/bam_list.txt)

# Extract the file name from the full path
file_name=$(basename ${BAM_FILE})

cp $BAM_FILE .
cp $BAM_FILE.bai .
cp $EXON_BED ./gene.bed

echo "Coverage for ${file_name}"
	
# Calculate coverage
bedtools coverage -a gene.bed -b ${file_name} > "${file_name}_coverage.txt"
	
# Output useful job stats
/usr/local/bin/showJobStats.scr 

# Get log file
LOG_FILE="$DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.out"

# Get memory usage
mem_max=$(awk '/Max memory available:/ {print $4}' "$LOG_FILE")
mem_used=$(awk '/Max memory used:/ {print $4}' "$LOG_FILE")

# Test for successful completion
if [[ ${mem_max} == ${mem_used} ]]; then

	# Array of numbers
	mem_array=(10.00 50.00 100.00 250.00 500.00)

	# Check if mem_max is the last number in the array
	if [[ "${mem_max}" == "${mem_array[-1]}" ]]; then
		echo "Failed with not enough memory!"
		echo "Coverage for ${file_name} failed with ${mem_max%.*}GB memory" >> $MAS_LOG
		cp $0 $OUT_DIR/Logs/$SLURM_JOB_NAME.sh
		cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.out $OUT_DIR/Logs/$SLURM_JOB_NAME
		mv $OUT_DIR/Logs/$SLURM_JOB_NAME/$SLURM_JOB_NAME.$SLURM_JOB_ID.out "$OUT_DIR/Logs/$SLURM_JOB_NAME/${file_name}.out"
		cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.err $OUT_DIR/Logs/$SLURM_JOB_NAME
		mv $OUT_DIR/Logs/$SLURM_JOB_NAME/$SLURM_JOB_NAME.$SLURM_JOB_ID.err "$OUT_DIR/Logs/$SLURM_JOB_NAME/${file_name}.err"
		exit 0  # Exit the script if mem_max is the last number

	fi

	# Find and set mem_next to the number after mem_max
	attempt=$(printf "%s\n" "${mem_array[@]}" | grep -n -F -x "${mem_max}" | cut -d: -f1)
	mem_next=${mem_array[$((attempt))]}
	mem_next="${mem_next%.*}GB"	
	echo "Failed with not enough memory! Starting again with ${mem_next} memory"	
	sbatch --export=ALL --array=${SLURM_ARRAY_TASK_ID} --mem-per-cpu=${mem_next} --output=$DEF_DIR/Logs/%x.%A/%x.%j.out --error=$DEF_DIR/Logs/%x.%A/%x.%j.err $0 | tee -a $MAS_LOG
	cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.out $OUT_DIR/Logs/$SLURM_JOB_NAME/Failed
	mv $OUT_DIR/Logs/$SLURM_JOB_NAME/Failed/$SLURM_JOB_NAME.$SLURM_JOB_ID.out "$OUT_DIR/Logs/$SLURM_JOB_NAME/Failed/${file_name}_attempt_${attempt}.out"
	cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.err $OUT_DIR/Logs/$SLURM_JOB_NAME/Failed
	mv $OUT_DIR/Logs/$SLURM_JOB_NAME/Failed/$SLURM_JOB_NAME.$SLURM_JOB_ID.err "$OUT_DIR/Logs/$SLURM_JOB_NAME/Failed/${file_name}_attempt_${attempt}.err"

else

	cp ${file_name}_coverage.txt $OUT_DIR/$SLURM_JOB_NAME
	echo "Completed!"
	echo "Coverage for ${file_name} completed with ${mem_max%.*}GB memory" >> $MAS_LOG
	cp $0 $OUT_DIR/Logs/$SLURM_JOB_NAME.sh
	cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.out $OUT_DIR/Logs/$SLURM_JOB_NAME
	mv $OUT_DIR/Logs/$SLURM_JOB_NAME/$SLURM_JOB_NAME.$SLURM_JOB_ID.out "$OUT_DIR/Logs/$SLURM_JOB_NAME/${file_name}.out"
	cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.err $OUT_DIR/Logs/$SLURM_JOB_NAME
	mv $OUT_DIR/Logs/$SLURM_JOB_NAME/$SLURM_JOB_NAME.$SLURM_JOB_ID.err "$OUT_DIR/Logs/$SLURM_JOB_NAME/${file_name}.err"

fi
