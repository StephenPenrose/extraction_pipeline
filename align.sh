#!/bin/bash

# SLURM job options
#SBATCH --job-name=align         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=20
#SBATCH --time=72:00:00
#SBATCH --account=pathogens

echo "Hello, World!"

# Go to a temporary directory
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

# Create output directory
mkdir $OUT_DIR/$SLURM_JOB_NAME
mkdir $OUT_DIR/Logs/$SLURM_JOB_NAME
mkdir $OUT_DIR/Logs/$SLURM_JOB_NAME/Failed

# Load necessary modules
module purge
module load BEDTools/2.30.0-GCC-11.2.0
module load MUSCLE/5.1-GCCcore-11.2.0

# Get the gene name from the genes.bed file
cp $OUT_DIR/genes.bed .
gene=$(echo $(sed -n "${SLURM_ARRAY_TASK_ID}p" genes.bed) | awk '{print $4}')

# Create output name for the region
outname=$(echo ${gene} "dp"${mindepth} | sed 's/ /_/g' | sed 's/:/_/g')

# Get a list of fasta files for that region
find $OUT_DIR/extract -type f -not -name "${gene}.fa" | grep ${gene} > ${gene}_fastas.txt

#Copy Refrence sequence
cp $OUT_DIR/extract/${gene}.fa ./${outname}.fa

# Get new names
ncol=$(awk -v 'FS=\t' '{print NF}' $SAMPLES | sort -nu | tail -n 1)
if [[ "$ncol" -eq 2 ]]; then
	cat $SAMPLES | awk '{ print $1, $2}' | sed 's/^.*\///g'| sed 's/.bam//g' > new_names.txt
elif [[ "$ncol" -eq 3 ]]; then
	cat $SAMPLES | awk '{ print $1, $3}' | sed 's/^.*\///g'| sed 's/.bam//g' > new_names.txt
fi

# Add fasta files together
while read f; do

    # Output record to multi-sample fasta
    current_name=$(head -n 1 $f | cut -d'_' -f1 | cut -d'>' -f2)
    new_name=$(grep -w "$current_name" new_names.txt | awk '{ print $2}' FS=' ')
    cat <(echo ">${new_name}") <(tail -n +2 $f) >> ${outname}.fa

done < ${gene}_fastas.txt


# Align records with muscle v5
muscle -align ${outname}.fa -output ${outname}_aligned.fa

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
	mem_array=(10.00 20.00 40.00 80.00 160.00 320.00 640.00)

	# Check if mem_max is the last number in the array
	if [[ "${mem_max}" == "${mem_array[-1]}" ]]; then
		echo "Failed with not enough memory!"
		echo "Align for ${gene} failed with ${mem_max%.*}GB memory" >> $MAS_LOG
		cp $0 $OUT_DIR/Logs/$SLURM_JOB_NAME.sh
		cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.out $OUT_DIR/Logs/$SLURM_JOB_NAME
		mv $OUT_DIR/Logs/$SLURM_JOB_NAME/$SLURM_JOB_NAME.$SLURM_JOB_ID.out "$OUT_DIR/Logs/$SLURM_JOB_NAME/${gene}.out"
		cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.err $OUT_DIR/Logs/$SLURM_JOB_NAME
		mv $OUT_DIR/Logs/$SLURM_JOB_NAME/$SLURM_JOB_NAME.$SLURM_JOB_ID.err "$OUT_DIR/Logs/$SLURM_JOB_NAME/${gene}.err"
		exit 0  # Exit the script if mem_max is the last number

	fi

	# Find and set mem_next to the number after mem_max
	attempt=$(printf "%s\n" "${mem_array[@]}" | grep -n -F -x "${mem_max}" | cut -d: -f1)
	mem_next=${mem_array[$((attempt))]}
	mem_next="$((${mem_next%.*} / 20))GB"	
	echo "Failed with not enough memory! Starting again with ${mem_array[$((attempt))]%.*}GB memory"		
	sbatch --export=ALL --array=${SLURM_ARRAY_TASK_ID} --mem-per-cpu=${mem_next} --output=$DEF_DIR/Logs/%x.%A/%x.%j.out --error=$DEF_DIR/Logs/%x.%A/%x.%j.err $0 | tee -a $MAS_LOG
	cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.out $OUT_DIR/Logs/$SLURM_JOB_NAME/Failed
	mv $OUT_DIR/Logs/$SLURM_JOB_NAME/Failed/$SLURM_JOB_NAME.$SLURM_JOB_ID.out "$OUT_DIR/Logs/$SLURM_JOB_NAME/Failed/${gene}_attempt_${attempt}.out"
	cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.err $OUT_DIR/Logs/$SLURM_JOB_NAME/Failed
	mv $OUT_DIR/Logs/$SLURM_JOB_NAME/Failed/$SLURM_JOB_NAME.$SLURM_JOB_ID.err "$OUT_DIR/Logs/$SLURM_JOB_NAME/Failed/${gene}_attempt_${attempt}.err"

else

	cp ${outname}_aligned.fa $OUT_DIR/$SLURM_JOB_NAME
	echo "Completed!"
	echo "Align for ${gene} completed with ${mem_max%.*}GB memory" >> $MAS_LOG
	cp $0 $OUT_DIR/Logs/$SLURM_JOB_NAME.sh
	cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.out $OUT_DIR/Logs/$SLURM_JOB_NAME
	mv $OUT_DIR/Logs/$SLURM_JOB_NAME/$SLURM_JOB_NAME.$SLURM_JOB_ID.out "$OUT_DIR/Logs/$SLURM_JOB_NAME/${gene}.out"
	cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.err $OUT_DIR/Logs/$SLURM_JOB_NAME
	mv $OUT_DIR/Logs/$SLURM_JOB_NAME/$SLURM_JOB_NAME.$SLURM_JOB_ID.err "$OUT_DIR/Logs/$SLURM_JOB_NAME/${gene}.err"

fi
