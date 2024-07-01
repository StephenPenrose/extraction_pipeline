#!/bin/bash

# SLURM job options
#SBATCH --job-name=stats
#SBATCH --account=pathogens
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB 
#SBATCH --time=1:00:00
#SBATCH --output=/dev/null

# Parse command line options
while getopts d: option
do
case "${option}"
in
d) DIR=${OPTARG};;
esac
done

echo "Hello, World!"

# Initialize variables
mem_avail=1000000


min_mem_used=1000000
max_mem_used=0

min_mem_used_files=()
max_mem_used_files=()


min_time=1000000
max_time=0

min_time_files=()
max_time_files=()

# Initialize variables
declare -A mem_used_files
declare -A time_files

# Loop over all .out files in the specified directory
for file in $DIR/*.out; do
	# Extract memory available
	mem_avail=$(grep "Max memory available:" $file | awk '{print $4}')
	
    # Extract memory used
    mem_used_line=$(grep "Max memory used:" $file)

    # Check if memory is in MiB or GiB and convert to GiB if necessary
    if [[ $mem_used_line == *"MiB"* ]]; then
        mem_used=$(echo $mem_used_line | awk '{print $4/1024}')
    else
        mem_used=$(echo $mem_used_line | awk '{print $4}')
    fi

    # Update min and max memory used and their corresponding files
    if awk 'BEGIN {exit !('$min_mem_used'>'$mem_used')}'; then
        min_mem_used=$mem_used
		min_mem_dsp=$(echo $mem_used_line | awk '{print $4}')
		min_mem_unt=$(echo $mem_used_line | awk '{print $5}')
        min_mem_used_files=(`basename $file .out`)
    elif awk 'BEGIN {exit !('$min_mem_used'=='$mem_used')}'; then
        min_mem_used_files+=(`basename $file .out`)
    fi

    if awk 'BEGIN {exit !('$max_mem_used'<'$mem_used')}'; then
        max_mem_used=$mem_used
		max_mem_dsp=$(echo $mem_used_line | awk '{print $4}')
		max_mem_unt=$(echo $mem_used_line | awk '{print $5}')
        max_mem_used_files=(`basename $file .out`)
    elif awk 'BEGIN {exit !('$max_mem_used'=='$mem_used')}'; then
        max_mem_used_files+=(`basename $file .out`)
    fi

    # Store memory used and file name
    mem_used_files[`basename $file .out`]="$mem_used $(echo $mem_used_line | awk '{print $4}') $(echo $mem_used_line | awk '{print $5}')"
	
	# Extract elapsed time
    time=$(grep "Elapsed walltime:" $file | awk '{print $3}' | awk -F: '{ print ($1 * 3600) + ($2 * 60) + $3 }')

    # Update min and max elapsed time and their corresponding files
    if awk 'BEGIN {exit !('$min_time'>'$time')}'; then
        min_time=$time
		min_time_dsp=$(grep "Elapsed walltime:" $file | awk '{print $3}')
        min_time_files=(`basename $file .out`)
    elif awk 'BEGIN {exit !('$min_time'=='$time')}'; then
        min_time_files+=(`basename $file .out`)
    fi

    if awk 'BEGIN {exit !('$max_time'<'$time')}'; then
        max_time=$time
		max_time_dsp=$(grep "Elapsed walltime:" $file | awk '{print $3}')
        max_time_files=(`basename $file .out`)
    elif awk 'BEGIN {exit !('$max_time'=='$time')}'; then
        max_time_files+=(`basename $file .out`)
    fi

    # Store elapsed time and file name
    time_files[`basename $file .out`]="$time $(grep "Elapsed walltime:" $file | awk '{print $3}')"	
done

# Print results
echo "Max memory available: $mem_avail GiB" >> $DIR/stats.txt
echo "" >> $DIR/stats.txt
echo "Min memory used: $min_mem_dsp $min_mem_unt, Files: ${min_mem_used_files[@]}" >> $DIR/stats.txt
echo "Max memory used: $max_mem_dsp $max_mem_unt, Files: ${max_mem_used_files[@]}" >> $DIR/stats.txt
echo "" >> $DIR/stats.txt
echo "Min elapsed time: $min_time_dsp, Files: ${min_time_files[@]}" >> $DIR/stats.txt
echo "Max elapsed time: $max_time_dsp, Files: ${max_time_files[@]}" >> $DIR/stats.txt
echo "" >> $DIR/stats.txt

# Print memory
echo "Memory Usage:" >> $DIR/stats.txt
for file in $(for key in "${!mem_used_files[@]}"; do printf "%s\t%s\n" "$key" "${mem_used_files[$key]}"; done | sort -rnk2 | cut -f1); do
    echo -e "$file\t$(echo ${mem_used_files[$file]} | awk '{print $2 $3}')" >> $DIR/stats.txt
done
echo "" >> $DIR/stats.txt

# Print time
echo "Elapsed Time:" >> $DIR/stats.txt
for file in $(for key in "${!time_files[@]}"; do printf "%s\t%s\n" "$key" "${time_files[$key]}"; done | sort -rnk2 | cut -f1); do
    echo -e "$file\t$(echo ${time_files[$file]} | awk '{print $2}')" >> $DIR/stats.txt
done
