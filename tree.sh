#!/bin/bash

# SLURM job options
#SBATCH --job-name=tree         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=102400KB 
#SBATCH --time=72:00:00
#SBATCH --account=pathogens

echo "Hello, World!"

# Go to a temporary directory
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

# Load necessary modules
module purge
module load IQ-TREE/2.2.2.6-gompi-2021b

# Get the gene name from the genes.bed file
cp $OUT_DIR/genes.bed .
gene=$(echo $(sed -n "${SLURM_ARRAY_TASK_ID}p" genes.bed) | awk '{print $4}')

# Find aligned gene file
cp $OUT_DIR/align/$(echo ${gene} "dp"${mindepth} | sed 's/ /_/g' | sed 's/:/_/g')_aligned.fa ./fasta.fa

# Infer maximum-likelihood tree
iqtree2 -s fasta.fa -alrt 1000 -B 1000 -T AUTO -pre ${gene}

cp -r *.treefile $OUT_DIR/$SLURM_JOB_NAME
mkdir $OUT_DIR/$SLURM_JOB_NAME/Other/${gene}
cp -r ./* $OUT_DIR/$SLURM_JOB_NAME/Other/${gene}
rm $OUT_DIR/$SLURM_JOB_NAME/Other/${gene}/fasta.fa
rm $OUT_DIR/$SLURM_JOB_NAME/Other/${gene}/genes.bed
rm $OUT_DIR/$SLURM_JOB_NAME/Other/${gene}/*.treefile

echo "Tree for ${gene} completed" >> $MAS_LOG
echo "Completed!"

# Output useful job stats
/usr/local/bin/showJobStats.scr 

cp $0 $OUT_DIR/Logs/$SLURM_JOB_NAME.sh
cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.out $OUT_DIR/Logs/$SLURM_JOB_NAME
mv $OUT_DIR/Logs/$SLURM_JOB_NAME/$SLURM_JOB_NAME.$SLURM_JOB_ID.out "$OUT_DIR/Logs/$SLURM_JOB_NAME/${gene}.out"
cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.err $OUT_DIR/Logs/$SLURM_JOB_NAME
mv $OUT_DIR/Logs/$SLURM_JOB_NAME/$SLURM_JOB_NAME.$SLURM_JOB_ID.err "$OUT_DIR/Logs/$SLURM_JOB_NAME/${gene}.err"
