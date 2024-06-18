#!/bin/bash

# SLURM job options
#SBATCH --job-name=variant         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=20GB 
#SBATCH --time=72:00:00
#SBATCH --account=pathogens

echo "Hello, World!"

# Go to a temporary directory
cd $TMPDIR
tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
cd ${tmp_dir}
pwd

# Create necessary directories
mkdir $OUT_DIR/$SLURM_JOB_NAME
mkdir $OUT_DIR/Logs/$SLURM_JOB_NAME

# Load necessary modules
module purge
module load BCFtools/1.17-GCC-11.2.0

# Extract variables from the bed file for SAMtools
cp $OUT_DIR/genes.bed .
region=$(echo $(sed -n "${SLURM_ARRAY_TASK_ID}p" genes.bed) | awk 'BEGIN { FS=" " } { print $1":"$2"-"$3}' )
gene=$(echo $(sed -n "${SLURM_ARRAY_TASK_ID}p" genes.bed) | awk '{print $4}')

cp $OUT_DIR/bam_lists/${gene}_bam_list.txt .
cat ${gene}_bam_list.txt | sed "s|^|$BAM_DIR\/|" > ${gene}_bam_cp_list.txt
xargs -a ${gene}_bam_cp_list.txt -I {} -P 0 cp -t . {} {}.bai
cp $REF_GNM ./reference.fna
cp $REF_GNM.fai ./reference.fna.fai
	
ncol=$(awk -v 'FS=\t' '{print NF}' $SAMPLES | sort -nu | tail -n 1)
if [[ "$ncol" -eq 3 ]]; then
	cat $SAMPLES | awk -v 'FS=\t' -v 'OFS=\t' '{ print $1, $2}' | sed 's/ /_/g' | sort | uniq > groups.txt
fi

# Check if groups are defined in the manifest file
if [ -f groups.txt ] ; then 

	# If true - Call variants, grouping samples by population/species
	bcftools mpileup -f reference.fna -b ${gene}_bam_list.txt -r ${region} -C 50 -E -q30 -Q20 -a FORMAT/AD,FORMAT/DP,FORMAT/QS |\
	bcftools call -mv -G groups.txt  |\
	bcftools norm -f reference.fna -Ou |\
	bcftools filter --IndelGap 5 -Ou |\
	bcftools sort - -Oz -o ${gene}.vcf.gz

else 

	# If false - Call variants considering all samples as the same population
	bcftools mpileup -f reference.fna -b ${gene}_bam_list.txt -r ${region} -C 50 -E -q30 -Q20 -a FORMAT/AD,FORMAT/DP,FORMAT/QS |\
	bcftools call -mv |\
	bcftools norm -f reference.fna -Ou |\
	bcftools filter --IndelGap 5 -Ou |\
	bcftools sort - -Oz -o ${gene}.vcf.gz

fi

# Index resulting VCF
tabix ${gene}.vcf.gz

# Copy VCF and index to output directory
cp ${gene}.vcf.gz $OUT_DIR/$SLURM_JOB_NAME
cp ${gene}.vcf.gz.tbi $OUT_DIR/$SLURM_JOB_NAME
		
# Decompress .vcf.gz and save to output
pigz -p ${SLURM_CPUS_PER_TASK} -d -c ${gene}.vcf.gz > $OUT_DIR/$SLURM_JOB_NAME/${gene}.vcf

echo "Variants for ${gene} called" >> $MAS_LOG

# Output useful job stats
/usr/local/bin/showJobStats.scr 

echo "Completed!"

cp $0 $OUT_DIR/Logs/$SLURM_JOB_NAME.sh
cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.out $OUT_DIR/Logs/$SLURM_JOB_NAME
mv $OUT_DIR/Logs/$SLURM_JOB_NAME/$SLURM_JOB_NAME.$SLURM_JOB_ID.out "$OUT_DIR/Logs/$SLURM_JOB_NAME/${gene}.out"
cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.err $OUT_DIR/Logs/$SLURM_JOB_NAME
mv $OUT_DIR/Logs/$SLURM_JOB_NAME/$SLURM_JOB_NAME.$SLURM_JOB_ID.err "$OUT_DIR/Logs/$SLURM_JOB_NAME/${gene}.err"
