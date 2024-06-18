#!/bin/bash

# SLURM job options
#SBATCH --job-name=extract         
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10GB 
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

# Load necessary modules
module purge
module load BCFtools/1.17-GCC-11.2.0
module load SAMtools/1.17-GCC-11.2.0
module load BEDTools/2.30.0-GCC-11.2.0

# Print coverage threshold
echo "Coverage threshold of $coverage_threshold"

# Extract variables from the bed file for SAMtools
cp $OUT_DIR/genes.bed .
region=$(echo $(sed -n "${SLURM_ARRAY_TASK_ID}p" genes.bed) | awk 'BEGIN { FS=" " } { print $1":"$2"-"$3}' )
gene=$(echo $(sed -n "${SLURM_ARRAY_TASK_ID}p" genes.bed) | awk '{print $4}')

cp $OUT_DIR/bam_lists/${gene}_bam_list.txt .
cat ${gene}_bam_list.txt | sed "s|^|$BAM_DIR\/|" > ${gene}_bam_cp_list.txt
xargs -a ${gene}_bam_cp_list.txt -I {} -P 0 cp -t . {} {}.bai
cp $OUT_DIR/variant/${gene}.vcf.gz .
cp $REF_GNM ./reference.fna
cp $REF_GNM.fai ./reference.fna.fai

# Extract the target region from the reference genome to align against
echo "$(sed -n "${SLURM_ARRAY_TASK_ID}p" genes.bed)" | awk 'BEGIN { FS=" "} {print $1 "\t" $2 "\t" $3}' > tmp.bed
bedtools getfasta -fi reference.fna -bed tmp.bed > ${gene}.fa
cp ${gene}.fa $OUT_DIR/$SLURM_JOB_NAME/${gene}.fa

# Loop through each sample id
for bam in $(bcftools query -l ${gene}.vcf.gz ); do {

	# Subset joint vcf file to just that individual and retain variant sites only 
	bcftools view ${gene}.vcf.gz -m2 -M2 -s ${bam} --min-ac 1 -Ob -o tmp_${bam}.bcf
			
	# Index the BCF file
	tabix tmp_${bam}.bcf
			
	# Subset bam to target regions to count coverage
	samtools view -b ${bam}.bam ${region} > subset_${bam}.bam
			
	# Find low coverage sites (those below mindepth) and output bedgraph format - these sites will be masked with N bases
	bedtools genomecov -bga -ibam subset_${bam}.bam | awk -v m=${mindepth} '$4 < m' > low_coverage_sites_${bam}.bed
			
	# Copy low coverage sites to mask file
	cat low_coverage_sites_${bam}.bed > mask_${bam}.bed
			
	# Call consensus for regions in the bed file - Masking low coverage sites
	samtools faidx reference.fna ${region} | bcftools consensus tmp_${bam}.bcf -s ${bam} -p ${bam}_ --mark-del '-' -m mask_${bam}.bed -H I > ${gene}_${bam}.fa
			
	# Copy the resulting fasta file to the output directory
	cp ${gene}_${bam}.fa $OUT_DIR/$SLURM_JOB_NAME

} & done

wait

echo "Extracted ${gene} from all bam files." >> $MAS_LOG

# Output useful job stats
/usr/local/bin/showJobStats.scr 

echo "Completed!"
cp $0 $OUT_DIR/Logs/$SLURM_JOB_NAME.sh
cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.out $OUT_DIR/Logs/$SLURM_JOB_NAME
mv $OUT_DIR/Logs/$SLURM_JOB_NAME/$SLURM_JOB_NAME.$SLURM_JOB_ID.out "$OUT_DIR/Logs/$SLURM_JOB_NAME/${gene}.out"
cp $DEF_DIR/Logs/$SLURM_JOB_NAME.$SLURM_ARRAY_JOB_ID/$SLURM_JOB_NAME.$SLURM_JOB_ID.err $OUT_DIR/Logs/$SLURM_JOB_NAME
mv $OUT_DIR/Logs/$SLURM_JOB_NAME/$SLURM_JOB_NAME.$SLURM_JOB_ID.err "$OUT_DIR/Logs/$SLURM_JOB_NAME/${gene}.err"
