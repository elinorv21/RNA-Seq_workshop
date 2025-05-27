#!/bin/bash
#SBATCH --job-name=fcounts
#SBATCH --output=fcounts_%j.out
#SBATCH --error=fcounts_%j.err
#SBATCH --time=01:00:00  
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4  
#SBATCH --mem=16G          # Adjust memory
#SBATCH --partition=savio2_htc  

# Activate conda environment
# source ~/.bashrc  # or source ~/.zshrc
# conda activate mouse_env

# Define variables 
BAM_FILES=(
/global/scratch/users/elinorvelasquez/mouse_data/star_results/renamed/base_rep1.bam
/global/scratch/users/elinorvelasquez/mouse_data/star_results/renamed/base_rep2.bam
/global/scratch/users/elinorvelasquez/mouse_data/star_results/renamed/base_rep3.bam
/global/scratch/users/elinorvelasquez/mouse_data/star_results/renamed/base_rep4.bam
/global/scratch/users/elinorvelasquez/mouse_data/star_results/renamed/base_rep5.bam
/global/scratch/users/elinorvelasquez/mouse_data/star_results/renamed/drug_rep1.bam
/global/scratch/users/elinorvelasquez/mouse_data/star_results/renamed/drug_rep2.bam
/global/scratch/users/elinorvelasquez/mouse_data/star_results/renamed/drug_rep3.bam
/global/scratch/users/elinorvelasquez/mouse_data/star_results/renamed/drug_rep4.bam
/global/scratch/users/elinorvelasquez/mouse_data/star_results/renamed/drug_rep5.bam
)
GTF_FILE=genome/annotation.gtf
OUTPUT_FILE=mouse_data/fcounts_results/gene_counts.txt

# Run featureCounts with explicit thread control
featureCounts -T $SLURM_CPUS_PER_TASK \
-a $GTF_FILE \
-o $OUTPUT_FILE \
-t exon -g gene_id \
-p -B -C \
--primary --countReadPairs \
"${BAM_FILES[@]}"

echo "featureCounts complete."
