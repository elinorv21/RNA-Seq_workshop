#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --partition=savio2_htc
#SBATCH --output=trimmomatic_%j.out
#SBATCH --error=trimmomatic_%j.err
#SBATCH --cpus-per-task=8  # Adjust as needed
#SBATCH --mem=16GB         # Adjust as needed
#SBATCH --time=1:00:00     # Adjust as needed

# Define variables (replace with your files) 
INPUT_R1=toy_data/fastq/ctrl_rep1_forward.fastq.gz
INPUT_R2=toy_data/fastq/ctrl_rep1_reverse.fastq.gz
OUTPUT_R1_PAIRED=toy_data/trim_results/trim_paired_ctrl_rep1_forward.fastq.gz
OUTPUT_R1_UNPAIRED=toy_data/trim_results/trim_unpaired_ctrl_rep1_forward.fastq.gz
OUTPUT_R2_PAIRED=toy_data/trim_results/trim_paired_ctrl_rep1_reverse.fastq.gz
OUTPUT_R2_UNPAIRED=toy_data/trim_results/trim_unpaired_ctrl_rep1_reverse.fastq.gz
ADAPTERS=software/Trimmomatic-0.39/adapters/TruSeq3-PE.fa  # Path to your adapter sequences; this file comes from github download

# Activate your conda environment
# conda activate mouse_env
# Load necessary modules
module load java/22.0.1

# Run Trimmomatic
trimmomatic PE \
    -threads $SLURM_CPUS_PER_TASK \
    $INPUT_R1 $INPUT_R2 \
    $OUTPUT_R1_PAIRED $OUTPUT_R1_UNPAIRED \
    $OUTPUT_R2_PAIRED $OUTPUT_R2_UNPAIRED \
    ILLUMINACLIP:$ADAPTERS:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:15 \
    MINLEN:36
