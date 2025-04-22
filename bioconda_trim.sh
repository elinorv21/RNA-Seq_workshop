#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --partition=savio2_htc
#SBATCH --output=trimmomatic_%j.out
#SBATCH --error=trimmomatic_%j.err
#SBATCH --cpus-per-task=8  # Adjust as needed
#SBATCH --mem=16GB         # Adjust as needed
#SBATCH --time=1:00:00     # Adjust as needed

# Define variables (replace with your files) 
INPUT_R1=reads_1.fastq.gz
INPUT_R2=reads_2.fastq.gz
OUTPUT_R1_PAIRED=reads_1_paired.fastq.gz
OUTPUT_R1_UNPAIRED=reads_1_unpaired.fastq.gz
OUTPUT_R2_PAIRED=reads_2_paired.fastq.gz
OUTPUT_R2_UNPAIRED=reads_2_unpaired.fastq.gz
ADAPTERS=adapters.fa  # Path to your adapter sequences

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
