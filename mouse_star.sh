#!/bin/bash
#SBATCH --job-name=mouse_star_sample1
#SBATCH --output=mouse_star_%j.out
#SBATCH --error=mouse_star_%j.err
#SBATCH --time=2:00:00  
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4 # Adjust threads
#SBATCH --mem=64G          # STAR can be memory-intensive
#SBATCH --partition=savio2_htc  # Change to appropriate partition

# Define variables; this code is for a single sample (represented by two files: forward and reverse) 
GENOME_DIR=genome/index_star  # Path to STAR index from scratch base directory
INPUT_R1=mouse_data/trim_results/trim_paired_sample1_1.fastq.gz # _1 = forward (input data is trimmomatic output)
INPUT_R2=mouse_data/trim_results/trim_paired_sample1_2.fastq.gz # _2 = reverse (note: we are using relative path not absolute path)
OUTPUT_DIR=mouse_data/star_results
SAMPLE_NAME=sample1  # For output file naming

# Run STAR
STAR \
    --runMode alignReads \
    --twopassMode Basic \
    --genomeDir $GENOME_DIR \
    --readFilesIn $INPUT_R1 $INPUT_R2 \
    --readFilesCommand 'gunzip -c' \
    --outFileNamePrefix $OUTPUT_DIR/$SAMPLE_NAME. \
    --outSAMtype BAM SortedByCoordinate \
    --runThreadN 4

echo "STAR alignment complete."
