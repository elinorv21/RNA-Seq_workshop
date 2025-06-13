#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --time=01:00:00 #Adjust time as needed
#SBATCH --output=fastqc_%j.out
#SBATCH --partition=savio2_htc
#SBATCH --error=fastqc_%j.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8  # Adjust threads for FastQC
#SBATCH --mem=8G           # Adjust memory as needed

# Load necessary modules 
module load java/22.0.1
module load perl/5.38.0-gcc-11.4.0
module load bio/picard/3.0.0-gcc-11.4.0
module load bio/fastqc/0.12.1-gcc-11.4.0

# Processes all gzipped FASTQ files in the directory; replace username/toy_data/fastq
INPUT_FILES="/global/scratch/users/your-username/mouse_data/fastq/*.fastq.gz"  # Process all gzipped FASTQ files in the directory
# Or specify individual files:
# INPUT_FILES="sample1.fastq.gz sample2.fastq.gz"

# Output directory; replace username/toy_data/fastqc_results_slurm
OUTPUT_DIR=/global/scratch/users/your-username/mouse_data/fastqc_results_slurm

# Creates output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Run FastQC
fastqc -t $SLURM_CPUS_PER_TASK \
       -o $OUTPUT_DIR \
       $INPUT_FILES

echo "FastQC analysis complete."

