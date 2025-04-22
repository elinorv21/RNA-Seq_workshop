#!/bin/bash
#SBATCH --partition=savio2_htc
#SBATCH --output=trimmomatic_%j.out
#SBATCH --error=trimmomatic_%j.err
#SBATCH --job-name=trimmomatic_multi
#SBATCH --array=1-6           # Create an array of 6 tasks (for 6 samples)
#SBATCH --cpus-per-task=8     # Adjust as needed per sample
#SBATCH --mem=16GB            # Adjust as needed per sample
#SBATCH --time=2:00:00        # Adjust total time if needed

# Define the base names of your input files (assuming a consistent naming scheme)
BASE_NAMES=("sample1" "sample2" "sample3" "sample4" "sample5" "sample6")
ADAPTERS=adapters.fa  # Path to your adapter sequences

# Define the location of your toy data
toy_data="toy_data"
OUTPUT_DIR="${toy_data}/results_folder"

# Get the sample index from the SLURM_ARRAY_TASK_ID environment variable
SAMPLE_INDEX=$((${SLURM_ARRAY_TASK_ID} - 1))
SAMPLE=${BASE_NAMES[$SAMPLE_INDEX]}

# Define input and output file names based on the sample
INPUT_R1="${toy_data}/samples_directory/${SAMPLE}_R1.fastq.gz"
INPUT_R2="${toy_data}/samples_directory/${SAMPLE}_R2.fastq.gz"
OUTPUT_R1_PAIRED="${OUTPUT_DIR}/${SAMPLE}_R1_paired.fq.gz"
OUTPUT_R1_UNPAIRED="${OUTPUT_DIR}/${SAMPLE}_R1_unpaired.fq.gz"
OUTPUT_R2_PAIRED="${OUTPUT_DIR}/${SAMPLE}_R2_paired.fq.gz"
OUTPUT_R2_UNPAIRED="${OUTPUT_DIR}/${SAMPLE}_R2_unpaired.fq.gz"

# Activate your conda environment
# conda activate mouse_env

# Load necessary modules
module load java/22.0.1

# Run Trimmomatic for the current sample
echo "Processing sample: $SAMPLE"
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

echo "Finished processing sample: $SAMPLE"
