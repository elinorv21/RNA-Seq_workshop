#!/bin/bash
#SBATCH --job-name=mouse_rmats								
#SBATCH --output=mouse_rmats_%j.out							
#SBATCH --error=mouse_rmats_%j.err							
#SBATCH --time=02:00:00									
#SBATCH --nodes=1										
#SBATCH --cpus-per-task=12								
#SBATCH --mem=32G									
#SBATCH --partition=savio2_htc  

# Define input and output paths: use absolute path not relative path
GROUP1_BAM_LIST="/global/scratch/users/your-username/mouse_data/star_results/renamed/base_bams.txt" 
GROUP2_BAM_LIST="/global/scratch/users/your-username/mouse_data/star_results/renamed/drug_bams.txt"  
ANNOTATION_GTF="/global/scratch/users/your-username/genome/annotation.gtf" 
LIB_TYPE="fr-unstranded"             
READ_LENGTH=150                    
NUM_THREADS=12                       
OUTPUT_DIR="/global/scratch/users/your-username/mouse_data/rmats_results" 
TMP_DIR="/global/scratch/users/your-username/mouse_data/tmp_rmats"

# Make sure output directories exist									
mkdir -p "$OUTPUT_DIR"										
mkdir -p "$TMP_DIR"

# Run rMATS
PYTHONPATH=. python ~/miniconda3/envs/rmats_env/bin/rmats.py \
  --b1 "$GROUP1_BAM_LIST" \
  --b2 "$GROUP2_BAM_LIST" \
  --gtf "$ANNOTATION_GTF" \
  --readLength "$READ_LENGTH" \
  --od "$OUTPUT_DIR" \
  --libType "$LIB_TYPE" \
  --nthread "$NUM_THREADS" \
  --tmp "$TMP_DIR" \
  --task both

echo "rMATS computation complete. Results are in: $OUTPUT_DIR"

