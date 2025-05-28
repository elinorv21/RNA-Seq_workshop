# RNA-Seq_workshop
## SLURM files
### mouse_star.sh: 
This slurm file runs the STAR algorithm on a single file to align reads to a reference genome.
### mouse_rmats.sh:
This slurm file runs the rMATS algorithm to produce text files representing a variety of alternative splicing events.

## Trimmomatic
### bioconda_trim.sh: 
This SLURM file describes a generic version of a slurm file for processing one paired-ended sample (two files: read1 and read2). 
### bioconda_trim_example.sh: 
This SLURM file is an example of a slurm file processing one paired-ended sample (two files: forward and reverse). 
### trim_loop.sh: 
This SLURM file describes a generic version of a slurm file for processing multiple paired-ended samples (each sample with two files).
### trim_loop_example.sh: 
This SLURM file is an example of a slurm file processing multiple paired-ended samples (each sample with two files: forward and reverse).

## R files
### mouse_dge.R:
This R file 
### mouse_cluster.R:
This R file
### mouse_rmats.R:
This R file
