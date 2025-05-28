# RNA-Seq_workshop
## SLURM files
### mouse_star.sh: 
This slurm file runs the STAR software on a single fastq file (Trimmomatic output) to align reads to a reference genome. Output is a bam file. Bai files need to be produced separately. The genome index (index_star) needs to be produced separately. The reference genome (genome.fa) and annotation file (annotation.gtf) need to be downloaded.
### mouse_rmats.sh:
This slurm file runs the rMATS software on bam files to produce text files representing a variety of alternative splicing events. The tmp file needs to be emptied for each run.
### mouse_fcounts.sh:
This slurm file runs the featureCounts software on bam files to produce a gene_counts text file for differential gene expression analysis.

### Trimmomatic files:
### bioconda_trim.sh: 
This slurm file describes a generic version of a slurm file for processing one paired-ended sample (two fastq files: read1 and read2). 
### bioconda_trim_example.sh: 
This slurm file is an example of a slurm file processing one paired-ended sample (two fastq files: forward and reverse). 
### trim_loop.sh: 
This slurm file describes a generic version of a slurm file for processing multiple paired-ended samples (each sample with two fastq files).
### trim_loop_example.sh: 
This slurm file is an example of a slurm file processing multiple paired-ended samples (each sample with two fastq files: forward and reverse).

## R files
### mouse_dge.R:
This R file produces a variety of plots helpful for analyzing differentially expressed genes, using the DESeq2 and edgeR R packages.
### mouse_cluster.R:
This R file produces a variety of plots helpful for performing functional analysis, using the clusterProfiler R package.
### mouse_rmats.R:
This R file produces a variety of plots helpful for analyzing alternative splicing events.
