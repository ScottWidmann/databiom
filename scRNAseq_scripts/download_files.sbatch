#!/bin/bash

#SBATCH -J GSE214967
#SBATCH -o GSE214967_sra.out
#SBATCH -N 1 -n 1

module load sratoolkit

# Loop through the list and download files
while read p; do
	prefetch "$p"
	fastq-dump --split-files --outdir ./"$p"/ "$p"
	pigz -p 32 ./"$p"/*.fastq
	rm ./"$p"/*.sra
done < SRR_Acc_List.txt

# Change the names of the files to be compatible with Cell Ranger

mv SRR21832248/SRR21832248_3.fastq.gz SRR21832248/GBM_P024_VORT_S4_L001_R1_001.fastq.gz

mv SRR21832248/SRR21832248_4.fastq.gz SRR21832248/GBM_P024_VORT_S4_L001_R2_001.fastq.gz

mv SRR21832249/SRR21832249_3.fastq.gz SRR21832249/GBM_P024_VORT_S4_L002_R1_001.fastq.gz

mv SRR21832249/SRR21832249_4.fastq.gz SRR21832249/GBM_P024_VORT_S4_L002_R2_001.fastq.gz

mv SRR21832250/SRR21832250_3.fastq.gz SRR21832250/GBM_P024_DMSO_S3_L001_R1_001.fastq.gz

mv SRR21832250/SRR21832250_4.fastq.gz SRR21832250/GBM_P024_DMSO_S3_L001_R2_001.fastq.gz

mv SRR21832251/SRR21832251_3.fastq.gz SRR21832251/GBM_P024_DMSO_S3_L002_R1_001.fastq.gz

mv SRR21832251/SRR21832251_4.fastq.gz SRR21832251/GBM_P024_DMSO_S3_L002_R2_001.fastq.gz
