#!/bin/bash
#SBATCH -J fastqc
#SBATCH -o fastqc.out
#SBATCH -N 1 -n 32

module load fastqc
module load multiqc

mkdir fastqc

for p in SRR*; do
        fastqc "$p"/*R1_001.fastq.gz "$p"/*R2_001.fastq.gz -t 32 -o ./fastqc
done

multiqc ./fastqc/
