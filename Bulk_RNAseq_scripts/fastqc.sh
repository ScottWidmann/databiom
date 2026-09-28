#!/bin/bash

#SBATCH -J fastqc
#SBATCH -o fastqc.out
#SBATCH -N 1 -n 32

module load fastqc

for p in SRR*; do
	fastqc "$p/${p}_1.fastq.gz" "$p/${p}_2.fastq.gz" -t 32 -o ./fastqc
done
