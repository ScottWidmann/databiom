#!/bin/bash

#SBATCH -J GSE184102
#SBATCH -o GSE184102_sra.out
#SBATCH -N 1 -n 32

module load sratoolkit

while read p; do
	prefetch "$p"
	fasterq-dump -e 32 --outdir ./"$p"/ "$p"
	pigz -p 32 ./"$p"/*.fastq
	rm ./"$p"/*.sra
done < SRR_Acc_List.txt
