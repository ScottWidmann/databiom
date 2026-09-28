#!/bin/bash

#SBATCH -J GSE270286
#SBATCH -o GSE270286_sra.out
#SBATCH -N 1 -n 32

module load sratoolkit

while read p; do
	prefetch "$p"
	fasterq-dump -e 32 --outdir ./"$p"/ "$p"
	pigz -p 32 ./"$p"/*.fastq
	rm ./"$p"/*.sra
done < GSE270286.txt
