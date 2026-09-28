#!/bin/bash
#SBATCH -J cellranger_count
#SBATCH -o cellranger_count_test.out
#SBATCH -N 1 -n 32

module load cellranger

cellranger count --id=treatment_1 --fastqs=./SRR21832248 --transcriptome=$REFERENCE_HUMAN --create-bam=true

cellranger count --id=treatment_2 --fastqs=./SRR21832249 --transcriptome=$REFERENCE_HUMAN --create-bam=true

cellranger count --id=control_1 --fastqs=./SRR21832250 --transcriptome=$REFERENCE_HUMAN --create-bam=true

cellranger count --id=control_2 --fastqs=./SRR21832251 --transcriptome=$REFERENCE_HUMAN --create-bam=true
