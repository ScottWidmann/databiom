#!/bin/bash
#SBATCH -J cellranger_count
#SBATCH -o cellranger_count.out
#SBATCH -N 1 -n 12

module load cellranger

cellranger count --id=run_count_treatment_1 --fastqs=./SRR21832248 --transcriptome=$REFERENCE_HUMAN --create-bam=true

#cellranger count --id=run_count_treatment_2 --fastqs=./SRR21832249 --sample=SRR21832249 --transcriptome=$REFERENCE_HUMAN --create-bam=true

#cellranger count --id=run_count_control_1 --fastqs=./SRR21832250 --sample=SRR21832250 --transcriptome=$REFERENCE_HUMAN --create-bam=true

#cellranger count --id=run_count_control_2 --fastqs=./SRR21832251 --sample=SRR21832251 --transcriptome=$REFERENCE_HUMAN --create-bam=true
