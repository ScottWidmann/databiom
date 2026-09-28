#!/bin/bash
#SBATCH -J cellranger_count
#SBATCH -o cellranger_count_aggr.out
#SBATCH -N 1 -n 32

module load cellranger

cellranger count --id=treatment_1 --fastqs=./SRR21832248 --transcriptome=$REFERENCE_HUMAN --create-bam=true

cellranger count --id=treatment_2 --fastqs=./SRR21832249 --transcriptome=$REFERENCE_HUMAN --create-bam=true

cellranger count --id=control_1 --fastqs=./SRR21832250 --transcriptome=$REFERENCE_HUMAN --create-bam=true

cellranger count --id=control_2 --fastqs=./SRR21832251 --transcriptome=$REFERENCE_HUMAN --create-bam=true

echo "sample_id,molecule_h5" > control_aggr.csv
echo "control_1,./control_1/outs/molecule_info.h5" >> control_aggr.csv
echo "control_2,./control_2/outs/molecule_info.h5" >> control_aggr.csv

cellranger aggr --id=control_aggr --csv=control_aggr.csv

echo "sample_id,molecule_h5" > treatment_aggr.csv
echo "treatment_1,./treatment_1/outs/molecule_info.h5" >> treatment_aggr.csv
echo "treatment_2,./treatment_2/outs/molecule_info.h5" >> treatment_aggr.csv

cellranger aggr --id=treatment_aggr --csv=treatment_aggr.csv

# Qualityy control with fastqc and multiqc
module load fastqc
module load multiqc

mkdir fastqc

for p in SRR*; do
        fastqc "$p"/*R1_001.fastq.gz "$p"/*R2_001.fastq.gz -t 32 -o ./fastqc
done

cd fastqc

multiqc .
