#!/bin/bash

#SBATCH -J alignment
#SBATCH -o alignment.out
#SBATCH -N 1 -n 32

module load star
module load RSEM

mkdir reference/star_index
mkdir reference/rsem_index
mkdir results
mkdir results/final

gunzip reference/GRCh38.primary_assembly.genome.fa.gz
gunzip reference/gencode.v47.primary_assembly.basic.annotation.gtf.gz

STAR --runMode genomeGenerate --runThreadN 32 --genomeDir reference/star_index --genomeFastaFiles reference/GRCh38.primary_assembly.genome.fa --sjdbGTFfile reference/gencode.v47.primary_assembly.basic.annotation.gtf

rsem-prepare-reference --gtf reference/gencode.v47.primary_assembly.basic.annotation.gtf reference/GRCh38.primary_assembly.genome.fa reference/rsem_index/rsem_index

for p in SRR*; do

        STAR --runThreadN 32 --genomeDir reference/star_index --readFilesIn "$p/${p}_1.fastq.gz" "$p/${p}_2.fastq.gz" --outFileNamePrefix results/"$p." --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --readFilesCommand gunzip -c

	rsem-calculate-expression -p 32 --paired-end --alignments results/"$p".Aligned.toTranscriptome.out.bam reference/rsem_index/rsem_index results/"$p"

	sed "s/\\.[0-9]*//" results/"$p".genes.results > results/final/"$p".genes.results.clean
done

pigz -p 32 results/*.bam
rm -r results/*.tmp
rm -r results/*.stat
rm -r results/*_STARtmp

# Quality control using fastqc and multiqc
module load fastqc

for p in SRR*; do
        fastqc "$p/${p}_1.fastq.gz" "$p/${p}_2.fastq.gz" -t 32 -o ./fastqc
done

