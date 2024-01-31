#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="hisat2_mapping"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --time=18:00:00
#SBATCH --mem=50G
#SBATCH --partition=pall
#SBATCH --array=1-12

source /data/users/jlan/rna_breast/module.sh

reads_path="/data/users/jlan/rna_breast/reads/"
ref_path="/data/users/jlan/rna_breast/reference/hisat2/"
result_path="/data/users/jlan/rna_breast/hisat2_samtools_result"
error_logs_path="${result_path}/error_logs"

mkdir -p ${result_path}
mkdir -p ${error_logs_path}

samples=(
    HER21
    HER22
    HER23
    NonTNBC1
    NonTNBC2
    NonTNBC3
    Normal1
    Normal2
    Normal3
    TNBC1
    TNBC2
    TNBC3
)

NAME=${samples[$SLURM_ARRAY_TASK_ID-1]}

sample_dir="${result_path}/${NAME}_results"
mkdir -p ${sample_dir}

{
    hisat2 -x ${ref_path}Homo_sapiens.GRCh38.dna.primary_assembly.fa -1 ${reads_path}${NAME}_R1.fastq.gz -2 ${reads_path}${NAME}_R2.fastq.gz -S ${sample_dir}/${NAME}.sam -p 4
    samtools view -hbS ${sample_dir}/${NAME}.sam > ${sample_dir}/${NAME}.bam
    samtools sort -m 4G -@ 4 -o ${sample_dir}/${NAME}_sorted.bam ${sample_dir}/${NAME}.bam
    samtools index ${sample_dir}/${NAME}_sorted.bam
} 2> ${error_logs_path}/${NAME}_error.log
