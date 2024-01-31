#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="featurecounts"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=16:00:00
#SBATCH --mem=4G
#SBATCH --partition=pall

source /data/users/jlan/rna_breast/module.sh

bam_path="/data/users/jlan/rna_breast/hisat2_samtools_result/"
ref_path="/data/users/jlan/rna_breast/reference/"
result_path="/data/users/jlan/rna_breast/featurecounts"

mkdir -p ${result_path}


featureCounts -p -t exon -g gene_id -a ${ref_path}Homo_sapiens.GRCh38.110.gtf -o ${result_path}/featurecounts.txt ${bam_path}*_results/*sorted.bam 2> ${result_path}/error.log
