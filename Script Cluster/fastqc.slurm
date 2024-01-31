#!/bin/bash

#SBATCH --job-name="Fastqc"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=3:00:00
#SBATCH --mem=8G
#SBATCH --partition=pall

source /data/users/jlan/rna_breast/module.sh

mkdir /data/users/jlan/rna_breast/qc_breast_cancer
cd /data/users/jlan/rna_breast/qc_breast_cancer

ln -s /data/users/jlan/rna_breast/reads/*R1*.fastq.gz .
ln -s /data/users/jlan/rna_breast/reads/*R2*.fastq.gz .

for k in `ls -1 *.fastq.gz`; 
do fastqc -t 2 ${k};
done

