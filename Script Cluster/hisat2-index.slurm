#!/bin/bash

#SBATCH --mail-type=fail
#SBATCH --job-name="hisat2_index"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem=8G
#SBATCH --partition=pall

module add UHTS/Aligner/hisat/2.2.1
mkdir ./reference/hisat2
cd ./reference/hisat2
hisat2-build /data/users/jlan/rna_breast/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly.fa