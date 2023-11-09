export PATH=/software/bin:/usr/local/bin:$PATH;
export LANG=en_US.UTF-8
module use /software/module/
module add UHTS/Quality_control/fastqc/0.11.9
module add UHTS/Analysis/MultiQC/1.8
module add UHTS/Aligner/bwa/0.7.17
module add UHTS/Analysis/samtools/1.10
module add Emboss/EMBOSS/6.6.0
module add UHTS/Analysis/SOAPdenovo2/2.04.241;
module add UHTS/Assembler/abyss/2.0.2;
module add UHTS/Assembler/SPAdes/3.15.4;
module add UHTS/Analysis/vcftools/0.1.15;
module add UHTS/Analysis/HTSeq/0.9.1;
module add UHTS/Analysis/mummer/4.0.0beta1
module add UHTS/Aligner/bowtie2/2.3.4.1
module add UHTS/Analysis/prokka/1.13
module add UHTS/Analysis/aragorn/1.2.38
module add SequenceAnalysis/MultipleSequenceAlignment/mafft/7.475
module add SequenceAnalysis/MultipleSequenceAlignment/prank/170427
module add SequenceAnalysis/StructurePrediction/mcl/14.137
module add UHTS/Analysis/cd-hit/4.6.8
module add UHTS/Analysis/BEDTools/2.29.2
module add Phylogeny/FastTree/2.1.10
module add UHTS/Analysis/GenomeAnalysisTK/4.2.0.0
module add UHTS/Analysis/picard-tools/2.21.8
module add UHTS/Analysis/trimmomatic/0.36
module add UHTS/Analysis/Roary/3.11.0
module add UHTS/Quality_control/fastp/0.19.5
#module add UHTS/Assembler/canu/2.1
module add UHTS/Analysis/minimap2/2.17
#module add UHTS/Analysis/busco/4.1.4
module load Blast/ncbi-blast/2.9.0+ #keep this for prokka...
#purge_haplotig use conda 
#ipa use conda
#hifiasm install in scripts
#Flye install in scripts
export EMBOSS_FILTER=1
module load Development/java/17.0.6
