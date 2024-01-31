# rnaseq_course_breast_cancer
**This is the analysis process and results of the RNA-Seq course！**

The samples are a subset from Eswaran et al. 2012, and the fastq files were downloaded through the Gene Expression Omnibus (GEO), accession GSE52194. **For each tumor type, the top three samples were selected.** The library preparation protocol did not preserve information on the transcribed strand (i.e., unstranded). The libraries were sequenced on an Illumina HiSeq 2000 in paired-end mode.

This study includes script files for processing the raw data for use on the Cluster. These files are located in the 'Script Cluster' directory.

Differential expression analysis was conducted using the R package DESeq2. Details of the process can be found in the DESeq2 folder.

Overrepresentation analysis was performed using the R package clusterProfiler. Details of the process can be found in the clusterProfiler folder.

The analysis with clusterProfiler relies on the results obtained from DESeq.
