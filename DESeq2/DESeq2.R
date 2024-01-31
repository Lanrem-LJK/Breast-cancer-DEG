library("DESeq2")
library("pheatmap")
library("RColorBrewer")
library("ggplot2")

#Input featurecounts and coldata
DESeq_featurecounts <- read.table("C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\DESeq_featurecounts.txt", header = TRUE, row.names = 1)
coldata <- transform(read.table("C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\coldata.txt", header = TRUE, row.names = 1), type = as.factor(type))
coldata$type <- as.factor(coldata$type)
coldata$type <- relevel(coldata$type, ref = "Normal")    #Changing the Default Reference Level to "Normal"

#Get Matrix
dds <- DESeqDataSetFromMatrix(countData = DESeq_featurecounts, colData = coldata, design = ~ type)
dds <- DESeq(dds)
resultsNames(dds)
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\counts.csv")

#Remove the dependence of the variance on the mean
ntd <- normTransform(dds)
vsd <- vst(dds, blind=TRUE)
rld <- rlog(dds, blind=TRUE)

#Get heatmap
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[, "type", drop = FALSE])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)

#Principal component plot
plotPCA(vsd, intgroup="type")

#Differential expression analysis [The results display comparisons of the earlier relative to the later.]
# 1.HER2 vs Normal
res_HER2_vs_Normal <- results(dds, contrast=c("type", "HER2", "Normal"))
res_HER2_Normal_Ordered <- res_HER2_vs_Normal[order(res_HER2_vs_Normal$pvalue),]

# 2.NonTNBC vs Normal
res_NonTNBC_vs_Normal <- results(dds, contrast=c("type", "NonTNBC", "Normal"))
res_NonTNBC_Normal_Ordered <- res_NonTNBC_vs_Normal[order(res_NonTNBC_vs_Normal$pvalue),]

# 3.TNBC vs Normal
res_TNBC_vs_Normal <- results(dds, contrast=c("type", "TNBC", "Normal"))
res_TNBC_Normal_Ordered <- res_TNBC_vs_Normal[order(res_TNBC_vs_Normal$pvalue),]

# 4.HER2 vs NonTNBC - Not utilized
res_HER2_vs_NonTNBC <- results(dds, contrast=c("type", "HER2", "NonTNBC"))
res_HER2_NonTNBC_Ordered <- res_HER2_vs_NonTNBC[order(res_HER2_vs_NonTNBC$pvalue),]

# 5.HER2 vs TNBC - Not utilized
res_HER2_vs_TNBC <- results(dds, contrast=c("type", "HER2", "TNBC"))
res_HER2_TNBC_Ordered <- res_HER2_vs_TNBC[order(res_HER2_vs_TNBC$pvalue),]

# 6.TNBC vs NonTNBC - Not utilized
res_TNBC_vs_NonTNBC <- results(dds, contrast=c("type", "TNBC", "NonTNBC"))
res_TNBC_NonTNBC_Ordered <- res_TNBC_vs_NonTNBC[order(res_TNBC_vs_NonTNBC$pvalue),]

#Get DE genes
diff_res_HER2_vs_Normal <-subset(res_HER2_vs_Normal,padj < 0.05)
write.csv(diff_res_HER2_vs_Normal,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\DE_HER2_vs_Normal.csv")

diff_res_NonTNBC_vs_Normal <-subset(res_NonTNBC_vs_Normal,padj < 0.05)
write.csv(diff_res_NonTNBC_vs_Normal,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\DE_NonTNBC_vs_Normal.csv")

diff_res_TNBC_vs_Normal <-subset(res_TNBC_vs_Normal,padj < 0.05)
write.csv(diff_res_TNBC_vs_Normal,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\DE_TNBC_vs_Normal.csv")

#Save up-regulated and down-regulated genes
diff_res_HER2_vs_Normal <-subset(res_HER2_vs_Normal,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
write.csv(diff_res_HER2_vs_Normal,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\HER2_vs_Normal.csv")
write.csv(res_HER2_vs_Normal,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\HER2_vs_Normal_universe.csv")

diff_res_NonTNBC_vs_Normal <-subset(res_NonTNBC_vs_Normal,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
write.csv(diff_res_NonTNBC_vs_Normal,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\NonTNBC_vs_Normal.csv")
write.csv(res_HER2_vs_Normal,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\NonTNBC_vs_Normal_universe.csv")

diff_res_TNBC_vs_Normal <-subset(res_TNBC_vs_Normal,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
write.csv(diff_res_TNBC_vs_Normal,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\TNBC_vs_Normal.csv")
write.csv(res_HER2_vs_Normal,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\TNBC_vs_Normal_universe.csv")

#Not utilized - up-regulated and down-regulated genes
diff_res_HER2_vs_NonTNBC <-subset(res_HER2_vs_NonTNBC,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
write.csv(diff_res_HER2_vs_NonTNBC,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\HER2_vs_NonTNBC.csv")
write.csv(res_HER2_vs_Normal,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\HER2_vs_NonTNBC_universe.csv")

diff_res_HER2_vs_TNBC <-subset(res_HER2_vs_TNBC,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
write.csv(diff_res_HER2_vs_TNBC,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\HER2_vs_TNBC.csv")
write.csv(res_HER2_vs_Normal,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\HER2_vs_TNBC_universe.csv")

diff_res_TNBC_vs_NonTNBC <-subset(res_TNBC_vs_NonTNBC,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
write.csv(diff_res_TNBC_vs_NonTNBC,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\TNBC_vs_NonTNBC.csv")
write.csv(res_HER2_vs_Normal,file= "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\TNBC_vs_NonTNBC_universe.csv")

#Create volcano plot.
resdata <- na.omit(as.data.frame(res_HER2_vs_Normal))
positive_genes <- subset(resdata, log2FoldChange >= 2 & padj < 0.05)
negative_genes <- subset(resdata, log2FoldChange <= -2 & padj < 0.05)
positive_genes <- positive_genes[order(positive_genes$padj), ]
negative_genes <- negative_genes[order(negative_genes$padj), ]
top_positive_genes <- head(positive_genes, 3)
top_negative_genes <- head(negative_genes, 3)
top_genes <- rbind(top_positive_genes, top_negative_genes)
resdata <- na.omit(as.data.frame(res_HER2_vs_Normal))
volcano_plot <- ggplot(resdata, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = ifelse(log2FoldChange >= 2 & padj < 0.05, "Upregulated", ifelse(log2FoldChange <= -2 & padj < 0.05, "Downregulated", "Not Significant"))), alpha = 0.6) +
  scale_color_manual(values = c('blue','gray','red'), labels = c( 'Downregulated','Not Significant', 'Upregulated')) +
  theme_minimal() +
  labs(x = 'log2 Fold Change', y = '-log10(padj)', color = 'Significance') +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "black") +  # 添加垂直虚线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") + # 添加水平虚线
  geom_text(data = top_genes, aes(label = rownames(top_genes)), 
            nudge_x = c(0, 0, 0), nudge_y = c(0, 0, 0), vjust = -1.2 ,size = 2.5)
print(volcano_plot)

#Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)

#=====================================================================================================================================================================================
#The following analyses were not reflected in the paper.

#Count outliers
par(mar=c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
W <- res_HER2_vs_Normal$stat
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic_HER2_Normal", ylab="maximum Cook's distance per gene", ylim=c(0,5), cex=.4, col=rgb(0,0,0,.3))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))

#Dispersion plot
plotDispEsts(dds)

#MA
par(mfrow=c(2,2),mar=c(2,2,1,1))
ylim <- c(-25,25)
resGA <- results(dds, lfcThreshold=.5, altHypothesis="greaterAbs")
resLA <- results(dds, lfcThreshold=.5, altHypothesis="lessAbs")
resG <- results(dds, lfcThreshold=.5, altHypothesis="greater")
resL <- results(dds, lfcThreshold=.5, altHypothesis="less")
drawLines <- function() abline(h=c(-.5,.5),col="dodgerblue",lwd=2)
plotMA(resGA, ylim=ylim); drawLines()
plotMA(resLA, ylim=ylim); drawLines()
plotMA(resG, ylim=ylim); drawLines()
plotMA(resL, ylim=ylim); drawLines()
