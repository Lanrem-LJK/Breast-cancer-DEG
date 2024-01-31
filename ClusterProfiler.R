library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(DOSE)
library(enrichplot)

#HER2_vs_Normal;TNBC_vs_Normal;NonTNBC_vs_Normal
#The generation of the corresponding file only requires changing the name at the input and output locations.

# 读取文件
Universe <- read.csv("C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\ego\\NonTNBC_vs_Normal_universe.csv", header = TRUE)
gene_data = read.csv("C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\DESeq2\\NonTNBC_vs_Normal.csv", header = TRUE)

# 提取第一列作为基因列表
genename <- gene_data[, 1]
UniList <- Universe[, 1]

#转换基因ID的索引
gene_ID_new <- bitr(genename, fromType = "ENSEMBL",
                    toType = c("ENTREZID", "SYMBOL"),
                    OrgDb = org.Hs.eg.db,drop = FALSE)
gene_ID_new <- gene_ID_new[!duplicated(gene_ID_new$ENSEMBL), ]

Universe_new <- bitr(UniList, fromType = "ENSEMBL",
                     toType = c("ENTREZID", "SYMBOL"),
                     OrgDb = org.Hs.eg.db,
                     drop = FALSE)
Universe_new <- distinct(Universe_new, ENSEMBL, .keep_all = TRUE)

#geneList
GeneList = gene_data[, 3]
names(GeneList) = as.character(gene_ID_new[,2])
GeneList = sort(GeneList, decreasing = TRUE)

ggo <- groupGO(gene     = gene_ID_new$ENTREZID,
               OrgDb    = org.Hs.eg.db,
               ont      = "CC",
               level    = 3,
               readable = TRUE)

ego <- enrichGO(gene          = gene_ID_new$ENSEMBL,
                universe      = Universe_new$ENSEMBL,
                keyType       = 'ENSEMBL',
                OrgDb         = org.Hs.eg.db,
                ont           = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

#Visualizing enriched terms
de <- names(GeneList)[abs(GeneList) > 2]
edo <- enrichDGN(de)
dotplot(edo, showCategory=10) + ggtitle("NonTNBC_vs_Normal")

#保存ego结果
write.table(ego, file = "C:\\Users\\16848\\Desktop\\RNA-Breast\\需要上传的数据\\ego\\ego_NonTNBC_vs_Normal.csv", sep = ",", row.names = FALSE)

