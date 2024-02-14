# GO解析
setwd("~/at-go") # 任意のディレクトリに移動
getwd()

library(clusterProfiler)
library(org.At.tair.db)

table <- read.table("result.txt", header = TRUE, sep = "\t")
head(table)
all.genes <- unique(toTable(org.At.tairGO)$gene_id)
is.degs <- table$gene

all.genes.entrez <- bitr(all.genes, fromType = "TAIR", toType = "ENTREZID", OrgDb = "org.At.tair.db")
is.degs.entrez <- bitr(is.degs, fromType = "TAIR", toType = "ENTREZID", OrgDb = "org.At.tair.db")

ego <- enrichGO(gene = is.degs.entrez[, 2],
                  universe = all.genes.entrez[, 2],
                  OrgDb = "org.At.tair.db",
                  ont = "MF",
                  pAdjustMethod = "BH", pvalueCutoff = 0.1, qvalueCutoff = 0.1,
                  readable = FALSE)
res <- summary(ego)
head(res[, -8])

dotplot(ego, showCategory = 20, font = 12)
