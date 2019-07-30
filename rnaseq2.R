## ----message=FALSE-------------------------------------------------------
source("functions.R")
library(DESeq2)
library(ggrepel)
library(pheatmap)
library(msigdbr)


## ------------------------------------------------------------------------

es <- read.gct("data/gse95755.gct")


## ------------------------------------------------------------------------
es <- es[, es$`cell type` == "Cardiomyocyte"]
es$condition <- gsub("_.*", "", es$title)
es


## ----deseq-myo, cache=TRUE, message=FALSE--------------------------------
dds <- DESeqDataSetFromMatrix(exprs(es), pData(es), ~condition)
rld <- rlog(dds)
es.norm <- es
exprs(es.norm) <- assay(rld)


## ----message=FALSE, fig.height=7-----------------------------------------
plotPCA(rld) + geom_text_repel(aes(label=name)) + theme_classic()


## ------------------------------------------------------------------------
apply(exprs(es), 2, sum)


## ------------------------------------------------------------------------
es.norm.top <- es.norm[head(order(apply(exprs(es.norm), 1, mean), 
                                  decreasing = TRUE), 6000), ]
scaledExprs <- t(scale(t(exprs(es.norm.top))))
set.seed(42)
km <- kmeans(scaledExprs, 12)
str(km)


## ----fig.height=4,dev='svg'----------------------------------------------

pheatmap(km$centers, 
      cluster_cols = FALSE, cluster_rows = FALSE)


## ----fig.height=4,dev='svg'----------------------------------------------
mat <- exprs(es.norm.top)[order(km$cluster, sample.int(nrow(es.norm.top))), ]
mat <- t(apply(mat, 1, scales::rescale))
grouping <- ceiling(seq_len(nrow(mat)) / nrow(mat) * 1000)
aggr <- Matrix.utils::aggregate.Matrix(mat, groupings=grouping, fun="mean")
rownames(aggr) <- rownames(mat)[!duplicated(grouping)]
annotation_row <- data.frame(cluster=paste0("c", km$cluster), row.names = names(km$cluster))



## ----fig.height=4, dev='svg'---------------------------------------------
pheatmap(aggr, 
         cluster_rows = FALSE, cluster_cols = FALSE, 
         annotation_row = annotation_row,
         show_rownames = FALSE)


## ----fig.height=4,dev='svg'----------------------------------------------
genes <- names(which(km$cluster==5))
heatmap_table <- exprs(es.norm)[genes, ]
heatmap_table <- t(apply(heatmap_table, 1, scales::rescale))
pheatmap(heatmap_table, show_rownames = FALSE,
         cluster_rows = FALSE, cluster_cols = FALSE)


## ------------------------------------------------------------------------
m_df <- msigdbr(species = "Mus musculus", category = "H")
pathways <- split(m_df$gene_symbol, m_df$gs_name)
pathways <- lapply(pathways, unique)
universe <- rownames(es.norm.top)


## ------------------------------------------------------------------------
hyp <- hyperTest(genes, pathways, universe)
head(hyp)

