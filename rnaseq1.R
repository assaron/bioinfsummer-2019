## ----message=FALSE-------------------------------------------------------
source("functions.R")
library(DESeq2)
library(ggplot2)
library(data.table)
library(ggrepel)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(pheatmap)
library(msigdbr)
library(fgsea)


## ------------------------------------------------------------------------
es <- read.gct("data/myo_p56_raw.gct")


## ------------------------------------------------------------------------
es


## ------------------------------------------------------------------------
head(exprs(es))


## ------------------------------------------------------------------------
head(fData(es))


## ------------------------------------------------------------------------
pData(es)


## ------------------------------------------------------------------------
str(es$treatment)


## ------------------------------------------------------------------------
es$condition <- c(rep("ShamP56", 4), rep("MIP56", 4))


## ------------------------------------------------------------------------
es$condition <- gsub("_.*", "", es$title)


## ------------------------------------------------------------------------
es$condition


## ------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(exprs(es), pData(es), ~condition)
dds


## ----fig.height=5--------------------------------------------------------
hist(exprs(es[,1]), breaks=100)


## ----fig.height=5--------------------------------------------------------
hist(log2(exprs(es))[,1], breaks=100)


## ----rlog,cache=TRUE-----------------------------------------------------
rld <- rlog(dds)
summary(assay(rld))


## ----fig.height=5--------------------------------------------------------
plotPCA(rld)


## ----ggplot2-1, fig.height=3.5, fig.width=7,dev='svg', message=FALSE-----
ggplot(data=mtcars, aes(x=mpg, y=hp, color = as.factor(gear))) + 
    geom_point()


## ----message=FALSE, fig.height=5-----------------------------------------
plotPCA(rld) + geom_text_repel(aes(label=name)) + theme_classic()



## ----filtering,cache=TRUE------------------------------------------------
outliers <- c("MIP56_Myo_3")

es.filtered <- es[, !es$title %in% outliers]

dds <- DESeqDataSetFromMatrix(exprs(es.filtered), 
                              pData(es.filtered),
                              ~condition)



## ----rlog-filtered,cache=TRUE,fig.height=5-------------------------------
rld <- rlog(dds)

plotPCA(rld) + geom_text_repel(aes(label=name)) + theme_classic()



## ------------------------------------------------------------------------
es.norm <- es.filtered
exprs(es.norm) <- assay(rld)


## ----deseq,cache=TRUE----------------------------------------------------
dds <- DESeq(dds)


## ----fig.height=5--------------------------------------------------------
plotDispEsts(dds)


## ------------------------------------------------------------------------
de <- results(dds, contrast = c("condition", "MIP56", "ShamP56"),
              cooksCutoff = F)
de


## ----message=FALSE-------------------------------------------------------
de <- data.table(ID=rownames(de), as.data.table(de))
de[order(pvalue), ]


## ----message=FALSE-------------------------------------------------------
de$entrez <- mapIds(org.Mm.eg.db, keys = de$ID,
                    keytype = "SYMBOL", column = "ENTREZID")
head(de)


## ------------------------------------------------------------------------
de.top12K <- head(de[order(baseMean, decreasing = TRUE), ][
                !duplicated(entrez) & !is.na(entrez), ], 12000)
de.top12K <- de.top12K[order(stat)]
de.top12K


## ------------------------------------------------------------------------
dir.create("work", showWarnings = FALSE)
fwrite(de.top12K, file="./work/ShamP56.vs.MIP56.de.tsv", sep="\t")


## ----pheatmap, fig.height=3.5, fig.width=7,dev='svg'---------------------
m <- matrix(rnorm(200), 20, 10)
pheatmap(m)


## ----top-heatmap, fig.height=3.5, fig.width=7,dev='svg'------------------
genes <- de.top12K[order(stat), c(head(ID), tail(ID))]
heatmap_table <- exprs(es.norm)[genes, ]
heatmap_table <- t(apply(heatmap_table, 1, scales::rescale))
pheatmap(heatmap_table, 
         cluster_rows = FALSE, cluster_cols = FALSE)



## ------------------------------------------------------------------------
stats <- setNames(de.top12K$stat, de.top12K$entrez)

str(stats)


## ----message=FALSE-------------------------------------------------------
# Hallmark pathways from MSigDB
m_df <- msigdbr(species = "Mus musculus", category = "H")
m_df
pathways <- split(m_df$entrez_gene, m_df$gs_name)
pathways <- lapply(pathways, unique)


## ----fgsea,cache=TRUE----------------------------------------------------
fr <- fgseaMultilevel(pathways, stats = stats, minSize=15, maxSize=500)
head(fr)


## ------------------------------------------------------------------------
frSig <- fr[padj < 0.01][order(pval)]
head(frSig)


## ----gseaPlot, fig.height=4, fig.width=8,dev='svg'-----------------------
plotEnrichment(pathways$`HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION`, stats) +
    ggtitle("Epythelial-mesenchymal transition")


## ----gseaTable, fig.height=4, fig.width=11,dev='svg'---------------------
plotGseaTable(pathways[fr[padj < 0.01][order(NES), pathway]], stats, fr)


## ------------------------------------------------------------------------
fr[, leadingEdge := mapIdsList(keys=leadingEdge,
                               x=org.Mm.eg.db, 
                               keytype="ENTREZID", 
                               column="SYMBOL")]
fwrite(fr, file="./work/ShamP56.vs.MIP56.gsea.tsv", 
       sep="\t", sep2=c("", " ", ""))

