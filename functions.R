library(Biobase)

# phantasus:::read.tsv
read.tsv <- function(file, header = TRUE, sep = "\t", quote = "", comment.char = "", 
          check.names = FALSE, ...) {
    args <- list(...)
    res <- utils::read.table(file, header = header, sep = sep, 
                             quote = quote, comment.char = comment.char, check.names = check.names, 
                             stringsAsFactors = FALSE, ...)
    if ((!"row.names" %in% names(args)) && (colnames(res)[1] == 
                                            "")) {
        rownames(res) <- res[, 1]
        res[[1]] <- NULL
    }
    res
}

# phantasus:::makeAnnotated
makeAnnotated <- function(data) 
{
    meta <- data.frame(labelDescription = colnames(data))
    rownames(meta) <- colnames(data)
    methods::new("AnnotatedDataFrame", data = data, varMeta = meta)
}

# phantasus::read.gct
read.gct <- function(gct, ...) {
    meta <- readLines(gct, n = 3)
    version <- meta[1]
    size <- as.numeric(unlist(strsplit(meta[2], "\t")))

    if (grepl("^#1.3", version)) {
        # number of column annotations = number of additional rows
        ann.col <- size[4]

        # number of row annotations = number of additional columns
        ann.row <- size[3]
    } else if (grepl("^#1.2", version)) {
        ann.col <- 0
        ann.row <- 1
    } else {
        stop("Unsupported version of gct: use 1.2 or 1.3")
    }

    colNames <- unlist(strsplit(meta[3], "\t"))
    if (grepl("/", colNames[1])) {
        rowIdField <- sub("(.*)/(.*)", "\\1", colNames[1])
        colIdField <- sub("(.*)/(.*)", "\\2", colNames[1])
    } else {
        rowIdField <- "id"
        colIdField <- "id"
    }

    colNames[1] <- rowIdField

    t <- read.tsv(gct, skip = 2 + 1 + ann.col, nrows = size[1],
                  col.names = colNames,
                  row.names = NULL, header = FALSE,  ...)

    if (any(duplicated(t[,1]))) {
        warning(sprintf("duplicated row IDs: %s; they were renamed",
                        paste0(t[head(which(duplicated(t[, 1]))),1], collapse = " ")))
        rownames(t) <- make.unique(t[,1])
    } else {
        rownames(t) <- t[,1]
    }


    exp <- as.matrix(t[, (ann.row + 2):ncol(t)])

    fdata <- makeAnnotated(t[, seq_len(ann.row + 1), drop = FALSE])


    if (ann.col > 0) {
        pdata.raw <- t(read.tsv(gct, skip = 2, nrows = ann.col + 1,
                                header = FALSE, row.names=NULL))
        pdata <- data.frame(pdata.raw[seq_len(ncol(exp)) + 1 + ann.row, ,
                                      drop = FALSE],
                            stringsAsFactors = FALSE)
        colnames(pdata) <- pdata.raw[1, ]
        colnames(pdata)[1] <- colIdField
        rownames(pdata) <- colnames(exp)
        pdata <- makeAnnotated(pdata)

        res <- ExpressionSet(exp, featureData = fdata, phenoData = pdata)
    } else {
        res <- ExpressionSet(exp, featureData = fdata)
    }

    res
}

mapIdsList <- function(keys, ...) {
    keysFlat <- unlist(keys)
    keysUnique <- unique(keysFlat)
    ansUnique <- mapIds(keys=keysUnique, ...)
    ansFlat <- ansUnique[keysFlat] 
    ans <- split(ansFlat, rep(seq_along(keys), sapply(keys, length)))
    ans
}


library(data.table)
hyperTest <- function(genes, pathways, universe) {
    stopifnot(require(data.table))
    pathways <- lapply(pathways, intersect, universe)
    if (!all(genes %in% universe)) {
        stop("There are genes not from the universe, please check the input")
    }
    overlaps <- data.frame(
        q=sapply(sapply(pathways, intersect, genes), length),
        m=sapply(pathways, length),
        n=length(universe)-sapply(pathways, length),
        k=length(genes))
    
    # q-1 because we want probability of having >=q white balls
    pathways.pvals <- with(overlaps,
                           mapply(phyper, q-1, m, n, k, lower.tail = FALSE))
    res <- data.table(pathway=names(pathways), 
                      pval=pathways.pvals,
                      k=overlaps$q,
                      K=overlaps$m)
    
    res[, padj := p.adjust(pval, method="BH")]

    res <- res[order(pval),]
    res
}