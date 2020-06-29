library(DESeq2)
library(data.table)
library(ggplot2)
library(rtracklayer)

dds <- readRDS("output/030_deseq/dds.Rds")
gff <- "data/ref/GCF_003254395.2_Amel_HAv3.1_genomic.gff"

discard <-  c("s1_A", "s0_A", "s0_B", "s0_C")

row_means <- rowMeans(counts(dds))
row_max <- rowMax(counts(dds))
keep_genes <- (row_means > 5 | row_max > 10)

dds_filtered <- dds[keep_genes,
                    !rownames(colData(dds)) %in% discard]
dds_filtered$group <- droplevels(dds_filtered$group)

dds_lrt <- DESeq(dds_filtered, test = "LRT", reduced = ~ 1)
res_lrt <- results(dds_lrt, alpha = 0.01, tidy = TRUE)

# annotate the LRT
gr <- import.gff3(gff,
                  feature.type = "gene")
annot <- as.data.table(mcols(gr))
annot[Name == "LOC724721"]

# explore the LRT
res_lrt[order(res_lrt$padj),]
x <- data.table(subset(res_lrt, padj < 0.01))
setorder(x, padj)
fwrite(x[, .(gene_id = row,
             baseMean,
             pvalue,
             padj)], "de_genes_lrt.csv")

# plotCounts(dds_lrt, "LOC100577359", intgroup = "group")

keep_lrt <- rownames(subset(results(dds_lrt), !is.na(padj)))

vst <- varianceStabilizingTransformation(dds_filtered[keep_lrt,],
                                         blind = TRUE)

vst_counts <- assay(vst)
pca <- prcomp(t(vst_counts), center = TRUE)
pca_dt <- data.table(pca$x, keep.rownames = TRUE)
pca_pd <- melt(pca_dt, id.vars = "rn")
pca_pd[, c("sample", "repl") := tstrsplit(rn, "_")]

pca18 <- pca_pd[variable %in% paste0("PC", 1:8)]

ggplot(pca18,
       aes(x = sample, y = value, colour = repl)) +
    facet_wrap(facets = vars(variable)) +
    geom_point()

# pull out "high activation" differences
dds_hl <- copy(dds_filtered)
dds_hl$hl <- factor(ifelse(colData(dds_hl)$group %in% c(2, 3),
       yes = "high",
       no = "low"))
design(dds_hl) <- ~ hl
dds_hl <- DESeq(dds_hl)
res_hl <- results(dds_hl,
                  name = "hl_low_vs_high",
                  alpha = 0.1,
                  lfcThreshold = 0.58)
res_hl[order(abs(res_hl$log2FoldChange), na.last = FALSE),]
subset(res_hl, padj < 0.1)
