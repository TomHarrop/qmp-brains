#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(DESeq2)
library(data.table)
library(ggplot2)

# dds <- readRDS("output/030_deseq/dds.Rds")
# alpha <- 0.01
dds_file <- snakemake@input[["dds"]]
alpha <- snakemake@params[["alpha"]]
plot_file <- snakemake@output[["pca_plot"]]
lrt_results_file <- snakemake@output[["lrt_results"]]

# read the data
dds <- readRDS(dds_file)

# exclude low-RIN sample (s1_A), and the confusing "response 0" group (since
# there is a separate QMP control)
discard <-  c("s1_A", "s0_A", "s0_B", "s0_C")

# filter genes with low expression
row_means <- rowMeans(counts(dds))
row_max <- rowMax(counts(dds))
keep_genes <- (row_means > 5 | row_max > 10)

# make a new dds object
dds_filtered <- dds[keep_genes,
                    !rownames(colData(dds)) %in% discard]
dds_filtered$group <- droplevels(dds_filtered$group)
design(dds_filtered) <- ~ group

# run a likelihood ratio test against the group variable
dds_lrt <- DESeq(dds_filtered,
                 test = "LRT",
                 reduced = ~ 1)
res_lrt <- results(dds_lrt,
                   alpha = alpha,
                   tidy = TRUE)

# extract the results
res_dt <- data.table(res_lrt)
setorder(res_dt, padj, na.last = TRUE)

# run a pca on the passed genes
keep_lrt <- res_dt[!is.na(padj), unique(row)]
vst <- varianceStabilizingTransformation(dds_filtered[keep_lrt,],
                                         blind = TRUE)
vst_counts <- assay(vst)
pca <- prcomp(t(vst_counts), center = TRUE)

# calculate variance
pct_var <- data.table(
    comp = paste0("PC", 1:length(pca$sdev)),
    pct = 100 * pca$sdev^2 / sum(pca$sdev^2))

# plot the PCA
pca_dt <- data.table(pca$x,
                     keep.rownames = TRUE)
pca_pd <- melt(pca_dt, id.vars = "rn", variable.name = "comp")
pca_pd2 <- merge(pca_pd, pct_var)

# format labels
pca_pd2[, fac_lab := paste0(comp, " (", round(pct, 1), "%)")]
pca_pd2[, c("sample", "repl") := tstrsplit(rn, "_")]
sample_order <- c("sQMP", "s1", "s2", "s3")
pca_pd2[, sample := factor(sample, levels = sample_order)]

# plot the interesting components
pca_to_plot <- pca_pd2[pct > 5]
gp <- ggplot(pca_to_plot,
       aes(x = sample, y = value, colour = repl)) +
    theme(legend.position = c(5/6, 1/4)) +
    xlab(NULL) + ylab("Value") +
    scale_color_brewer(palette = "Set1",
                       guide = guide_legend(title = "Replicate")) +
    facet_wrap(facets = vars(fac_lab)) +
    geom_point()


# write the output
ggsave(plot_file,
       gp,
       width = 10,
       height = 7.5,
       units = "in",
       device = cairo_pdf)

fwrite(res_dt,
       lrt_results_file)

# log
sessionInfo()
