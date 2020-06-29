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

BiocParallel::register(
    BiocParallel::MulticoreParam(workers = snakemake@threads))

#############
# FUNCTIONS # 
#############

ExcludeAndRunDeseq <- function(dds, keep_genes, discard) {
    # make a new dds object
    dds_filtered <- dds[keep_genes,
                        !rownames(colData(dds)) %in% discard]
    dds_filtered$group <- droplevels(dds_filtered$group)
    design(dds_filtered) <- ~ group
    
    # run DESeq2
    dds_filtered <- DESeq(dds_filtered, parallel = TRUE)
    
    # extract result combos
    group_fact <- dds_filtered$group
    group_fact <- relevel(group_fact, "control")
    groups <- rev(levels(group_fact))
    group_pairs <- combn(groups, 2)
    all_res <- apply(group_pairs, 2, function(pair) 
        results(dds_filtered,
                contrast = c("group", pair[[1]], pair[[2]]),
                alpha = alpha,
                lfcThreshold = lfccutoff,
                tidy = TRUE))
    
    # combine results
    names(all_res) <- apply(group_pairs, 2, function(pair)
        paste(pair[[2]], pair[[1]], sep = "_vs_"))
    
    res_dt <- rbindlist(all_res, idcol = "contrast")
    setorder(res_dt, row, contrast, na.last = TRUE)
    return(res_dt)
}

###########
# GLOBALS #
###########

# dds_file <- "output/030_deseq/dds.Rds"
# alpha <- 0.1
# lfccutoff <- log(1.5, 2)
dds_file <- snakemake@input[["dds"]]
alpha <- snakemake@params[["alpha"]]
lfccutoff <- snakemake@params[["lfc_cutoff"]]

########
# MAIN #
########

# read the data
dds <- readRDS(dds_file)

# exclude low-RIN sample (s1_A), and the confusing "response 0" group (since
# there is a separate QMP control)
discard_s1 <-  c("s1_A", "s1_B", "s1_C")
discard_rin <-  c("s1_A")

# filter genes with low expression
row_means <- rowMeans(counts(dds))
row_max <- rowMax(counts(dds))
keep_genes <- (row_means > 5 | row_max > 10)

no_s1 <- ExcludeAndRunDeseq(dds, keep_genes, discard = discard_s1)
no_s1a <- ExcludeAndRunDeseq(dds, keep_genes, discard = discard_rin)

fwrite(no_s1a, snakmake@output[["res_with_s1"]])
fwrite(no_s1, snakmake@output[["res_no_s1"]])

sessionInfo()

