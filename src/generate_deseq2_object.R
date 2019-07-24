#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]],
            open = "wt")
sink(log,
     type = "message")
sink(log,
     append = TRUE,
     type = "output")

library(data.table)
library(DESeq2)
library(rtracklayer)
library(tximport)


file_list <- snakemake@input[["quant_files"]]
gff <- snakemake@input[["gff"]]
dds_file <- snakemake@output[[1]]

# read gff using GenomicRanges
gr <- import.gff3(gff,
                  feature.type = c("exons", "CDS", "mRNA", "gene"))

# extract a data.frame of tx to gene
mrnas <- gr[gr$type == "mRNA",]
mrna_dt <- as.data.table(mcols(mrnas))
tx2gene <- data.frame(mrna_dt[, .(TXNAME = Name,
                                  GENEID = as.character(gene))])

# get file list
# file_list <- list.files("output/020_salmon",
#                         full.names = TRUE,
#                         recursive = TRUE,
#                         pattern = "quant.sf")
names(file_list) <- gsub(".*/", "", dirname(file_list))

# import quant files
txi <- tximport(file_list, type = "salmon", tx2gene = tx2gene)

# which stuff isn't in tx2gene?
read_list <- lapply(file_list, fread)
reads <- rbindlist(read_list, idcol = "sample")
missing_genes <- reads[!Name %in% tx2gene$TXNAME, unique(Name)]

# generate col_data
col_data <- data.table(samplename = names(file_list))
col_data[, replicate := gsub(".*([[:upper:]])$", "\\1", samplename)]
col_data[, group := as.character(as.numeric(
    gsub("^s([[:digit:]]+)_.*", "\\1", samplename)))]
col_data[grepl("QMP", samplename), group := "control"]

# generate DESeq object
dds <- DESeqDataSetFromTximport(
    txi,
    colData = data.frame(col_data, row.names = 'samplename'),
    design = ~ group
)

saveRDS(dds, dds_file)

# log
sessionInfo()

