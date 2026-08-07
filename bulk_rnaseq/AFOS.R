# AFOS bulk RNA-seq analysis
# DESeq2 differential expression in EGFP- vs AFOS-transduced mouse pulmonary fibroblasts.
# Inputs: sampleTable_AFOS.txt  sample metadata
#         counts_AFOS.csv       gene-level counts; col 1 = ensembl_id, col 2 = gene symbol, cols 3-8 = samples

# ── 0. Load Packages ──────────────────────────────────────────────────────────

library(DESeq2)
library(ggplot2)

# ── 1. Run DESeq2 ─────────────────────────────────────────────────────────────

#load and format necessary files
#file directories
sampleTableFile <- "~/path/to/sampleTable_AFOS.txt"
countsFile      <- "~/path/to/counts_AFOS.csv"

#read meta data
col_data <- read.delim(sampleTableFile, header = TRUE)

#read and format counts data
counts         <- read.csv(countsFile, header = TRUE)
mouseGeneNames <- counts[, 2]
count_data     <- counts[, 3:8]
count_data     <- as.matrix(sapply(count_data, as.numeric))
rownames(count_data) <- counts$ensembl_id
colnames(count_data) <- col_data$sampleID

#run DESeq2
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Condition)
dds <- DESeq(dds)
resultsNames(dds)


# ── 2. Principal component analysis ───────────────────────────────────────────

vsd <- varianceStabilizingTransformation(dds, blind = TRUE)
plotPCA(vsd, intgroup = "Condition")


# ── 3. Differential expression analysis ───────────────────────────────────────

#EGFP vs AFOS
EGFP_vs_AFOS <- results(dds, contrast = c("Condition", "EGFP", "AFOS"))
EGFP_vs_AFOS <- cbind(EGFP_vs_AFOS, mouseGeneNames)
