# Small molecule screen bulk RNA-seq analysis
# DESeq2 processing of the small molecule screen, GSVA scoring of MSigDB pathway and custom
# ECM signatures against matched controls, and comparison of MEK inhibition with AFOS.
# Inputs: sampleTable_screen.txt  sample metadata
#         counts_screen.csv       gene-level counts; col 1 = ensembl_id, col 2 = gene_symbol, cols 3+ = samples
#         sampleTable_AFOS.txt and counts_AFOS.csv for the MEKi vs AFOS comparison in section 5

# ── 0. Load Packages ──────────────────────────────────────────────────────────

library(DESeq2)
library(GSVA)
library(msigdbr)
library(org.Mm.eg.db)  
library(limma) 
library(matrixStats)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)

# ── 1. Run DESeq2 ─────────────────────────────────────────────────────────────

#load and format necessary files
#file directories
sampleTableFile <- "~/path/to/sampleTable_screen.txt"
countsFile <- "~/path/to/counts_screen.csv"

#read and format meta data
col_data <- read.delim(sampleTableFile, header = TRUE)

#read and format counts data
counts <- read.csv(countsFile, header = TRUE)
mouseGeneNames <- counts[,2]
count_data <- counts[, 3:50]
count_data <- as.matrix(sapply(count_data, as.numeric))
rownames(count_data) <- counts$ensembl_id
colnames(count_data) <- col_data$sampleID

#run DESeq2
dds_screen <- DESeqDataSetFromMatrix(countData = count_data, colData = col_data, design = ~ Experiment + Condition)
dds_screen <- DESeq(dds_screen)
resultsNames(dds_screen)

# ── 2. GSVA scoring of pathway signatures ─────────────────────────────────────

#pull signatures from msigdb
get_sets <- function(coll, sub = NULL) tryCatch(
  msigdbr(species = "Mus musculus", collection = coll, subcollection = sub),
  error = function(e) msigdbr(species = "Mus musculus", category = coll, subcategory = sub))
to_list <- function(df) { df <- df[!is.na(df$ensembl_gene) & df$ensembl_gene != "", ]; split(df$ensembl_gene, df$gs_name) }

hallmark <- to_list(get_sets("H"))
kegg     <- to_list(get_sets("C2", "CP:KEGG_LEGACY"))
pid      <- to_list(get_sets("C2", "CP:PID"))
c2all    <- to_list(get_sets("C2"))
c6       <- to_list(get_sets("C6"))

target_sets <- list(
  "PI3K/AKT/mTOR"    = hallmark[["HALLMARK_PI3K_AKT_MTOR_SIGNALING"]],
  "JAK-STAT"         = kegg[["KEGG_JAK_STAT_SIGNALING_PATHWAY"]],
  "JNK"              = c2all[["HAN_JNK_SINGALING_UP"]],
  "TGF-beta"         = c2all[["PLASARI_TGFB1_TARGETS_1HR_UP"]],
  "WNT"              = hallmark[["HALLMARK_WNT_BETA_CATENIN_SIGNALING"]],
  "YAP"              = c6[["CORDENONSI_YAP_CONSERVED_SIGNATURE"]],
  "p38 MAPK"         = pid[["PID_P38_ALPHA_BETA_DOWNSTREAM_PATHWAY"]],
  "AP-1"             = c2all[["MATTHEWS_AP1_TARGETS"]],
  "MEK"              = hallmark[["HALLMARK_KRAS_SIGNALING_UP"]],
  "Notch"            = c2all[["VILIMAS_NOTCH1_TARGETS_UP"]])
target_sets <- target_sets[!sapply(target_sets, is.null)]

#vst matrix, strip version, dedup
vst_mat <- assay(vst(dds_screen, blind = FALSE))
rownames(vst_mat) <- sub("\\.\\d+$", "", rownames(vst_mat))
vst_mat <- vst_mat[order(-rowMeans(vst_mat)), ]
vst_mat <- vst_mat[!duplicated(rownames(vst_mat)) & !is.na(rownames(vst_mat)), ]
target_sets <- lapply(target_sets, function(g) intersect(g, rownames(vst_mat)))
target_sets <- target_sets[lengths(target_sets) >= 5]

#run gsva
scores <- gsva(gsvaParam(vst_mat, target_sets))


#tidy and normalize each condition to its own experiment's control
scores_df <- as.data.frame(t(scores)) %>%
  tibble::rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "GeneSet", values_to = "Score") %>%
  left_join(data.frame(SampleID = col_data$sampleID, Treatment = col_data$Treatment,
                       Condition = as.character(col_data$Condition),
                       Experiment = col_data$Experiment), by = "SampleID") %>%
  group_by(GeneSet, Experiment) %>%
  mutate(Score = Score - mean(Score[Condition == "Control"])) %>% ungroup() %>%
  filter(Condition != "Control")

#plot function
plot_set <- function(set_name, group_var = "Condition") {
  df <- scores_df %>% filter(GeneSet == set_name)
  summ <- df %>% group_by(.data[[group_var]]) %>%
    summarise(m = mean(Score), sd = sd(Score), .groups = "drop") %>% arrange(.data[[group_var]])
  lvls <- unique(summ[[group_var]])
  summ[[group_var]] <- factor(summ[[group_var]], levels = lvls); df[[group_var]] <- factor(df[[group_var]], levels = lvls)
  ggplot(summ, aes(.data[[group_var]], m)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_jitter(data = df, aes(.data[[group_var]], Score),
                width = 0.15, size = 1, colour = "grey30", inherit.aes = FALSE) +
    geom_errorbar(aes(ymin = m - sd, ymax = m + sd), width = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    labs(title = set_name, x = NULL, y = "score vs matched control") +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


#draw plots
wrap_plots(lapply(names(target_sets), plot_set, group_var = "Condition"), ncol = 3)
wrap_plots(lapply(names(target_sets), plot_set, group_var = "Treatment"), ncol = 3)

# ── 3. GSVA scoring of custom signatures ──────────────────────────────────────

#custom signatures
signatures <- list(
  "Developmental ECM" = c("Adam12","Adam19","Adamts10","Adamts17","Adamts6","Adamts9","Aspn","Bmp2",
                          "Cask","Col14a1","Col24a1","Col25a1","Col27a1","Dst","Eln","Eng","Fbln1","Fbln5","Fbn2","Hmcn1",
                          "Hpse2","Htra1","Itga1","Itga9","Itgav","Itgb1","Lrp12","Ltbp2","Mfap2","Ntn4","P3h2","P4ha3",
                          "Ppib","Reck","Sdc2","Sparc","Tgfb2","Tgfbi"),
  "Fibrotic ECM" = c("Adamts12","Adamts15","Bgn","Bmp1","Col16a1","Col1a1","Col1a2","Col23a1","Col5a2",
                     "Col5a3","Col7a1","Ctsd","Fbn1","Fgf2","Fn1","Has2","Hspg2","Itga5","Lama5","Ndnf","Plec",
                     "Serpine1","Spp1","Thbs1","Timp1","Tnc"),
  "AP-1 signature" = c("Bcl3","Ca3","Cav1","Cd44","Hbegf","Stmn1","Il1rl1","Mmp10","Plaur","Mvd",
                       "Ptgs2","Spp1","Vim","Spry1","Mgst1","Flnc","Lbh"))

#convert symbols to ensembl
ecm_sets <- lapply(signatures, function(genes) {
  ids <- mapIds(org.Mm.eg.db, keys = genes, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
  unname(ids[!is.na(ids)]) })

#vst matrix, strip version, dedup
#built from dds_screen directly so this section does not depend on section 2 having run
expr_mat <- assay(vst(dds_screen, blind = FALSE))
rownames(expr_mat) <- sub("\\.\\d+$", "", rownames(expr_mat))
expr_mat <- expr_mat[order(-rowMeans(expr_mat)), ]
expr_mat <- expr_mat[!duplicated(rownames(expr_mat)) & !is.na(rownames(expr_mat)), ]
ecm_sets <- lapply(ecm_sets, function(g) intersect(g, rownames(expr_mat)))

#run gsva
ecm_scores <- gsva(gsvaParam(expr_mat, ecm_sets))

#tidy and normalize each condition to its own experiment's control
ecm_df <- as.data.frame(t(ecm_scores)) %>%
  tibble::rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "Signature", values_to = "Score") %>%
  left_join(data.frame(SampleID = col_data$sampleID, Treatment = col_data$Treatment,
                       Condition = as.character(col_data$Condition),
                       Experiment = col_data$Experiment), by = "SampleID") %>%
  group_by(Signature, Experiment) %>%
  mutate(Score = Score - mean(Score[Condition == "Control"])) %>% ungroup() %>%
  filter(Condition != "Control")

#plot function
plot_ecm <- function(sig_name, group_var = "Condition") {
  df <- ecm_df %>% filter(Signature == sig_name)
  summ <- df %>% group_by(.data[[group_var]]) %>%
    summarise(m = mean(Score), sd = sd(Score), .groups = "drop") %>% arrange(.data[[group_var]])
  lvls <- unique(summ[[group_var]])
  summ[[group_var]] <- factor(summ[[group_var]], levels = lvls); df[[group_var]] <- factor(df[[group_var]], levels = lvls)
  ggplot(summ, aes(.data[[group_var]], m)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_jitter(data = df, aes(.data[[group_var]], Score),
                width = 0.15, size = 1, colour = "grey30", inherit.aes = FALSE) +
    geom_errorbar(aes(ymin = m - sd, ymax = m + sd), width = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    labs(title = sig_name, x = NULL, y = "score vs matched control") +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#write one csv per signature for statistical testing elsewhere

#draw plots
wrap_plots(lapply(names(ecm_sets), plot_ecm, group_var = "Condition"), ncol = 2)


# ── 4. MEKi vs control analysis ───────────────────────────────────────────────

#comparison of just trametinib and the matched control
#gene symbol lookup
gene_sym <- setNames(counts$gene_symbol, counts$ensembl_id)

#subset to experiment 3 control + MEK
keep       <- col_data$Experiment == 3 & col_data$Condition %in% c("Control", "MEK")
sub_meta   <- col_data[keep, ]
sub_counts <- count_data[, sub_meta$sampleID, drop = FALSE]
sub_meta$Condition <- relevel(factor(as.character(sub_meta$Condition)), ref = "Control")

#deseq2, design ~ condition (single experiment, no batch term)
dds_mek <- DESeqDataSetFromMatrix(sub_counts, sub_meta, design = ~ Condition)
dds_mek <- dds_mek[rowSums(counts(dds_mek) >= 10) >= 3, ]
dds_mek <- DESeq(dds_mek)

#MEK vs control results
res <- results(dds_mek, contrast = c("Condition", "MEK", "Control"), alpha = 0.05)
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("ensembl_id") %>%
  mutate(gene_symbol = gene_sym[ensembl_id]) %>%
  arrange(padj)


# ── 5. MEKi vs AFOS analysis ──────────────────────────────────────────────────

#file directories
screen_sampleTable <- "~/path/to/sampleTable_screen.txt"
screen_counts      <- "~/path/to/counts_screen.csv"
afos_sampleTable   <- "~/path/to/sampleTable_AFOS.txt"
afos_counts        <- "~/path/to/counts_AFOS.csv"

#load drug screen counts (ensembl_id, gene_symbol, then samples)
screen_meta <- read.delim(screen_sampleTable, header = TRUE)
screen_raw  <- read.csv(screen_counts, header = TRUE)
screen_mat  <- as.matrix(sapply(screen_raw[, 3:ncol(screen_raw)], as.numeric))
rownames(screen_mat) <- sub("\\.\\d+$", "", screen_raw$ensembl_id)
colnames(screen_mat) <- screen_meta$sampleID

#load AFOS counts (ensembl_id, symbol, then 6 samples: 3 EGFP, 3 AFOS)
afos_meta <- read.delim(afos_sampleTable, header = TRUE)
afos_raw  <- read.csv(afos_counts, header = TRUE)
afos_mat  <- as.matrix(sapply(afos_raw[, 3:8], as.numeric))
rownames(afos_mat) <- sub("\\.\\d+$", "", afos_raw$ensembl_id)
colnames(afos_mat) <- afos_meta$sampleID

#gene symbol lookup (from screen file)
gene_sym <- setNames(screen_raw$gene_symbol, sub("\\.\\d+$", "", screen_raw$ensembl_id))

#merge on shared genes
common <- intersect(rownames(screen_mat), rownames(afos_mat))
merged_counts <- cbind(screen_mat[common, ], afos_mat[common, ])
merged_counts <- round(merged_counts); storage.mode(merged_counts) <- "integer"

#combined metadata; pool controls (drug Control + AFOS EGFP) and label EGFP -> Control
afos_cond <- ifelse(afos_meta$Condition == "EGFP", "Control", afos_meta$Condition)
merged_meta <- data.frame(
  Condition  = c(as.character(screen_meta$Condition), afos_cond),
  Experiment = c(as.character(screen_meta$Experiment), rep("AFOS", nrow(afos_meta))),
  row.names  = c(screen_meta$sampleID, afos_meta$sampleID),
  stringsAsFactors = FALSE)
merged_meta$Condition  <- relevel(factor(merged_meta$Condition), ref = "Control")
merged_meta$Experiment <- factor(merged_meta$Experiment)
stopifnot(all(rownames(merged_meta) == colnames(merged_counts)))

#deseq2 -> vst -> remove experiment batch (preserve condition)
dds_merged <- DESeqDataSetFromMatrix(merged_counts, merged_meta, design = ~ Experiment + Condition)
dds_merged <- estimateSizeFactors(dds_merged)
vsd <- vst(dds_merged, blind = FALSE)
M <- removeBatchEffect(assay(vsd), batch = merged_meta$Experiment,
                       design = model.matrix(~ Condition, merged_meta))

#condition-averaged matrix
conditions <- levels(merged_meta$Condition)
cond_mat <- sapply(conditions, function(cd) rowMeans(M[, merged_meta$Condition == cd, drop = FALSE]))
colnames(cond_mat) <- conditions

#gene sets: top 500 variable, reactome ECM org, custom dev+fibrotic ECM
rv <- rowVars(cond_mat)
row_symbols <- gene_sym[rownames(cond_mat)]

msig <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
reactome_ecm <- unique(msig$gene_symbol[msig$gs_name == "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"])

custom_ecm <- c(
  "Adam12","Adam19","Adamts10","Adamts17","Adamts6","Adamts9","Aspn","Bmp2","Cask","Col14a1",
  "Col24a1","Col25a1","Col27a1","Dst","Eln","Eng","Fbln1","Fbln5","Fbn2","Hmcn1","Hpse2","Htra1",
  "Itga1","Itga9","Itgav","Itgb1","Lrp12","Ltbp2","Mfap2","Ntn4","P3h2","P4ha3","Ppib","Reck",
  "Sdc2","Sparc","Tgfb2","Tgfbi","Adamts12","Adamts15","Bgn","Bmp1","Col16a1","Col1a1","Col1a2",
  "Col23a1","Col5a2","Col5a3","Col7a1","Ctsd","Fbn1","Fgf2","Fn1","Has2","Hspg2","Itga5","Lama5",
  "Ndnf","Plec","Serpine1","Spp1","Thbs1","Timp1","Tnc")

views <- list(
  "top 500 variable genes" = rownames(cond_mat)[order(rv, decreasing = TRUE)][1:min(500, sum(rv > 0))],
  "Reactome ECM genes"     = rownames(cond_mat)[!is.na(row_symbols) & row_symbols %in% reactome_ecm & rv > 0],
  "custom ECM signature"   = rownames(cond_mat)[!is.na(row_symbols) & row_symbols %in% custom_ecm & rv > 0])

#custom ECM genes to highlight in red on every plot
ecm_highlight <- views[["custom ECM signature"]]

#scatter: AFOS vs MEK log2FC from control, with correlations
afos_vs_mek <- function(genes, label) {
  df <- data.frame(
    AFOS_FC = cond_mat[genes, "AFOS"] - cond_mat[genes, "Control"],
    MEK_FC  = cond_mat[genes, "MEK"]  - cond_mat[genes, "Control"],
    is_ecm  = genes %in% ecm_highlight)
  pear <- cor.test(df$AFOS_FC, df$MEK_FC, method = "pearson")
  spr  <- suppressWarnings(cor.test(df$AFOS_FC, df$MEK_FC, method = "spearman"))
  lim  <- max(abs(c(df$AFOS_FC, df$MEK_FC)), na.rm = TRUE)
  ggplot(df, aes(AFOS_FC, MEK_FC)) +
    geom_hline(yintercept = 0, colour = "grey80") + geom_vline(xintercept = 0, colour = "grey80") +
    geom_point(aes(colour = is_ecm, alpha = is_ecm), size = 1.6) +
    geom_smooth(method = "lm", se = FALSE, colour = "black", linewidth = 0.5) +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick"),
                        labels = c("other", "ECM"), name = NULL) +
    scale_alpha_manual(values = c("FALSE" = 0.4, "TRUE" = 1), guide = "none") +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(title = paste0("AFOS vs MEK — ", label),
         subtitle = sprintf("Pearson r = %.3f (p = %.2e)   Spearman rho = %.3f (p = %.2e)",
                            pear$estimate, pear$p.value, spr$estimate, spr$p.value),
         x = "AFOS log2FC vs Control", y = "MEK log2FC vs Control") +
    theme_bw()
}

#draw the three scatter plots
afos_vs_mek(views[["top 500 variable genes"]], "top 500 variable genes")
afos_vs_mek(views[["Reactome ECM genes"]],     "Reactome ECM genes")
afos_vs_mek(views[["custom ECM signature"]],   "custom ECM signature")
