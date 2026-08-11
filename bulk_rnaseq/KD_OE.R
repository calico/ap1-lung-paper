# AP-1 knockdown / overexpression bulk RNA-seq analysis
# DESeq2 differential expression for AP-1 family knockdowns and CREB5/ATF7 overexpression
# in mouse pulmonary fibroblasts, plus MuSiC deconvolution and ssGSEA/GSVA signature scoring.
# Inputs: sampleTable_KD_OE.txt  sample metadata
#         counts_KD_OE.csv       gene-level counts; col 1 = ensembl_id, col 2 = gene _ymbol, cols 3-45 = samples
#         fibroblasts_final.rds  annotated multiome fibroblast object used as the single-cell reference

# ── 0. Load Packages ──────────────────────────────────────────────────────────

library(DESeq2)
library(GSVA)
library(Seurat)
library(MuSiC)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyr)
library(tibble)


# ── 1. Run DESeq2 ─────────────────────────────────────────────────────────────

#load and format necessary files
#file directories
sampleTableFile <- "~/path/to/sampleTable_KD_OE.txt"
countsFile      <- "~/path/to/counts_KD_OE.csv"

#read and format meta data
col_data <- read.delim(sampleTableFile, header = TRUE)

#read and format counts data
counts         <- read.csv(countsFile, header = TRUE)
mouseGeneNames <- counts[, 2]
count_data     <- counts[, 3:45]
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

#differential expression analysis versus the matched control
contrast_list <- list(
  EGFP_vs_Creb5         = c("EGFP",      "Creb5 OE"),
  EGFP_vs_Atf7          = c("EGFP",      "Atf7 OE"),
  siControl_vs_siJun    = c("siControl", "siJun"),
  siControl_vs_siJunb   = c("siControl", "siJunb"),
  siControl_vs_siJund   = c("siControl", "siJund"),
  siControl_vs_siFos    = c("siControl", "siFos"),
  siControl_vs_siFosl1  = c("siControl", "siFosl1"),
  siControl_vs_siFosl2  = c("siControl", "siFosl2"),
  siControl_vs_siFosb   = c("siControl", "siFosb"),
  EGFP_vs_siJunb_Creb5  = c("EGFP",      "siJunb Creb5 OE"),
  EGFP_vs_siJunb_Atf7   = c("EGFP",      "siJunb Atf7 OE"),
  EGFP_vs_siFosl1_Creb5 = c("EGFP",      "siFosl1 Creb5 OE"),
  EGFP_vs_siFosl1_Atf7  = c("EGFP",      "siFosl1 Atf7 OE")
)

#run each contrast and append gene symbols
#results are held in a named list, e.g. de_results$siControl_vs_siJunb
de_results <- lapply(contrast_list, function(comparison) {
  res <- results(dds, contrast = c("Condition", comparison[1], comparison[2]), tidy = TRUE)
  cbind(res, mouseGeneNames)
})


# ── 4. MuSiC analysis ─────────────────────────────────────────────────────────

#analysis to predict the cell type identity of bulk RNA-seq data using single-cell expression profiles
#single cell expression profiles are defined by fibroblast types in fibroblasts_final
cell_type_column <- "cellType_specific"

#load in bulk rna-seq data, look at completely untreated cells
#column identities: 1 = ensembl_id, 2 = gene_id (symbol), 3+ = samples
counts        <- read.csv(countsFile, header = TRUE)
bulk_rownames <- make.unique(counts[, 2])
bulk_mtx      <- as.matrix(sapply(counts[, 3:5], as.numeric))   #use control samples
rownames(bulk_mtx) <- bulk_rownames
bulk_mtx      <- ceiling(bulk_mtx)

#load in single cell reference
sc_seurat_obj <- readRDS(file = "/path/to/fibroblasts_final.rds")

#drop cells with missing cell-type labels
keep_cells    <- !is.na(sc_seurat_obj@meta.data[[cell_type_column]])
sc_seurat_obj <- sc_seurat_obj[, keep_cells]

sc_counts  <- as.matrix(GetAssayData(sc_seurat_obj, assay = "SCT", slot = "counts"))
sc_coldata <- DataFrame(
  cellType  = sc_seurat_obj@meta.data[[cell_type_column]],
  sampleID  = colnames(sc_counts),
  row.names = colnames(sc_counts)
)
sc_sce <- SingleCellExperiment(assays = list(counts = sc_counts), colData = sc_coldata)

#match genes
common_genes <- intersect(rownames(sc_sce), rownames(bulk_mtx))
sc_sce   <- sc_sce[common_genes, ]
bulk_mtx <- bulk_mtx[common_genes, ]

#music
music_results <- music_prop(
  bulk.mtx  = bulk_mtx,
  sc.sce    = sc_sce,
  clusters  = "cellType",
  samples   = "sampleID",
  select.ct = unique(sc_sce$cellType)
)

estimated_proportions <- music_results$Est.prop.weighted


# ── 5. ssGSEA scoring ─────────────────────────────────────────────────────────

#performing ssGSEA scoring to analyze top markers of in vivo fibroblast
#populations in the wild-type, untreated mouse pulmonary fibroblasts in culture

#define cell type signatures
#top ten markers by p-adjusted value from running FindAllMarkers() on the multiome fibroblast dataset from this manuscript
cell_signatures <- list(
  "Dev_alveolar_fibroblasts"    = c("Frem1", "Slc27a6", "Col25a1", "Mfap4", "Hs6st3", "Adamts17", "Hbb-bs", "Nav2", "Adh1", "Opcml"),
  "Alveolar_myofibroblasts"     = c("Egfem1", "Nrxn3", "Rbfox1", "Tenm4", "P2ry14", "Ank2", "Zfp536", "Plxna4", "Tgfbi", "Fndc1"),
  "P7_ductal_myofibroblasts"    = c("Aspn", "Nrxn3", "Egfem1", "Rgs6", "Tgfbi", "Csmd1", "Kcnma1", "Kctd16", "Hpse2", "Nlgn1"),
  "Alveolar_fibroblasts"        = c("Inmt", "Aldh1a1", "Fmo2", "Neat1", "9530026P05Rik", "Pid1", "Abca8a", "Hsd11b1", "Lncpint", "Zbtb16"),
  "Matrix_fibroblasts"          = c("Thbs1", "Errfi1", "Egr1", "mt-Co2", "Mt1", "mt-Atp6", "mt-Nd1", "mt-Cytb", "Casp4", "mt-Co3"),
  "Adventitial_fibroblasts"     = c("Dcn", "Gpc6", "Col14a1", "Ebf1", "Adamtsl1", "Pi16", "Cdon", "Adam23", "Kcnt2", "Scara5"),
  "Adult_ductal_myofibroblasts" = c("Csmd1", "Enpp2", "Cdh4", "Lgr5", "Kctd16", "Hpse2", "Mapk4", "Cacna2d3", "Grem2", "Pcdh7"),
  "Activated_myofibroblasts"    = c("Spp1", "Runx1", "Chl1", "Slc39a14", "Mt2", "Timp1", "Kif26b", "Itga5", "Mt1", "Slc20a1")
)

#prepare the VST matrix (control samples only)
vsd     <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
unique_names      <- make.unique(mouseGeneNames)
rownames(vst_mat) <- unique_names
vst_mat_controls  <- vst_mat[, 1:3]

#run ssgsea
ssgsea_param  <- ssgseaParam(exprData = vst_mat_controls, geneSets = cell_signatures)
ssgsea_scores <- gsva(ssgsea_param)


# ── 6. GSVA scoring ───────────────────────────────────────────────────────────

#prepare the VST matrix
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)
rownames(vst_mat) <- mouseGeneNames

#keep the row with highest mean expression per gene symbol
row_means <- rowMeans(vst_mat)
keep      <- order(row_means, decreasing = TRUE)
vst_mat_dedup <- vst_mat[keep, ]
vst_mat_dedup <- vst_mat_dedup[!duplicated(rownames(vst_mat_dedup)), ]
vst_mat_dedup <- vst_mat_dedup[!is.na(rownames(vst_mat_dedup)), ]

#signatures to score
signatures <- list(
  #developmental ecm signature
  "Developmental ECM" = c("Adam12", "Adam19", "Adamts10", "Adamts17", "Adamts6", "Adamts9", "Aspn", "Bmp2", "Cask", "Col14a1", "Col24a1", "Col25a1", "Col27a1", "Dst", "Eln", "Eng", "Fbln1", "Fbln5", "Fbn2", "Hmcn1", "Hpse2", "Htra1", "Itga1", "Itga9", "Itgav", "Itgb1", "Lrp12", "Ltbp2", "Mfap2", "Ntn4", "P3h2", "P4ha3", "Ppib", "Reck", "Sdc2", "Sparc", "Tgfb2", "Tgfbi"),
  #fibrotic ecm signature
  "Fibrotic ECM"      = c("Adamts12","Adamts15","Bgn","Bmp1","Col16a1","Col1a1","Col1a2","Col23a1","Col5a2","Col5a3","Col7a1","Ctsd","Fbn1","Fgf2","Fn1","Has2","Hspg2","Itga5","Lama5","Ndnf","Plec","Serpine1","Spp1","Thbs1","Timp1","Tnc"),
  #AP-1 signature
  "AP-1 signature"    = c("Bcl3","Ca3","Cav1","Cd44","Hbegf","Stmn1","Il1rl1","Mmp10","Plaur","Mvd","Ptgs2","Spp1","Vim","Spry1","Mgst1","Flnc","Lbh"))

#run gsva on all three signatures at once
gsva_results <- gsva(gsvaParam(vst_mat_dedup, signatures))

#reference group for each condition: single knockdowns are compared to siControl,
#overexpression and combined knockdown/overexpression conditions to EGFP.
#the order of this vector also sets the left-to-right order of the bars.
ref_map <- c(
  "siJun"            = "siControl",
  "siJunb"           = "siControl",
  "siJund"           = "siControl",
  "siFos"            = "siControl",
  "siFosb"           = "siControl",
  "siFosl1"          = "siControl",
  "siFosl2"          = "siControl",
  "Creb5 OE"         = "EGFP",
  "Atf7 OE"          = "EGFP",
  "siJunb Creb5 OE"  = "EGFP",
  "siJunb Atf7 OE"   = "EGFP",
  "siFosl1 Creb5 OE" = "EGFP",
  "siFosl1 Atf7 OE"  = "EGFP")
condition_order <- names(ref_map)

#tidy the scores and attach the condition of each sample
gsva_long <- as.data.frame(t(gsva_results)) %>%
  tibble::rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "Signature", values_to = "Score") %>%
  left_join(data.frame(SampleID  = as.character(col_data$sampleID),
                       Condition = as.character(col_data$Condition)), by = "SampleID")

#mean score of each reference group, per signature
ref_means <- gsva_long %>%
  dplyr::filter(Condition %in% unique(ref_map)) %>%
  group_by(Signature, Reference = Condition) %>%
  summarise(ref_mean = mean(Score), .groups = "drop")

#center every condition on its own reference; the reference conditions and the
#untreated Control samples are not plotted
gsva_df <- gsva_long %>%
  dplyr::filter(Condition %in% names(ref_map)) %>%
  mutate(Reference = unname(ref_map[Condition])) %>%
  left_join(ref_means, by = c("Signature", "Reference")) %>%
  mutate(Score     = Score - ref_mean,
         Condition = factor(Condition, levels = condition_order)) %>%
  dplyr::select(SampleID, Condition, Reference, Signature, Score)

#plot function
plot_signature <- function(sig_name) {
  df <- gsva_df %>% dplyr::filter(Signature == sig_name)
  summ <- df %>% group_by(Condition) %>%
    summarise(m = mean(Score), sd = sd(Score), .groups = "drop")
  ggplot(summ, aes(Condition, m)) +
    geom_col(fill = "steelblue", width = 0.7) +
    geom_jitter(data = df, aes(Condition, Score),
                width = 0.15, size = 1, colour = "grey30", inherit.aes = FALSE) +
    geom_errorbar(aes(ymin = m - sd, ymax = m + sd), width = 0.25) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    labs(title = sig_name, x = NULL, y = "score vs matched control") +
    theme_bw(base_size = 9) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#draw plots, one panel per signature actually scored
wrap_plots(lapply(rownames(gsva_results), plot_signature), ncol = 3)
