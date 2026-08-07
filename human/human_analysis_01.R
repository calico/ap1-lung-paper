# Analysis of published human lung datasets
# ECM signature scoring and AP-1 transcription factor activity in
# published human scRNA-seq (Habermann GSE135893, Adams GSE136831, Tsukui GSE132771),
# sn-multiome (Valenzi GSE214085) and whole-lung microarray (GSE47460) datasets.
# Inputs: Seurat objects saved as .rds for each published dataset; GSE47460 is downloaded from GEO.

# ── 0. Load Packages ──────────────────────────────────────────────────────────

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(Seurat)
library(Signac)
library(harmony)
library(dorothea)
library(decoupleR)
library(viper)
library(BSgenome.Hsapiens.UCSC.hg19)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(chromVAR)
library(GEOquery)
library(Biobase)
library(limma)
library(GSVA)

# ── 1. GSE135893 scRNA-seq dataset analysis ───────────────────────────────────

#Habermann et al dataset
#initialize Seurat object from publicly available data from GSE135893
Habermann_seurat <- readRDS(file = "/path/to/Habermann_seurat") 

#subset the mesenchymal clusters
Idents(Habermann_seurat) <- "population"
mesenchymal_Habermann <- subset(x = Habermann_seurat, idents = "Mesenchymal")

#subset the fibroblast clusters
Idents(mesenchymal_Habermann) <- "celltype"
fibroblasts_Habermann <- subset(x = mesenchymal_Habermann, idents = c("PLIN2+ Fibroblasts","Fibroblasts","HAS1 High Fibroblasts","Myofibroblasts"))

#process RNA assay
DefaultAssay(fibroblasts_Habermann) <- "RNA"
fibroblasts_Habermann <- NormalizeData(fibroblasts_Habermann, scale.factor = 10000)
fibroblasts_Habermann <- FindVariableFeatures(fibroblasts_Habermann, selection.method = "mean.var.plot", nfeatures = 2000)
fibroblasts_Habermann <- ScaleData(fibroblasts_Habermann, features = rownames(fibroblasts_Habermann))

#add ECM signatures converted to human gene names
matrix_sig_final_human <- c("ADAMTS12","ADAMTS15","BGN","BMP1","COL16A1","COL1A1","COL1A2","COL23A1","COL5A2","COL5A3","COL7A1","CTSD","FBN1","FGF2","FN1","HAS2","HSPG2","ITGA5","LAMA5","NDNF","PLEC","SERPINE1","SPP1","THBS1","TIMP1","TNC")
all_mfb_ecm_genes_human <- c("ADAM12","ADAM19","ADAMTS10","ADAMTS17","ADAMTS6","ADAMTS9","ASPN","BMP2","CASK","COL14A1","COL24A1","COL25A1","COL27A1","DST","ELN","ENG","FBLN1","FBLN5","FBN2","HMCN1","HPSE2","HTRA1","ITGA1","ITGA9","ITGAV","ITGB1","LRP12","LTBP2","MFAP2","NTN4","P3H2","P4HA3","PPIB","RECK","SDC2","SPARC","TGFB2","TGFBI")

fibroblasts_Habermann <- AddModuleScore(fibroblasts_Habermann, features = list(matrix_sig_final_human), name = "matrix_sig_final_human")
fibroblasts_Habermann <- AddModuleScore(fibroblasts_Habermann, features = list(all_mfb_ecm_genes_human), name = "all_mfb_ecm_genes_human")

#pseudobulk fibrotic and developmental ecm scores by sample + diagnosis
meta_df <- fibroblasts_Habermann@meta.data
robust_scores <- meta_df %>%
  group_by(orig.ident, Diagnosis) %>%
  summarize(
    mean_matrix_sig = mean(matrix_sig_final_human1, na.rm = TRUE),
    mean_mfb_ecm    = mean(all_mfb_ecm_genes_human1, na.rm = TRUE),
    n_cells         = n(),
    .groups         = "drop"
  ) %>%
  filter(n_cells >= 20)

#subset control and IPF samples only
Idents(fibroblasts_Habermann) <- "Diagnosis"
fibroblasts2_Habermann <- subset(x = fibroblasts_Habermann, idents = c("Control","IPF"))

#process RNA assay
DefaultAssay(fibroblasts2_Habermann) <- "RNA"
fibroblasts2_Habermann <- NormalizeData(fibroblasts2_Habermann, scale.factor = 10000)
fibroblasts2_Habermann <- FindVariableFeatures(fibroblasts2_Habermann, selection.method = "mean.var.plot", nfeatures = 2000)
fibroblasts2_Habermann <- ScaleData(fibroblasts2_Habermann, features = rownames(fibroblasts2_Habermann))

#run DoRothEA analysis
data("dorothea_hs", package = "dorothea")
regulons <- dorothea_hs %>% filter(confidence %in% c("A", "B", "C"))
expr <- as.matrix(GetAssayData(fibroblasts2_Habermann, slot = "scale.data", assay = "RNA"))
tf_activities <- run_viper(mat = expr, network = regulons, .source = "tf", .target = "target", .mor = "mor")

tf_matrix <- tf_activities %>%
  dplyr::select(condition, source, score) %>%
  tidyr::pivot_wider(names_from = source, values_from = score) %>%
  tibble::column_to_rownames("condition")

tf_matrix <- tf_matrix[colnames(fibroblasts2_Habermann), ]
fibroblasts2_Habermann <- AddMetaData(fibroblasts2_Habermann, metadata = tf_matrix)

#pseudobulk DoRothEA scores by sample + diagnosis
meta_df <- fibroblasts2_Habermann@meta.data
robust_ap1_scores <- meta_df %>%
  group_by(orig.ident, Diagnosis) %>%
  summarize(
    mean_JUN   = mean(JUN, na.rm = TRUE),
    mean_JUNB  = mean(JUNB, na.rm = TRUE),
    mean_JUND  = mean(JUND, na.rm = TRUE),
    mean_FOS   = mean(FOS, na.rm = TRUE),
    mean_FOSL1 = mean(FOSL1, na.rm = TRUE),
    mean_FOSL2 = mean(FOSL2, na.rm = TRUE),
    n_cells    = n(),
    .groups    = "drop"
  ) %>%
  filter(n_cells >= 20)

#pseudobulk fibrotic ecm signature by sample + cell type + diagnosis
meta_df <- fibroblasts2_Habermann@meta.data
robust_scores_celltype <- meta_df %>%
  group_by(orig.ident, celltype, Diagnosis) %>%
  summarize(
    mean_matrix_sig = mean(matrix_sig_final_human1, na.rm = TRUE),
    n_cells         = n(),
    .groups         = "drop"
  ) %>%
  filter(n_cells >= 20)


# ── 2. GSE136831 scRNA-seq dataset analysis ───────────────────────────────────

#Adams et al dataset
#initialize Seurat object from publicly available data from GSE136831
Adams_seurat <- readRDS(file = "/path/to/Adams_seurat") 

#subset the stromal clusters
Idents(Adams_seurat) <- "CellType_Category"
stromal_Adams <- subset(x = Adams_seurat, idents = "Stromal")

#subset the fibroblast clusters
Idents(stromal_Adams) <- "Manuscript_Identity"
fibroblast_Adams <- subset(x = stromal_Adams, idents = c("Fibroblast","Myofibroblast"))

#process RNA data
fibroblast_Adams <- SCTransform(fibroblast_Adams, vars.to.regress = c("nCount_RNA"), verbose = FALSE)
DefaultAssay(fibroblast_Adams) <- "SCT"
fibroblast_Adams <- RunPCA(fibroblast_Adams, npcs = 50)

#SCT assay integration
DefaultAssay(fibroblast_Adams) <- "SCT"
fibroblast_Adams <- RunHarmony(fibroblast_Adams, group.by.vars = c("Library_Identity"), reduction = "pca",
                               dims.use = NULL, theta = 2, lambda = NULL, sigma = 0.1,
                               nclust = NULL, tau = 0, block.size = 0.05, max.iter.harmony = 10,
                               max.iter.cluster = 20, epsilon.cluster = 1e-05,
                               epsilon.harmony = 1e-04, plot_convergence = FALSE, verbose = TRUE,
                               reference_values = NULL, reduction.save = "harmonyBatch",
                               assay.use = "SCT", project.dim = TRUE)

#clustering
fibroblast_Adams <- RunUMAP(fibroblast_Adams, reduction = "pca", dims = 1:7)
fibroblast_Adams <- FindNeighbors(fibroblast_Adams, reduction = "pca", dims = 1:7, k.param = 50)
fibroblast_Adams <- FindClusters(fibroblast_Adams, resolution = 0.4)

#process RNA assay
DefaultAssay(fibroblast_Adams) <- "RNA"
fibroblast_Adams <- NormalizeData(fibroblast_Adams, scale.factor = 10000)
fibroblast_Adams <- FindVariableFeatures(fibroblast_Adams, selection.method = "mean.var.plot", nfeatures = 2000)
fibroblast_Adams <- ScaleData(fibroblast_Adams, features = rownames(fibroblast_Adams))

#add ecm signatures
fibroblast_Adams <- AddModuleScore(fibroblast_Adams, features = list(matrix_sig_final_human), name = "matrix_sig_final_human")
fibroblast_Adams <- AddModuleScore(fibroblast_Adams, features = list(all_mfb_ecm_genes_human), name = "all_mfb_ecm_genes_human")

#pseudobulk fibrotic ecm signature by sample + diagnosis
meta_df_adams <- fibroblast_Adams@meta.data
robust_scores_adams <- meta_df_adams %>%
  group_by(Subject_Identity, Disease_Identity) %>% 
  summarize(
    mean_matrix_sig = mean(matrix_sig_final_human1, na.rm = TRUE),
    n_cells         = n(),
    .groups         = "drop"
  ) %>%
  filter(n_cells >= 20)

#pseudobulk fibrotic ecm signature by sample + cell type + diagnosis
meta_df_adams <- fibroblast_Adams@meta.data
robust_scores_adams <- meta_df_adams %>%
  group_by(Subject_Identity, Disease_Identity, Manuscript_Identity) %>% 
  summarize(
    mean_matrix_sig = mean(matrix_sig_final_human1, na.rm = TRUE),
    n_cells         = n(),
    .groups         = "drop"
  ) %>%
  filter(n_cells >= 20)


# ── 3. GSE132771 scRNA-seq dataset analysis ───────────────────────────────────

#Tsukui et al dataset
#initialize Seurat object from publicly available data from GSE132771
Tsukui_human_seurat <- readRDS(file = "/path/to/Tsukui_human_seurat") 

#add disease labels
#add metadata to combine samples
Idents(Tsukui_human_seurat) <- "orig.ident"
Ctrl1 <- "SampleGSM3891620"
Ctrl2 <- "SampleGSM3891622"
Ctrl3 <- "SampleGSM3891624"
IPF1 <- "SampleGSM3891626"
IPF2 <- "SampleGSM3891628"
IPF3 <- "SampleGSM3891630"

Ctrl <- c(Ctrl1, Ctrl2, Ctrl3)
IPF <- c(IPF1, IPF2, IPF3)

Tsukui_human_seurat@meta.data[which(Tsukui_human_seurat@meta.data$orig.ident %in% Ctrl1), "sampleID"] <- "Ctrl1"
Tsukui_human_seurat@meta.data[which(Tsukui_human_seurat@meta.data$orig.ident %in% Ctrl2), "sampleID"] <- "Ctrl2"
Tsukui_human_seurat@meta.data[which(Tsukui_human_seurat@meta.data$orig.ident %in% Ctrl3), "sampleID"] <- "Ctrl3"
Tsukui_human_seurat@meta.data[which(Tsukui_human_seurat@meta.data$orig.ident %in% IPF1), "sampleID"] <- "IPF1"
Tsukui_human_seurat@meta.data[which(Tsukui_human_seurat@meta.data$orig.ident %in% IPF2), "sampleID"] <- "IPF2"
Tsukui_human_seurat@meta.data[which(Tsukui_human_seurat@meta.data$orig.ident %in% IPF3), "sampleID"] <- "IPF3"
Tsukui_human_seurat@meta.data[which(Tsukui_human_seurat@meta.data$orig.ident %in% Ctrl), "condition"] <- "Ctrl"
Tsukui_human_seurat@meta.data[which(Tsukui_human_seurat@meta.data$orig.ident %in% IPF), "condition"] <- "IPF"

#process RNA dataset
Tsukui_human_seurat_scaled <- SCTransform(Tsukui_human_seurat, vars.to.regress = c("nCount_RNA"), verbose = FALSE)
DefaultAssay(Tsukui_human_seurat_scaled) <- "SCT"
Tsukui_human_seurat_scaled <- RunPCA(Tsukui_human_seurat_scaled, npcs = 50)

#SCT assay integration
Tsukui_human_seurat_scaled <- RunHarmony(Tsukui_human_seurat_scaled, group.by.vars = c("sampleID"), reduction = "pca",
                                         dims.use = NULL, theta = 2, lambda = NULL, sigma = 0.1,
                                         nclust = NULL, tau = 0, block.size = 0.05, max.iter.harmony = 10,
                                         max.iter.cluster = 20, epsilon.cluster = 1e-05,
                                         epsilon.harmony = 1e-04, plot_convergence = FALSE, verbose = TRUE,
                                         reference_values = NULL, reduction.save = "harmonyBatch",
                                         assay.use = "SCT", project.dim = TRUE)

#clustering
Tsukui_human_seurat_scaled <- RunUMAP(Tsukui_human_seurat_scaled, reduction = "harmonyBatch", dims = 1:50)
Tsukui_human_seurat_scaled <- FindNeighbors(Tsukui_human_seurat_scaled, reduction = "harmonyBatch", dims = 1:50, k.param = 10)
Tsukui_human_seurat_scaled <- FindClusters(Tsukui_human_seurat_scaled, resolution = 0.6)

#subset the COL1A1+ clusters
Idents(Tsukui_human_seurat_scaled) <- "seurat_clusters"
fibroblasts_human_Tsukui <- subset(x = Tsukui_human_seurat_scaled, idents = c("0","9","18","6","12","2","10","5","13"))

#process RNA assay
DefaultAssay(fibroblasts_human_Tsukui) <- "RNA"
fibroblasts_human_Tsukui <- NormalizeData(fibroblasts_human_Tsukui, scale.factor = 10000)
fibroblasts_human_Tsukui <- FindVariableFeatures(fibroblasts_human_Tsukui, selection.method = "mean.var.plot", nfeatures = 2000)
fibroblasts_human_Tsukui <- ScaleData(fibroblasts_human_Tsukui, features = rownames(fibroblasts_human_Tsukui))

#add ecm signatures
fibroblasts_human_Tsukui <- AddModuleScore(fibroblasts_human_Tsukui, features = list(matrix_sig_final_human), name = "matrix_sig_final_human")
fibroblasts_human_Tsukui <- AddModuleScore(fibroblasts_human_Tsukui, features = list(all_mfb_ecm_genes_human), name = "all_mfb_ecm_genes_human")

#pseudobulk by sample + disease
meta_df_tsukui <- fibroblasts_human_Tsukui@meta.data
robust_scores_tsukui <- meta_df_tsukui %>%
  group_by(sampleID) %>% 
  summarize(
    mean_matrix_sig = mean(matrix_sig_final_human1, na.rm = TRUE),
    mean_mfb_ecm    = mean(all_mfb_ecm_genes_human1, na.rm = TRUE),
    n_cells         = n(),
    .groups         = "drop"
  ) %>%
  filter(n_cells >= 20)


# ── 4. GSE214085 snRNA-seq and snATAC-seq dataset analysis ────────────────────

#Valenzi et al dataset
#initialize RDS objects from publicly available data from GSE214085
Valenzi_seurat <- readRDS("/path/to/Valenzi_seurat")    #snRNA-seq experiment
Valenzi_atac <- readRDS("/path/to/Valenzi_atac")     #snATAC-seq experiment

#process snRNA-seq dataset
#qc and filtering
Valenzi_seurat[["pct_mt"]] <- PercentageFeatureSet(Valenzi_seurat, pattern = "^MT-")
Valenzi_seurat <- subset(Valenzi_seurat, nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 25000 & pct_mt < 20)

#process RNA assay
Valenzi_seurat <- NormalizeData(Valenzi_seurat)
Valenzi_seurat <- FindVariableFeatures(Valenzi_seurat, nfeatures = 3000)
Valenzi_seurat <- ScaleData(Valenzi_seurat)
Valenzi_seurat <- RunPCA(Valenzi_seurat, npcs = 50)
Valenzi_seurat <- RunHarmony(Valenzi_seurat, group.by.vars = "sample_id", dims.use = 1:30)

#clustering
Valenzi_seurat <- FindNeighbors(Valenzi_seurat, reduction = "harmony", dims = 1:30)
Valenzi_seurat <- FindClusters(Valenzi_seurat, resolution = 0.5)
Valenzi_seurat <- RunUMAP(Valenzi_seurat, reduction = "harmony", dims = 1:30, reduction.name = "umap_rna")

#annotate fibroblast cluster
Valenzi_seurat$cell_type <- ifelse(Valenzi_seurat$seurat_clusters == "3",
                               "Fibroblast", "Other")

#process snATAC-seq dataset
#qc and filtering
Valenzi_atac <- NucleosomeSignal(Valenzi_atac)
Valenzi_atac <- TSSEnrichment(Valenzi_atac, fast = FALSE)
Valenzi_atac$FRiP <- Valenzi_atac$peak_region_fragments / Valenzi_atac$passed_filters
Valenzi_atac <- subset(Valenzi_atac, nCount_ATAC > 500 & nCount_ATAC < 30000 & nucleosome_signal < 2 & TSS.enrichment > 2 & blacklist_region_fragments < 500)

#process ATAC assay
Valenzi_atac <- RunTFIDF(Valenzi_atac)
Valenzi_atac <- FindTopFeatures(Valenzi_atac, min.cutoff = "q5")
Valenzi_atac <- RunSVD(Valenzi_atac)
lsi_dims <- 2:30
Valenzi_atac <- RunHarmony(Valenzi_atac, group.by.vars = "sample_id",
                           reduction.use = "lsi", dims.use = lsi_dims,
                           assay.use = "ATAC", project.dim = FALSE)

#clustering
Valenzi_atac <- FindNeighbors(Valenzi_atac, reduction = "harmony", dims = 1:length(lsi_dims))
Valenzi_atac <- FindClusters(Valenzi_atac, resolution = 0.5, algorithm = 3)
Valenzi_atac <- RunUMAP(Valenzi_atac, reduction = "harmony", dims = 1:length(lsi_dims), reduction.name = "umap_atac")

#compute gene activity
gene_activity <- GeneActivity(Valenzi_atac)
Valenzi_atac[["RNA_activity"]] <- CreateAssayObject(counts = gene_activity)
Valenzi_atac <- NormalizeData(Valenzi_atac, assay = "RNA_activity", normalization.method = "LogNormalize", scale.factor = median(Valenzi_atac$nCount_ATAC))

#find transfer anchors
DefaultAssay(Valenzi_seurat) <- "RNA"
DefaultAssay(Valenzi_atac)   <- "RNA_activity"

transfer_anchors <- FindTransferAnchors(
  reference           = Valenzi_seurat,
  query               = Valenzi_atac,
  features            = VariableFeatures(Valenzi_seurat),
  reference.reduction = "pca",
  query.assay         = "RNA_activity",
  reduction           = "pcaproject",
  dims                = 1:29
)

predictions <- TransferData(
  anchorset        = transfer_anchors,
  refdata          = Valenzi_seurat$cell_type,
  weight.reduction = Valenzi_atac[["lsi"]],
  k.weight         = 10,
  dims             = lsi_dims 
)

Valenzi_atac <- AddMetaData(Valenzi_atac, metadata = predictions)
Valenzi_atac$predicted_cell_type <- Valenzi_atac$predicted.id
Valenzi_atac$transfer_confident  <- Valenzi_atac$prediction.score.max > 0.5

#visualize
DimPlot(Valenzi_atac, group.by = "predicted_cell_type", reduction = "umap_atac",
        label = TRUE, repel = TRUE) + NoLegend() + ggtitle("ATAC — predicted cell type")
FeaturePlot(Valenzi_atac, features = "prediction.score.max", reduction = "umap_atac",
            min.cutoff = 0.2) + scale_color_viridis_c() + ggtitle("Label transfer confidence")

#assign fibroblast label to atac dataset
Valenzi_atac$cell_type <- ifelse(Valenzi_atac$predicted_cell_type %in% c("7"), "Fibroblast", "Other")

#subset fibroblasts
Idents(Valenzi_atac) <- "predicted_cell_type"
Valenzi_atac_fibroblasts <- subset(Valenzi_atac, idents = "Fibroblast")

#process atac dataset
DefaultAssay(Valenzi_atac_fibroblasts) <- "ATAC"
Valenzi_atac_fibroblasts <- RunTFIDF(Valenzi_atac_fibroblasts)
Valenzi_atac_fibroblasts <- FindTopFeatures(Valenzi_atac_fibroblasts, min.cutoff = "q5")
Valenzi_atac_fibroblasts <- RunSVD(Valenzi_atac_fibroblasts)
Valenzi_atac_fibroblasts <- RunHarmony(Valenzi_atac_fibroblasts, group.by.vars = "sample_id", reduction.use = "lsi", dims.use = 2:10, assay.use = "ATAC", project.dim = FALSE)
Valenzi_atac_fibroblasts <- FindNeighbors(Valenzi_atac_fibroblasts, reduction = "harmony", dims = 2:9)
Valenzi_atac_fibroblasts <- FindClusters(Valenzi_atac_fibroblasts, resolution = 0.3)
Valenzi_atac_fibroblasts <- RunUMAP(Valenzi_atac_fibroblasts, reduction = "harmony", dims = 2:9, reduction.name = "umap_fib")

#run chromvar on atac fibroblasts
DefaultAssay(Valenzi_atac_fibroblasts) <- "ATAC"
genome <- BSgenome.Hsapiens.UCSC.hg19
pwm_set <- getMatrixSet(JASPAR2020, opts = list(collection = "CORE", tax_group = "vertebrates", matrixtype = "PWM"))

Valenzi_atac_fibroblasts <- AddMotifs(object = Valenzi_atac_fibroblasts, genome = genome, pfm = pwm_set)
Valenzi_atac_fibroblasts <- RunChromVAR(object = Valenzi_atac_fibroblasts, genome = genome)

#differential peak accessibility analysis
Idents(Valenzi_atac_fibroblasts) <- "disease"
DefaultAssay(Valenzi_atac_fibroblasts) <- "ATAC"
da_peaks <- FindMarkers(Valenzi_atac_fibroblasts, ident.1 = "IPF", ident.2 = "Control", test.use   = "LR", latent.vars = "nCount_ATAC", min.pct    = 0.05)


# ── 5. Save objects ───────────────────────────────────────────────────────────

saveRDS(fibroblasts_Habermann,  "/path/to/fibroblasts_Habermann.rds")
saveRDS(fibroblasts2_Habermann,  "/path/to/fibroblasts2_Habermann.rds")
saveRDS(fibroblast_Adams,    "/path/to/fibroblast_Adams.rds")
saveRDS(fibroblasts_human_Tsukui,  "/path/to/fibroblasts_human_Tsukui.rds")
saveRDS(Valenzi_atac_fibroblasts, "/path/to/Valenzi_atac_fibroblasts.rds")


# ── 6. GSE47460 microarray dataset analysis ───────────────────────────────────

#analysis of publicly available whole-lung microarray data from GSE47460

#define gene signatures
fibrosis_human <- toupper(c("Adamts12","Adamts15","Bgn","Bmp1","Col16a1","Col1a1","Col1a2","Col23a1","Col5a2","Col5a3","Col7a1","Ctsd","Fbn1","Fgf2","Fn1","Has2","Hspg2","Itga5","Lama5","Ndnf","Plec","Serpine1","Spp1","Thbs1","Timp1","Tnc"))
afos_human <- toupper(c("Ephx1","Col11a1","Gpx3","Mgp","Timp3","Mylk","Ehd3","Cdo1","Ifit2","Tgtp1"))
fibroblast_genes <- c("PDGFRA","NPNT","CES1D","SLC7A10","PI16","CCL11","IL33","ADH7","DCN","HHIP","ASPN","FGF18")

#download and extract GSE47460 (GPL14550)
gse   <- getGEO("GSE47460", GSEMatrix = TRUE, getGPL = TRUE)
eset  <- gse[[which(sapply(gse, function(x) annotation(x) == "GPL14550"))]]
expr  <- exprs(eset)
pheno <- pData(eset)

#parse clinical metadata
extract_text <- function(pheno, pattern) {
  char_cols <- grep("^characteristics_ch1", colnames(pheno), value = TRUE)
  result <- rep(NA_character_, nrow(pheno))
  for (cc in char_cols) {
    vals <- pheno[[cc]]; matches <- grepl(pattern, vals, ignore.case = TRUE)
    if (any(matches)) result[matches] <- trimws(gsub(".*:\\s*", "", vals[matches]))
  }
  result
}
extract_num <- function(pheno, pattern) {
  char_cols <- grep("^characteristics_ch1", colnames(pheno), value = TRUE)
  result <- rep(NA_real_, nrow(pheno))
  for (cc in char_cols) {
    vals <- pheno[[cc]]; matches <- grepl(pattern, vals, ignore.case = TRUE)
    if (any(matches)) result[matches] <- suppressWarnings(as.numeric(gsub(".*:\\s*", "", vals[matches])))
  }
  result
}

pheno$disease_state <- extract_text(pheno, "^disease state")
pheno$ild_subtype   <- extract_text(pheno, "^ild subtype")
pheno$smoker        <- extract_text(pheno, "^smoker\\?")
pheno$sex           <- extract_text(pheno, "^Sex")
pheno$age           <- extract_num(pheno,  "^age")
pheno$fvc_pct       <- extract_num(pheno,  "%predicted fvc \\(pre-bd\\)")
pheno$dlco_pct      <- extract_num(pheno,  "%predicted dlco")

pheno$ild_subtype_clean <- case_when(
  grepl("UIP/IPF",  pheno$ild_subtype, ignore.case = TRUE) ~ "IPF/UIP",
  grepl("NSIP",     pheno$ild_subtype, ignore.case = TRUE) ~ "NSIP",
  grepl("HP",       pheno$ild_subtype, ignore.case = TRUE) ~ "HP",
  grepl("RB-ILD",   pheno$ild_subtype, ignore.case = TRUE) ~ "RB-ILD",
  grepl("Unclass",  pheno$ild_subtype, ignore.case = TRUE) ~ "Unclassifiable ILD",
  grepl("DIP",      pheno$ild_subtype, ignore.case = TRUE) ~ "DIP",
  !is.na(pheno$ild_subtype) & pheno$ild_subtype != ""      ~ pheno$ild_subtype,
  TRUE ~ NA_character_
)
pheno$diagnosis_clean <- case_when(
  pheno$disease_state == "Control"                          ~ "Control",
  !is.na(pheno$ild_subtype_clean)                          ~ pheno$ild_subtype_clean,
  pheno$disease_state == "Interstitial lung disease"        ~ "ILD (unspecified)",
  pheno$disease_state == "Chronic Obstructive Lung Disease" ~ "COPD",
  TRUE ~ NA_character_
)
pheno$diagnosis_broad <- case_when(
  pheno$diagnosis_clean == "Control" ~ "Control",
  pheno$diagnosis_clean == "IPF/UIP" ~ "IPF",
  pheno$diagnosis_clean == "COPD"    ~ "COPD",
  pheno$disease_state   == "Interstitial lung disease" ~ "Other ILD",
  TRUE ~ NA_character_
)

#probe to gene mapping (GPL14550 Agilent)
fdata           <- fData(eset)
probe_genes_raw <- as.character(fdata[["GENE_SYMBOL"]])
probe_genes     <- sapply(strsplit(probe_genes_raw, " /// "), function(x) trimws(x[1]))
names(probe_genes) <- rownames(fdata)

full_expr <- expr
rownames(full_expr) <- probe_genes[rownames(full_expr)]
full_expr <- full_expr[!is.na(rownames(full_expr)) & rownames(full_expr) != "", ]
full_expr_collapsed <- avereps(full_expr, ID = rownames(full_expr))

#run ssgsea
gene_sets <- list(
  Fibrosis_Signature = fibrosis_human[fibrosis_human %in% rownames(full_expr_collapsed)],
  AFOS_Signature     = afos_human[afos_human     %in% rownames(full_expr_collapsed)],
  Fibroblast_Score   = fibroblast_genes[fibroblast_genes %in% rownames(full_expr_collapsed)]
)

ssgsea_param  <- ssgseaParam(exprData = full_expr_collapsed, geneSets = gene_sets)
ssgsea_scores <- gsva(ssgsea_param)

pheno$fibrosis_ssgsea <- as.numeric(ssgsea_scores["Fibrosis_Signature", rownames(pheno)])
pheno$afos_ssgsea     <- as.numeric(ssgsea_scores["AFOS_Signature",     rownames(pheno)])
pheno$fib_score       <- as.numeric(ssgsea_scores["Fibroblast_Score",   rownames(pheno)])

sig_list <- list(
  list(label = "Fibrosis Signature", col = "fibrosis_ssgsea", prefix = "fibrosis"),
  list(label = "AFOS Signature",     col = "afos_ssgsea",     prefix = "afos")
)

#signature score boxplots, broad groups
broad_data <- pheno[!is.na(pheno$diagnosis_broad), ]
broad_data$diagnosis_broad <- factor(broad_data$diagnosis_broad,
                                     levels = c("Control", "IPF", "Other ILD", "COPD"))
broad_colors <- c("Control" = "#4DAF4A", "IPF" = "#E41A1C",
                  "Other ILD" = "#FF7F00", "COPD" = "#377EB8")

for (sig in sig_list) {
  stat_test <- broad_data %>%
    tukey_hsd(as.formula(paste(sig$col, "~ diagnosis_broad"))) %>%
    add_xy_position(x = "diagnosis_broad", step.increase = 0.06)
  aov_label <- paste0("One-way ANOVA p = ", format.pval(
    summary(aov(as.formula(paste(sig$col, "~ diagnosis_broad")), data = broad_data))[[1]][["Pr(>F)"]][1],
    digits = 3))
  
  print(
    ggplot(broad_data, aes(x = diagnosis_broad, y = .data[[sig$col]], fill = diagnosis_broad)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) +
      geom_jitter(width = 0.15, alpha = 0.25, size = 0.8) +
      scale_fill_manual(values = broad_colors) +
      stat_pvalue_manual(stat_test, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE) +
      labs(title = paste(sig$label, "— ssGSEA Score by Disease Group"),
           subtitle = aov_label, y = paste(sig$label, "ssGSEA"), x = "") +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(size = 12))
  )
}

#signature score boxplots, detailed groups
detailed_data <- pheno[!is.na(pheno$diagnosis_clean), ]
group_counts  <- table(detailed_data$diagnosis_clean)
detailed_data <- detailed_data[detailed_data$diagnosis_clean %in% names(group_counts[group_counts >= 3]), ]

group_medians <- detailed_data %>%
  group_by(diagnosis_clean) %>%
  summarise(med = median(fibrosis_ssgsea, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med))
detailed_data$diagnosis_clean <- factor(detailed_data$diagnosis_clean,
                                        levels = group_medians$diagnosis_clean)

n_groups   <- length(levels(detailed_data$diagnosis_clean))
detail_pal <- setNames(
  c("#E41A1C","#FF7F00","#984EA3","#4DAF4A","#377EB8",
    "#A65628","#F781BF","#999999","#66C2A5","#FC8D62")[1:n_groups],
  levels(detailed_data$diagnosis_clean)
)
if ("Control" %in% names(detail_pal)) detail_pal["Control"] <- "#4DAF4A"
if ("IPF/UIP" %in% names(detail_pal)) detail_pal["IPF/UIP"] <- "#E41A1C"
if ("COPD"    %in% names(detail_pal)) detail_pal["COPD"]    <- "#377EB8"

for (sig in sig_list) {
  stat_test_ctrl <- detailed_data %>%
    tukey_hsd(as.formula(paste(sig$col, "~ diagnosis_clean"))) %>%
    filter(group1 == "Control" | group2 == "Control") %>%
    add_xy_position(x = "diagnosis_clean", step.increase = 0.05)
  aov_label <- paste0("One-way ANOVA p = ", format.pval(
    summary(aov(as.formula(paste(sig$col, "~ diagnosis_clean")), data = detailed_data))[[1]][["Pr(>F)"]][1],
    digits = 3))
  
  print(
    ggplot(detailed_data, aes(x = diagnosis_clean, y = .data[[sig$col]], fill = diagnosis_clean)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
      scale_fill_manual(values = detail_pal) +
      stat_pvalue_manual(stat_test_ctrl, label = "p.adj.signif", tip.length = 0.01, hide.ns = FALSE) +
      labs(title = paste(sig$label, "ssGSEA — All ILD Subtypes"),
           subtitle = paste("vs Control |", aov_label),
           y = paste(sig$label, "ssGSEA"), x = "") +
      theme_minimal() +
      theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1, size = 10))
  )
}

#forest plots for control and IPF
forest_df <- pheno[pheno$diagnosis_broad %in% c("Control", "IPF"), ]
forest_df$sex_female   <- as.integer(forest_df$sex    == "2-Female")
forest_df$smoker_never <- as.integer(forest_df$smoker == "3-Never")
forest_df <- forest_df[complete.cases(forest_df[, c("age","sex_female","smoker_never","fib_score")]), ]
forest_df$age_z       <- as.numeric(scale(forest_df$age))
forest_df$fib_score_z <- as.numeric(scale(forest_df$fib_score))

for (sig in sig_list) {
  model_data <- forest_df[!is.na(forest_df[[sig$col]]), ]
  m          <- lm(as.formula(paste(sig$col, "~ age_z + sex_female + smoker_never + fib_score_z")),
                   data = model_data)
  coef_df           <- as.data.frame(summary(m)$coefficients)
  coef_df$Term      <- rownames(coef_df)
  colnames(coef_df) <- c("Estimate","SE","t","P","Term")
  ci                <- confint(m)
  coef_df$Lower     <- ci[, 1]
  coef_df$Upper     <- ci[, 2]
  coef_df$Label     <- coef_df$Term
  coef_df$Label     <- gsub("sex_female",      "Sex: Female",            coef_df$Label)
  coef_df$Label     <- gsub("smoker_never",    "Smoking: Never",         coef_df$Label)
  coef_df$Label     <- gsub("fib_score_z",     "Fibroblast ssGSEA (SD)", coef_df$Label)
  coef_df$Label     <- gsub("age_z",           "Age (SD)",               coef_df$Label)
  
  plot_df             <- coef_df[coef_df$Term != "(Intercept)", ]
  plot_df$Significant <- plot_df$P < 0.05
  plot_df             <- plot_df[order(plot_df$Estimate), ]
  plot_df$Label       <- factor(plot_df$Label, levels = plot_df$Label)
  
  print(
    ggplot(plot_df, aes(x = Estimate, y = Label, color = Significant)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
      geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2, linewidth = 0.8) +
      geom_point(size = 3.5) +
      scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#999999"),
                         labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
                         name = "Significance") +
      labs(title    = paste("Covariate Model:", sig$label),
           subtitle = paste0("Continuous predictors standardised to SD units | Control + IPF | n = ",
                             nrow(model_data)),
           x = "Coefficient (95% CI)", y = "") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 11))
  )
}

#unadjusted control and IPF FVC and DLCO scatter plots 
ipf_ctrl_colors <- c("Control" = "#4DAF4A", "IPF" = "#E41A1C")

for (sig in sig_list) {
  for (outcome in list(
    list(col = "fvc_pct",  label = "FVC % Predicted"),
    list(col = "dlco_pct", label = "DLCO % Predicted")
  )) {
    sub <- broad_data[broad_data$diagnosis_broad %in% c("Control","IPF") &
                        !is.na(broad_data[[outcome$col]]), ]
    sub$diagnosis_broad <- factor(sub$diagnosis_broad, levels = c("Control","IPF"))
    
    print(
      ggplot(sub, aes(x = .data[[sig$col]], y = .data[[outcome$col]])) +
        geom_point(aes(color = diagnosis_broad), alpha = 0.5, size = 1.8) +
        geom_smooth(method = "lm", color = "black", se = TRUE, linewidth = 0.8) +
        scale_color_manual(values = ipf_ctrl_colors) +
        stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",
                 size = 3.5, color = "black") +
        facet_wrap(~diagnosis_broad, scales = "free", nrow = 1) +
        labs(title    = paste(outcome$label, "vs", sig$label, "ssGSEA"),
             subtitle = "Unadjusted — Control and IPF only",
             x = paste(sig$label, "ssGSEA"), y = outcome$label) +
        theme_minimal() + theme(legend.position = "none")
    )
  }
}

#partial correlations (adjusted) control and IPF FVC and DLCO scatter plots
ipf_ctrl <- pheno[pheno$diagnosis_broad %in% c("Control","IPF"), ]
ipf_ctrl$diagnosis_broad <- factor(ipf_ctrl$diagnosis_broad, levels = c("Control","IPF"))
ipf_ctrl$smoker_binary <- case_when(
  ipf_ctrl$smoker == "3-Never"                              ~ "Never",
  ipf_ctrl$smoker %in% c("1-Current","2-Ever (>100)")      ~ "Ever/Current",
  TRUE ~ NA_character_
)
ipf_ctrl$smoker_binary <- factor(ipf_ctrl$smoker_binary, levels = c("Never","Ever/Current"))
ipf_ctrl_cc <- ipf_ctrl[complete.cases(ipf_ctrl[, c("age","sex","smoker_binary","fib_score")]), ]

for (sig in sig_list) {
  for (outcome in list(
    list(col = "fvc_pct",  label = "FVC % Predicted"),
    list(col = "dlco_pct", label = "DLCO % Predicted")
  )) {
    sub       <- ipf_ctrl_cc[!is.na(ipf_ctrl_cc[[outcome$col]]), ]
    resid_sig <- resid(lm(as.formula(paste(sig$col,       "~ age + sex + smoker_binary + fib_score")), data = sub))
    resid_out <- resid(lm(as.formula(paste(outcome$col,   "~ age + sex + smoker_binary + fib_score")), data = sub))
    ct        <- cor.test(resid_sig, resid_out, method = "spearman")
    
    plot_df           <- data.frame(resid_sig, resid_out, diagnosis = sub$diagnosis_broad)
    plot_df$diagnosis <- factor(plot_df$diagnosis, levels = c("Control","IPF"))
    subtitle          <- paste0("Overall partial rho = ", round(ct$estimate, 3),
                                ", p = ", format.pval(ct$p.value, digits = 3),
                                "\nResiduals after adjusting for age, sex, smoking, fibroblast score")
    
    print(
      ggplot(plot_df, aes(x = resid_sig, y = resid_out)) +
        geom_point(aes(color = diagnosis), alpha = 0.5, size = 1.8) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
        scale_color_manual(values = ipf_ctrl_colors) +
        stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",
                 size = 3.5, color = "black") +
        facet_wrap(~diagnosis, scales = "free", nrow = 1) +
        labs(title    = paste(sig$label, "vs", outcome$label, "— Adjusted"),
             subtitle = subtitle,
             x = paste(sig$label, "ssGSEA (residual)"),
             y = paste(outcome$label, "(residual)")) +
        theme_minimal() + theme(legend.position = "none")
    )
  }
}
