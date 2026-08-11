# Mouse lung multiome - part 3 of 4: published mouse scRNA-seq datasets
# Processes the published mouse lung fibroblast datasets used for comparison: Curras-Alonso et al
# (GSE211713), Tsukui et al (GSE132771), Zepp et al (GSE149563), Strunz et al (GSE141259) and Narvaez de Pilar et al
# (GSE180822). Each is clustered, subset to fibroblasts, and scored for the developmental
# and fibrotic ECM signatures.
#
# Inputs: a Seurat object saved as .rds for each published dataset 
# created from data deposited to GEO
# Outputs: fibroblasts_Curras, fibroblasts_Tsukui, fibroblasts_Zepp, fibroblasts_Strunz
#          and fibroblasts_Narvaez, each with ECM module scores added

# ── 0. Load Packages ──────────────────────────────────────────────────────────

library(Seurat)
library(harmony)
library(dplyr)

# ── 1. ECM signatures ─────────────────────────────────────────────────────────

#the same vectors defined in part 2, repeated here so this script runs on its own

#developmental ECM signature
all_mfb_ecm_genes <- c("Adam12", "Adam19", "Adamts10", "Adamts17", "Adamts6", "Adamts9", "Aspn", "Bmp2", "Cask", "Col14a1", "Col24a1", "Col25a1", "Col27a1", "Dst", "Eln", "Eng", "Fbln1", "Fbln5", "Fbn2", "Hmcn1", "Hpse2", "Htra1", "Itga1", "Itga9", "Itgav", "Itgb1", "Lrp12", "Ltbp2", "Mfap2", "Ntn4", "P3h2", "P4ha3", "Ppib", "Reck", "Sdc2", "Sparc", "Tgfb2", "Tgfbi")

#fibrotic ECM signature
matrix_genes_final <- c("Adamts12","Adamts15","Bgn","Bmp1","Col16a1","Col1a1","Col1a2","Col23a1","Col5a2","Col5a3","Col7a1","Ctsd","Fbn1","Fgf2","Fn1","Has2","Hspg2","Itga5","Lama5","Ndnf","Plec","Serpine1","Spp1","Thbs1","Timp1","Tnc")


# ── 2. GSE211713 scRNA-seq dataset analysis ───────────────────────────────────

#Curras-Alonso et al dataset
#initialize Seurat object from publicly available data from GSE211713
Curras_seurat <- readRDS(file = "/path/to/Curras_seurat") 

#add treatment labels
Ctrl <- c("SampleGSM6499593","SampleGSM6499594","SampleGSM6499595","SampleGSM6499596","SampleGSM6499597")
IR_1mo <- c("SampleGSM6499603","SampleGSM6499604")
IR_2mo <- c("SampleGSM6499605","SampleGSM6499606")
IR_3mo <- c("SampleGSM6499607","SampleGSM6499608")
IR_4mo <- c("SampleGSM6499609","SampleGSM6499610")
IR_5mo <- c("SampleGSM6499611","SampleGSM6499612")

Curras_seurat@meta.data[which(Curras_seurat@meta.data$orig.ident %in% Ctrl), "6group"] <- "Ctrl"
Curras_seurat@meta.data[which(Curras_seurat@meta.data$orig.ident %in% IR_1mo), "6group"] <- "IR_1mo"
Curras_seurat@meta.data[which(Curras_seurat@meta.data$orig.ident %in% IR_2mo), "6group"] <- "IR_2mo"
Curras_seurat@meta.data[which(Curras_seurat@meta.data$orig.ident %in% IR_3mo), "6group"] <- "IR_3mo"
Curras_seurat@meta.data[which(Curras_seurat@meta.data$orig.ident %in% IR_4mo), "6group"] <- "IR_4mo"
Curras_seurat@meta.data[which(Curras_seurat@meta.data$orig.ident %in% IR_5mo), "6group"] <- "IR_5mo"

Curras_seurat@meta.data[which(Curras_seurat@meta.data$orig.ident %in% Ctrl), "early_late"] <- "Ctrl"
Curras_seurat@meta.data[which(Curras_seurat@meta.data$orig.ident %in% IR_1mo), "early_late"] <- "Radiation_early"
Curras_seurat@meta.data[which(Curras_seurat@meta.data$orig.ident %in% IR_2mo), "early_late"] <- "Radiation_early"
Curras_seurat@meta.data[which(Curras_seurat@meta.data$orig.ident %in% IR_3mo), "early_late"] <- "Radiation_early"
Curras_seurat@meta.data[which(Curras_seurat@meta.data$orig.ident %in% IR_4mo), "early_late"] <- "Radiation_late"
Curras_seurat@meta.data[which(Curras_seurat@meta.data$orig.ident %in% IR_5mo), "early_late"] <- "Radiation_late"

#process RNA dataset
Curras_seurat_scaled <- SCTransform(Curras_seurat, vars.to.regress = c("nCount_RNA"), verbose = FALSE)
DefaultAssay(Curras_seurat_scaled) <- "SCT"
Curras_seurat_scaled <- RunPCA(Curras_seurat_scaled, npcs = 50)

#SCT assay integration
Curras_seurat_scaled <- RunHarmony(Curras_seurat_scaled, group.by.vars = c("orig.ident"), reduction = "pca",
                                   dims.use = NULL, theta = 2, lambda = NULL, sigma = 0.1,
                                   nclust = NULL, tau = 0, block.size = 0.05, max.iter.harmony = 10,
                                   max.iter.cluster = 20, epsilon.cluster = 1e-05,
                                   epsilon.harmony = 1e-04, plot_convergence = FALSE, verbose = TRUE,
                                   reference_values = NULL, reduction.save = "harmonyBatch",
                                   assay.use = "SCT", project.dim = TRUE)

#clustering
Curras_seurat_scaled <- RunUMAP(Curras_seurat_scaled, reduction = "harmonyBatch", dims = 1:50)
Curras_seurat_scaled <- FindNeighbors(Curras_seurat_scaled, reduction = "harmonyBatch", dims = 1:50, k.param = 10)
Curras_seurat_scaled <- FindClusters(Curras_seurat_scaled, resolution = 0.6)

#subset the Col1a1+ clusters
Idents(Curras_seurat_scaled) <- "seurat_clusters"
fibroblasts_Curras <- subset(x = Curras_seurat_scaled, idents = c("6","23","38"))

#process RNA dataset
fibroblasts_Curras <- SCTransform(fibroblasts_Curras, verbose = FALSE)
DefaultAssay(fibroblasts_Curras) <- "SCT"
fibroblasts_Curras <- RunPCA(fibroblasts_Curras, npcs = 50)

#process RNA assay
DefaultAssay(fibroblasts_Curras) <- "RNA"
fibroblasts_Curras <- NormalizeData(fibroblasts_Curras, scale.factor = 10000)
fibroblasts_Curras <- FindVariableFeatures(fibroblasts_Curras, selection.method = "mean.var.plot", nfeatures = 2000)
fibroblasts_Curras <- ScaleData(fibroblasts_Curras, features = rownames(fibroblasts_Curras))

#identify differentially expressed genes between fibroblasts in control and 4 month irradiated mice
Idents(fibroblasts_Curras) <- "6group"
DefaultAssay(fibroblasts_Curras) <- "RNA"
IR_4mo_v_ctrl <- FindMarkers(fibroblasts_Curras, ident.1 = "IR_4mo", ident.2 = "Ctrl")

#add ECM signatures
DefaultAssay(fibroblasts_Curras) <- "RNA"
fibroblasts_Curras <- AddModuleScore(fibroblasts_Curras, features = list(all_mfb_ecm_genes), name = "all_mfb_ecm_genes")
fibroblasts_Curras <- AddModuleScore(fibroblasts_Curras, features = list(matrix_genes_final), name = "matrix_genes_final")

#pseudobulk fibroblasts_Curras
#pseudobulking gene signature module scores
#choose pseudobulking variables
group_var     <- "orig.ident"
min_cells     <- 30
module_scores <- c(
  "all_mfb_ecm_genes1",
  "matrix_genes_final1"
)

DefaultAssay(fibroblasts_Curras) <- "RNA"

#create cell counts per group and filter
group_counts <- table(fibroblasts_Curras@meta.data[[group_var]])
keep_groups  <- names(group_counts[group_counts >= min_cells])

#subset object
obj_filtered <- fibroblasts_Curras[, fibroblasts_Curras@meta.data[[group_var]] %in% keep_groups]

#pseudobulk module scores
group_meta <- obj_filtered@meta.data %>%
  select(all_of(c(group_var, module_scores))) %>%
  group_by(.data[[group_var]]) %>%
  summarise(
    n_cells = n(),
    across(all_of(module_scores), mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(.data[[group_var]])


# ── 3. GSE132771 scRNA-seq dataset analysis ───────────────────────────────────

#Tsukui et al dataset
#initialize Seurat object from publicly available data from GSE132771
Tsukui_seurat <- readRDS(file = "/path/to/Tsukui_seurat") 

#combine the Col1a1+ and Col1a1- from each mouse together
Idents(Tsukui_seurat) <- "orig.ident"
UT1 <- c("SampleGSM3891616","SampleGSM3891618")
UT2 <- c("SampleGSM3891617","SampleGSM3891619")
Bleo1 <- c("SampleGSM3891612","SampleGSM3891614")
Bleo2 <- c("SampleGSM3891613","SampleGSM3891615")
Tsukui_seurat@meta.data[which(Tsukui_seurat@meta.data$orig.ident %in% UT1), "sampleID"] <- "UT1"
Tsukui_seurat@meta.data[which(Tsukui_seurat@meta.data$orig.ident %in% UT2), "sampleID"] <- "UT2"
Tsukui_seurat@meta.data[which(Tsukui_seurat@meta.data$orig.ident %in% Bleo1), "sampleID"] <- "Bleo1"
Tsukui_seurat@meta.data[which(Tsukui_seurat@meta.data$orig.ident %in% Bleo2), "sampleID"] <- "Bleo2"

#process RNA data
Tsukui_seurat_scaled <- SCTransform(Tsukui_seurat, vars.to.regress = c("nCount_RNA"), verbose = FALSE)
DefaultAssay(Tsukui_seurat_scaled) <- "SCT"
Tsukui_seurat_scaled <- RunPCA(Tsukui_seurat_scaled, npcs = 50)

#SCT assay integration
Tsukui_seurat_scaled <- RunHarmony(Tsukui_seurat_scaled, group.by.vars = c("sampleID"), reduction = "pca",
                                   dims.use = NULL, theta = 2, lambda = NULL, sigma = 0.1,
                                   nclust = NULL, tau = 0, block.size = 0.05, max.iter.harmony = 10,
                                   max.iter.cluster = 20, epsilon.cluster = 1e-05,
                                   epsilon.harmony = 1e-04, plot_convergence = FALSE, verbose = TRUE,
                                   reference_values = NULL, reduction.save = "harmonyBatch",
                                   assay.use = "SCT", project.dim = TRUE)

#clustering
Tsukui_seurat_scaled <- RunUMAP(Tsukui_seurat_scaled, reduction = "harmonyBatch", dims = 1:50)
Tsukui_seurat_scaled <- FindNeighbors(Tsukui_seurat_scaled, reduction = "harmonyBatch", dims = 1:50, k.param = 10)
Tsukui_seurat_scaled <- FindClusters(Tsukui_seurat_scaled, resolution = 0.6)

#subset the Col1a1+ clusters
Idents(Tsukui_seurat_scaled) <- "seurat_clusters"
fibroblasts_Tsukui <- subset(x = Tsukui_seurat_scaled, idents = c("0","2","3","5","6","9","16","31"))

#process RNA data
fibroblasts_Tsukui <- SCTransform(fibroblasts_Tsukui, verbose = FALSE)
DefaultAssay(fibroblasts_Tsukui) <- "SCT"
fibroblasts_Tsukui <- RunPCA(fibroblasts_Tsukui, npcs = 50)

#process RNA assay
DefaultAssay(fibroblasts_Tsukui) <- "RNA"
fibroblasts_Tsukui <- NormalizeData(fibroblasts_Tsukui, scale.factor = 10000)
fibroblasts_Tsukui <- FindVariableFeatures(fibroblasts_Tsukui, selection.method = "mean.var.plot", nfeatures = 2000)
fibroblasts_Tsukui <- ScaleData(fibroblasts_Tsukui, features = rownames(fibroblasts_Tsukui))

#add general untreated and bleomycin labels
Idents(fibroblasts_Tsukui) <- "sampleID"
control <- c("UT1","UT2")
bleo <- c("Bleo1","Bleo2")
fibroblasts_Tsukui@meta.data[which(fibroblasts_Tsukui@meta.data$sampleID %in% control), "condition"] <- "control"
fibroblasts_Tsukui@meta.data[which(fibroblasts_Tsukui@meta.data$sampleID %in% bleo), "condition"] <- "bleo"

#identify differentially expressed genes between fibroblasts in untreated and bleomycin treated mice
Idents(fibroblasts_Tsukui) <- "condition"
DefaultAssay(fibroblasts_Tsukui) <- "RNA"
bleo_vs_ctrl <- FindMarkers(fibroblasts_Tsukui, ident.1 = "bleo", ident.2 = "control")

#add ECM signatures
DefaultAssay(fibroblasts_Tsukui) <- "RNA"
fibroblasts_Tsukui <- AddModuleScore(fibroblasts_Tsukui, features = list(all_mfb_ecm_genes), name = "all_mfb_ecm_genes")
fibroblasts_Tsukui <- AddModuleScore(fibroblasts_Tsukui, features = list(matrix_genes_final), name = "matrix_genes_final")

#pseudobulk fibroblasts_Tsukui
#pseudobulking gene signature module scores
#choose pseudobulking variables
group_var     <- "sampleID"
min_cells     <- 30
module_scores <- c(
  "all_mfb_ecm_genes1",
  "matrix_genes_final1"
)

DefaultAssay(fibroblasts_Tsukui) <- "RNA"

#create cell counts per group and filter
group_counts <- table(fibroblasts_Tsukui@meta.data[[group_var]])
keep_groups  <- names(group_counts[group_counts >= min_cells])

#subset object
obj_filtered <- fibroblasts_Tsukui[, fibroblasts_Tsukui@meta.data[[group_var]] %in% keep_groups]

#pseudobulk module scores
group_meta <- obj_filtered@meta.data %>%
  select(all_of(c(group_var, module_scores))) %>%
  group_by(.data[[group_var]]) %>%
  summarise(
    n_cells = n(),
    across(all_of(module_scores), mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(.data[[group_var]])


# ── 4. GSE149563 scRNA-seq dataset analysis ───────────────────────────────────

#Zepp et al dataset
#initialize Seurat object from publicly available data from GSE149563
Zepp_seurat <- readRDS(file = "/path/to/Zepp_seurat") 

#process RNA data
Zepp_seurat <- subset(Zepp_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & nCount_RNA > 200 & nCount_RNA < 75000 & percent.mt < 15)
Zepp_seurat <- SCTransform(Zepp_seurat, vars.to.regress = c("nCount_RNA"), verbose = FALSE)
DefaultAssay(Zepp_seurat) <- "SCT"
Zepp_seurat <- RunPCA(Zepp_seurat, npcs = 50)

#SCT assay integration
Zepp_seurat_scaled <- RunHarmony(Zepp_seurat, group.by.vars = c("orig.ident"), reduction = "pca",
                                 dims.use = NULL, theta = 2, lambda = NULL, sigma = 0.1,
                                 nclust = NULL, tau = 0, block.size = 0.05, max.iter.harmony = 10,
                                 max.iter.cluster = 20, epsilon.cluster = 1e-05,
                                 epsilon.harmony = 1e-04, plot_convergence = FALSE, verbose = TRUE,
                                 reference_values = NULL, reduction.save = "harmonyBatch",
                                 assay.use = "SCT", project.dim = TRUE)

#clustering
Zepp_seurat_scaled <- RunUMAP(Zepp_seurat_scaled, reduction = "harmonyBatch", dims = 1:40)
Zepp_seurat_scaled <- FindNeighbors(Zepp_seurat_scaled, reduction = "harmonyBatch", dims = 1:40, k.param = 10)
Zepp_seurat_scaled <- FindClusters(Zepp_seurat_scaled, resolution = 0.6)

#subset the Col1a1+ clusters
Idents(Zepp_seurat_scaled) <- "seurat_clusters"
fibroblasts <- subset(x = Zepp_seurat_scaled, idents = c("1","2","3","8","11","13","15","18","19","23","24","29","30","32","34"))

#process RNA data
fibroblasts <- SCTransform(fibroblasts, verbose = FALSE)
DefaultAssay(fibroblasts) <- "SCT"
fibroblasts <- RunPCA(fibroblasts, npcs = 50)

#clustering
fibroblasts <- RunUMAP(fibroblasts, reduction = "pca", dims = 1:40)
fibroblasts <- FindNeighbors(fibroblasts, reduction = "pca", dims = 1:40, k.param = 10)
fibroblasts <- FindClusters(fibroblasts, resolution = 0.6)

#re-subset the Col1a1+ clusters
Idents(fibroblasts) <- "seurat_clusters"
fibroblasts_Zepp <- subset(x = fibroblasts, idents = c("1","2","3","4","5","6","7","8","9","10","11","12","13","16","17","20","21","22","23","24","25","26","29","31","32","33","34","36","37"))

#process RNA data
fibroblasts_Zepp <- SCTransform(fibroblasts_Zepp, verbose = FALSE)
DefaultAssay(fibroblasts_Zepp) <- "SCT"
fibroblasts_Zepp <- RunPCA(fibroblasts_Zepp, npcs = 50)

#process RNA assay
DefaultAssay(fibroblasts_Zepp) <- "RNA"
fibroblasts_Zepp <- NormalizeData(fibroblasts_Zepp, scale.factor = 10000)
fibroblasts_Zepp <- FindVariableFeatures(fibroblasts_Zepp, selection.method = "mean.var.plot", nfeatures = 2000)
fibroblasts_Zepp <- ScaleData(fibroblasts_Zepp, features = rownames(fibroblasts_Zepp))

#add ECM signatures
DefaultAssay(fibroblasts_Zepp) <- "RNA"
fibroblasts_Zepp <- AddModuleScore(fibroblasts_Zepp, features = list(all_mfb_ecm_genes), name = "all_mfb_ecm_genes")
fibroblasts_Zepp <- AddModuleScore(fibroblasts_Zepp, features = list(matrix_genes_final), name = "matrix_genes_final")


# ── 5. GSE141259 scRNA-seq dataset analysis ───────────────────────────────────

#Strunz et al dataset
#initialize Seurat object from publicly available data from GSE141259
Strunz_seurat <- readRDS(file = "/path/to/Strunz_seurat") 

#process RNA data
Strunz_RDS_scaled <- SCTransform(Strunz_seurat, vars.to.regress = c("nCount_RNA"), verbose = FALSE)
DefaultAssay(Strunz_RDS_scaled) <- "SCT"
Strunz_RDS_scaled <- RunPCA(Strunz_RDS_scaled, npcs = 50)

#SCT assay integration
Strunz_RDS_scaled <- RunHarmony(Strunz_RDS_scaled, group.by.vars = c("identifier"), reduction = "pca",
                                dims.use = NULL, theta = 2, lambda = NULL, sigma = 0.1,
                                nclust = NULL, tau = 0, block.size = 0.05, max.iter.harmony = 10,
                                max.iter.cluster = 20, epsilon.cluster = 1e-05,
                                epsilon.harmony = 1e-04, plot_convergence = FALSE, verbose = TRUE,
                                reference_values = NULL, reduction.save = "harmonyBatch",
                                assay.use = "SCT", project.dim = TRUE)

#clustering
Strunz_RDS_scaled <- RunUMAP(Strunz_RDS_scaled, reduction = "harmonyBatch", dims = 1:50)
Strunz_RDS_scaled <- FindNeighbors(Strunz_RDS_scaled, reduction = "harmonyBatch", dims = 1:50, k.param = 20)
Strunz_RDS_scaled <- FindClusters(Strunz_RDS_scaled, resolution = 0.7)

#subset the Col1a1+ clusters
Idents(Strunz_RDS_scaled) <- "seurat_clusters"
fibroblasts_Strunz <- subset(x = Strunz_RDS_scaled, idents = c("9","13"))

#process RNA data
fibroblasts_Strunz <- SCTransform(fibroblasts_Strunz, verbose = FALSE)
DefaultAssay(fibroblasts_Strunz) <- "SCT"
fibroblasts_Strunz <- RunPCA(fibroblasts_Strunz, npcs = 50)

#process RNA assay
DefaultAssay(fibroblasts_Strunz) <- "RNA"
fibroblasts_Strunz <- NormalizeData(fibroblasts_Strunz, scale.factor = 10000)
fibroblasts_Strunz <- FindVariableFeatures(fibroblasts_Strunz, selection.method = "mean.var.plot", nfeatures = 2000)
fibroblasts_Strunz <- ScaleData(fibroblasts_Strunz, features = rownames(fibroblasts_Strunz))

#add ECM signatures
DefaultAssay(fibroblasts_Strunz) <- "RNA"
fibroblasts_Strunz <- AddModuleScore(fibroblasts_Strunz, features = list(all_mfb_ecm_genes), name = "all_mfb_ecm_genes")
fibroblasts_Strunz <- AddModuleScore(fibroblasts_Strunz, features = list(matrix_genes_final), name = "matrix_genes_final")

#pseudobulk fibroblasts_Strunz
#pseudobulking gene signature module scores
#choose pseudobulking variables
group_var     <- "identifier"
min_cells     <- 30
module_scores <- c(
  "all_mfb_ecm_genes1",
  "matrix_genes_final1"
)

DefaultAssay(fibroblasts_Strunz) <- "RNA"

#create cell counts per group and filter
group_counts <- table(fibroblasts_Strunz@meta.data[[group_var]])
keep_groups  <- names(group_counts[group_counts >= min_cells])

#subset object
obj_filtered <- fibroblasts_Strunz[, fibroblasts_Strunz@meta.data[[group_var]] %in% keep_groups]

#pseudobulk module scores
group_meta <- obj_filtered@meta.data %>%
  select(all_of(c(group_var, module_scores))) %>%
  group_by(.data[[group_var]]) %>%
  summarise(
    n_cells = n(),
    across(all_of(module_scores), mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(.data[[group_var]])


# ── 6. GSE180822 scRNA-seq dataset analysis ───────────────────────────────────

#Narvaez del Pilar et al dataset
#initialize Seurat object from publicly available data from GSE180822
Narvaez_seurat <- readRDS(file = "/path/to/Narvaez_seurat") 

#add timepoint meta data
Idents(Narvaez_seurat) <- "orig.ident"
E17 <- "SampleGSM5471469"
E19 <- "SampleGSM5471470"
P7 <- "SampleGSM5471471"
P13 <- "SampleGSM5471472"
P20 <- "SampleGSM5471473"
P70 <- "SampleGSM5471474"

Narvaez_seurat@meta.data[which(Narvaez_seurat@meta.data$orig.ident %in% E17), "timepoint"] <- "E17"
Narvaez_seurat@meta.data[which(Narvaez_seurat@meta.data$orig.ident %in% E19), "timepoint"] <- "E19"
Narvaez_seurat@meta.data[which(Narvaez_seurat@meta.data$orig.ident %in% P7), "timepoint"] <- "P7"
Narvaez_seurat@meta.data[which(Narvaez_seurat@meta.data$orig.ident %in% P13), "timepoint"] <- "P13"
Narvaez_seurat@meta.data[which(Narvaez_seurat@meta.data$orig.ident %in% P20), "timepoint"] <- "P20"
Narvaez_seurat@meta.data[which(Narvaez_seurat@meta.data$orig.ident %in% P70), "timepoint"] <- "P70"

#process RNA data
Narvaez_seurat_scaled <- SCTransform(Narvaez_seurat, vars.to.regress = c("nCount_RNA"), verbose = FALSE)
DefaultAssay(Narvaez_seurat_scaled) <- "SCT"
Narvaez_seurat_scaled <- RunPCA(Narvaez_seurat_scaled, npcs = 50)

#SCT assay integration
Narvaez_seurat_scaled <- RunHarmony(Narvaez_seurat_scaled, group.by.vars = c("timepoint"), reduction = "pca",
                             dims.use = NULL, theta = 2, lambda = NULL, sigma = 0.1,
                             nclust = NULL, tau = 0, block.size = 0.05, max.iter.harmony = 10,
                             max.iter.cluster = 20, epsilon.cluster = 1e-05,
                             epsilon.harmony = 1e-04, plot_convergence = FALSE, verbose = TRUE,
                             reference_values = NULL, reduction.save = "harmonyBatch",
                             assay.use = "SCT", project.dim = TRUE)

#clustering
Narvaez_seurat_scaled <- RunUMAP(Narvaez_seurat_scaled, reduction = "harmonyBatch", dims = 1:40)
Narvaez_seurat_scaled <- FindNeighbors(Narvaez_seurat_scaled, reduction = "harmonyBatch", dims = 1:40, k.param = 10)
Narvaez_seurat_scaled <- FindClusters(Narvaez_seurat_scaled, resolution = 0.6)

#subset the Col1a1+ clusters
Idents(Narvaez_seurat_scaled) <- "seurat_clusters"
fibroblasts <- subset(x = Narvaez_seurat_scaled, idents = c("1","2","11","15","19","20","26","32","35"))

#process RNA data
fibroblasts <- SCTransform(fibroblasts, verbose = FALSE)
DefaultAssay(fibroblasts) <- "SCT"
fibroblasts <- RunPCA(fibroblasts, npcs = 50)

#clustering
fibroblasts <- RunUMAP(fibroblasts, reduction = "pca", dims = 1:50)
fibroblasts <- FindNeighbors(fibroblasts, reduction = "pca", dims = 1:50, k.param = 30)
fibroblasts <- FindClusters(fibroblasts, resolution = 0.6)

#re-subset Col1a1+ cells
Idents(fibroblasts) <- "seurat_clusters"
fibroblasts_Narvaez <- subset(x = fibroblasts, idents = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","17","18"))

#process RNA data
fibroblasts_Narvaez <- SCTransform(fibroblasts_Narvaez, verbose = FALSE)
DefaultAssay(fibroblasts_Narvaez) <- "SCT"
fibroblasts_Narvaez <- RunPCA(fibroblasts_Narvaez, npcs = 50)

#process RNA assay
DefaultAssay(fibroblasts_Narvaez) <- "RNA"
fibroblasts_Narvaez <- NormalizeData(fibroblasts_Narvaez, scale.factor = 10000)
fibroblasts_Narvaez <- FindVariableFeatures(fibroblasts_Narvaez, selection.method = "mean.var.plot", nfeatures = 2000)
fibroblasts_Narvaez <- ScaleData(fibroblasts_Narvaez, features = rownames(fibroblasts_Narvaez))

#add ECM signatures
DefaultAssay(fibroblasts_Narvaez) <- "RNA"
fibroblasts_Narvaez <- AddModuleScore(fibroblasts_Narvaez, features = list(all_mfb_ecm_genes), name = "all_mfb_ecm_genes")
fibroblasts_Narvaez <- AddModuleScore(fibroblasts_Narvaez, features = list(matrix_genes_final), name = "matrix_genes_final")


# ── 7. Save objects ───────────────────────────────────────────────────────────

saveRDS(fibroblasts_Curras,  "/path/to/fibroblasts_Curras.rds")
saveRDS(fibroblasts_Tsukui,  "/path/to/fibroblasts_Tsukui.rds")
saveRDS(fibroblasts_Zepp,    "/path/to/fibroblasts_Zepp.rds")
saveRDS(fibroblasts_Strunz,  "/path/to/fibroblasts_Strunz.rds")
saveRDS(fibroblasts_Narvaez, "/path/to/fibroblasts_Narvaez.rds")
