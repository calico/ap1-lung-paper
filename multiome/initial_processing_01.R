# Mouse lung multiome - part 1 of 4: initial processing of all cells and subsetting to fibroblasts
# Builds the multiome objects from the Cell Ranger ARC output: QC, RNA and ATAC processing,
# Harmony integration, WNN clustering, and the iterative mesenchymal and fibroblast
# subsetting that produces fibroblasts_final.
#
# Inputs: filtered_feature_bc_matrix_<n>.h5 and atac_fragments_<n>.tsv.gz for samples 1-20
# Outputs: multiome_integrated and fibroblasts_final

# ── 0. Load Packages ──────────────────────────────────────────────────────────

library(Seurat)
library(Signac)
library(harmony)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(homologene)
library(ggplot2)
library(dplyr)

# ── 1. Create Seurat Object ───────────────────────────────────────────────────

#define samples and base directory
samples <- 1:20
directory <- "/path/to/files/"
seurat_list <- list()

#get gene annotations for mm10
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"

for (sample in samples) {
  #define file paths and object name
  h5_file <- file.path(directory, paste0("filtered_feature_bc_matrix_", sample, ".h5"))
  atac_file <- file.path(directory, paste0("atac_fragments_", sample, ".tsv.gz"))
  
  #load the RNA and ATAC data
  counts <- Read10X_h5(h5_file)
  fragpath <- atac_file
  
  #create a Seurat object containing the RNA data
  seurat_obj <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA"
  )
  
  #create ATAC assay and add it to the object
  seurat_obj[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    genome = "mm10",
    annotation = annotation
  ) 
  
  #store in list
  seurat_list[[paste0("Sample_", sample)]] <- seurat_obj
}

#merge Seurat objects
multiome_seurat <- merge(seurat_list[[1]], y = seurat_list[-1], add.cell.ids = samples)

#add sample ID metadata
multiome_seurat$sampleID <- sapply(strsplit(Cells(multiome_seurat), "_"), `[`, 1)
multiome_seurat$sampleID <- factor(multiome_seurat$sampleID)


# ── 2. QC filtering ───────────────────────────────────────────────────────────

#add percent mitochondrial DNA, nucleosome signal, and transcription start site enrichment
multiome_seurat[["percent.mt"]] <- PercentageFeatureSet(multiome_seurat, pattern = "^mt-")
multiome_seurat <- NucleosomeSignal(multiome_seurat)
multiome_seurat <- TSSEnrichment(multiome_seurat)

#filter cells
multiome_seurat_filtered <- subset(
  x = multiome_seurat,
  subset = nCount_RNA > 200 & nCount_RNA < 25000 &
    nFeature_RNA > 200 & nFeature_RNA < 2500 &
    percent.mt < 10 &
    nCount_ATAC > 700 & nCount_ATAC < 100000 &
    nucleosome_signal < 2 & TSS.enrichment > 2
)


# ── 3. Process RNA and ATAC in whole dataset ──────────────────────────────────

#process RNA assay
DefaultAssay(multiome_seurat_filtered) <- "RNA"
multiome_filtered_SCT <- SCTransform(multiome_seurat_filtered, vars.to.regress = c("percent.mt","nCount_RNA"))
multiome_filtered_SCT <- RunPCA(multiome_filtered_SCT)
ElbowPlot(multiome_filtered_SCT)

#process ATAC assay
DefaultAssay(multiome_filtered_SCT) <- "ATAC"
multiome_filtered_SCT_ATAC <- FindTopFeatures(multiome_filtered_SCT, min.cutoff = "q5")
multiome_filtered_SCT_ATAC <- RunTFIDF(multiome_filtered_SCT_ATAC)
multiome_filtered_SCT_ATAC <- RunSVD(multiome_filtered_SCT_ATAC)

#run harmony for batch integration
#SCT assay integration
DefaultAssay(multiome_filtered_SCT_ATAC) <- "SCT"
multiome_filtered_SCT_ATAC <- RunHarmony(multiome_filtered_SCT_ATAC, group.by.vars = "sampleID", 
                                         reduction = "pca", theta = 2, tau = 0, max.iter.harmony = 10, 
                                         max.iter.cluster = 20, epsilon.cluster = 1e-05, 
                                         epsilon.harmony = 1e-04, plot_convergence = TRUE, verbose = TRUE, 
                                         reduction.save = "harmonySCT", assay.use = "SCT", project.dim = FALSE)

#ATAC assay integration
DefaultAssay(multiome_filtered_SCT_ATAC) <- "ATAC"
multiome_filtered_SCT_ATAC <- RunHarmony(multiome_filtered_SCT_ATAC, group.by.vars = "sampleID", 
                                         reduction = "lsi", theta = 2, tau = 0, max.iter.harmony = 10, 
                                         max.iter.cluster = 20, epsilon.cluster = 1e-05, 
                                         epsilon.harmony = 1e-04, plot_convergence = TRUE, verbose = TRUE, 
                                         reduction.save = "harmonyATAC", assay.use = "ATAC", project.dim = FALSE)

multiome_integrated <- multiome_filtered_SCT_ATAC

#WNN analysis
multiome_integrated <- FindMultiModalNeighbors(multiome_integrated, reduction.list = list("harmonySCT", "harmonyATAC"), 
                                               dims.list = list(1:40, 2:40))
multiome_integrated <- RunUMAP(multiome_integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                               reduction.key = "wnnUMAP_")
multiome_integrated <- FindClusters(multiome_integrated, graph.name = "wsnn", algorithm = 3, resolution = 0.4)

#add sample metadata
#condition
P7_UT <- c("1","3","7","11")
adult_saline <- c("2","4","8","12")
adult_bleo <- c("5","9","13","15")
adult_UT <- c("6","10","14","16")
aged_UT <- c("17","18","19","20")

#broad condition
#combines adult_saline and adult_UT samples into adult_healthy
P7 <- c("1","3","7","11")
adult_healthy <- c("2","4","8","12","6","10","14","16")
adult_bleo <- c("5","9","13","15")
aged <- c("17","18","19","20")

multiome_integrated@meta.data[which(multiome_integrated@meta.data$sampleID %in% P7_UT), "condition"] <- "P7_UT"
multiome_integrated@meta.data[which(multiome_integrated@meta.data$sampleID %in% adult_saline), "condition"] <- "adult_saline"
multiome_integrated@meta.data[which(multiome_integrated@meta.data$sampleID %in% adult_bleo), "condition"] <- "adult_bleo"
multiome_integrated@meta.data[which(multiome_integrated@meta.data$sampleID %in% adult_UT), "condition"] <- "adult_UT"
multiome_integrated@meta.data[which(multiome_integrated@meta.data$sampleID %in% aged_UT), "condition"] <- "aged_UT"

multiome_integrated@meta.data[which(multiome_integrated@meta.data$sampleID %in% P7), "condition_broad"] <- "P7"
multiome_integrated@meta.data[which(multiome_integrated@meta.data$sampleID %in% adult_healthy), "condition_broad"] <- "adult_healthy"
multiome_integrated@meta.data[which(multiome_integrated@meta.data$sampleID %in% adult_bleo), "condition_broad"] <- "adult_bleo"
multiome_integrated@meta.data[which(multiome_integrated@meta.data$sampleID %in% aged), "condition_broad"] <- "aged"


# ── 4. Subset mesenchymal cells ───────────────────────────────────────────────

#strategy for subsetting mesenchymal cells: be liberal with including clusters,
#then re-cluster and identify the true mesenchymal clusters. iterate until
#confident in mesenchymal cell identity

#mesenchymal subset 1 & processing
Idents(multiome_integrated) <- "seurat_clusters"
mesenchymal_cells <- subset(multiome_integrated, idents = c("0","8","12","14","15","16","17","24","33","37","38","43","45"))
#RNA processing
DefaultAssay(mesenchymal_cells) <- "RNA"
mesenchymal_cells_SCT <- SCTransform(mesenchymal_cells, vars.to.regress = c("percent.mt","nCount_RNA"))
mesenchymal_cells_SCT <- RunPCA(mesenchymal_cells_SCT)
ElbowPlot(mesenchymal_cells_SCT)
#ATAC processing
DefaultAssay(mesenchymal_cells_SCT) <- "ATAC"
mesenchymal_cells_SCT_ATAC <- FindTopFeatures(mesenchymal_cells_SCT, min.cutoff = "q5")
mesenchymal_cells_SCT_ATAC <- RunTFIDF(mesenchymal_cells_SCT_ATAC)
mesenchymal_cells_SCT_ATAC <- RunSVD(mesenchymal_cells_SCT_ATAC)
#SCT assay integration
DefaultAssay(mesenchymal_cells_SCT_ATAC) <- "SCT"
mesenchymal_cells_SCT_ATAC <- RunHarmony(mesenchymal_cells_SCT_ATAC, group.by.vars = "sampleID", 
                                         reduction = "pca", theta = 2, tau = 0, max.iter.harmony = 10, 
                                         max.iter.cluster = 20, epsilon.cluster = 1e-05, 
                                         epsilon.harmony = 1e-04, plot_convergence = TRUE, verbose = TRUE, 
                                         reduction.save = "harmonySCT", assay.use = "SCT", project.dim = FALSE)
#ATAC assay integration
DefaultAssay(mesenchymal_cells_SCT_ATAC) <- "ATAC"
mesenchymal_cells_SCT_ATAC <- RunHarmony(mesenchymal_cells_SCT_ATAC, group.by.vars = "sampleID", 
                                         reduction = "lsi", theta = 2, tau = 0, max.iter.harmony = 10, 
                                         max.iter.cluster = 20, epsilon.cluster = 1e-05, 
                                         epsilon.harmony = 1e-04, plot_convergence = TRUE, verbose = TRUE, 
                                         reduction.save = "harmonyATAC", assay.use = "ATAC", project.dim = FALSE)

mesenchymal_integrated <- mesenchymal_cells_SCT_ATAC
  
#WNN analysis
mesenchymal_integrated <- FindMultiModalNeighbors(mesenchymal_integrated, reduction.list = list("harmonySCT", "harmonyATAC"), 
                                                  dims.list = list(1:40, 2:40), k.nn = 50)
mesenchymal_integrated <- RunUMAP(mesenchymal_integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                                  reduction.key = "wnnUMAP_")
mesenchymal_integrated <- FindClusters(mesenchymal_integrated, graph.name = "wsnn", algorithm = 3, resolution = 0.5)
DimPlot(mesenchymal_integrated, reduction = "wnn.umap", label = TRUE, repel = TRUE, label.size = 2.5, raster = FALSE) + NoLegend()

#mesenchymal subset 2 & processing
Idents(mesenchymal_integrated) <- "seurat_clusters"
mesenchymal_cells2 <- subset(mesenchymal_integrated, idents = c("0","1","2","3","6","7","8","9","16","17","18","19"))
#RNA processing
DefaultAssay(mesenchymal_cells2) <- "RNA"
mesenchymal_cells2_SCT <- SCTransform(mesenchymal_cells2, vars.to.regress = c("percent.mt","nCount_RNA"))
mesenchymal_cells2_SCT <- RunPCA(mesenchymal_cells2_SCT)
ElbowPlot(mesenchymal_cells2_SCT)
#ATAC processing
DefaultAssay(mesenchymal_cells2_SCT) <- "ATAC"
multiome_filtered2_SCT_ATAC <- FindTopFeatures(mesenchymal_cells2_SCT, min.cutoff = "q5")
multiome_filtered2_SCT_ATAC <- RunTFIDF(multiome_filtered2_SCT_ATAC)
multiome_filtered2_SCT_ATAC <- RunSVD(multiome_filtered2_SCT_ATAC)
#SCT assay integration
DefaultAssay(multiome_filtered2_SCT_ATAC) <- "SCT"
multiome_filtered2_SCT_ATAC <- RunHarmony(multiome_filtered2_SCT_ATAC, group.by.vars = "sampleID", 
                                          reduction = "pca", theta = 2, tau = 0, max.iter.harmony = 10, 
                                          max.iter.cluster = 20, epsilon.cluster = 1e-05, 
                                          epsilon.harmony = 1e-04, plot_convergence = TRUE, verbose = TRUE, 
                                          reduction.save = "harmonySCT", assay.use = "SCT", project.dim = FALSE)
#ATAC assay integration
DefaultAssay(multiome_filtered2_SCT_ATAC) <- "ATAC"
multiome_filtered2_SCT_ATAC <- RunHarmony(multiome_filtered2_SCT_ATAC, group.by.vars = "sampleID", 
                                          reduction = "lsi", theta = 2, tau = 0, max.iter.harmony = 10, 
                                          max.iter.cluster = 20, epsilon.cluster = 1e-05, 
                                          epsilon.harmony = 1e-04, plot_convergence = TRUE, verbose = TRUE, 
                                          reduction.save = "harmonyATAC", assay.use = "ATAC", project.dim = FALSE)

mesenchymal2_integrated <- multiome_filtered2_SCT_ATAC

#WNN analysis
mesenchymal2_integrated <- FindMultiModalNeighbors(mesenchymal2_integrated, reduction.list = list("harmonySCT", "harmonyATAC"), 
                                                   dims.list = list(1:50, 2:40), k.nn = 50)
mesenchymal2_integrated <- RunUMAP(mesenchymal2_integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                                   reduction.key = "wnnUMAP_")
mesenchymal2_integrated <- FindClusters(mesenchymal2_integrated, graph.name = "wsnn", algorithm = 3, resolution = 0.5)

#mesenchymal subset 3 & processing
Idents(mesenchymal2_integrated) <- "seurat_clusters"
mesenchymal_cells3 <- subset(mesenchymal2_integrated, idents = c("0","1","2","3","4","5","6","7","8","9","10","13","14","15","16","17","19"))
#RNA processing
DefaultAssay(mesenchymal_cells3) <- "RNA"
mesenchymal_cells3_SCT <- SCTransform(mesenchymal_cells3, vars.to.regress = c("percent.mt","nCount_RNA"))
mesenchymal_cells3_SCT <- RunPCA(mesenchymal_cells3_SCT)
ElbowPlot(mesenchymal_cells3_SCT)
#ATAC processing
DefaultAssay(mesenchymal_cells3_SCT) <- "ATAC"
mesenchymal_cells3_SCT_ATAC <- FindTopFeatures(mesenchymal_cells3_SCT, min.cutoff = "q5")
mesenchymal_cells3_SCT_ATAC <- RunTFIDF(mesenchymal_cells3_SCT_ATAC)
mesenchymal_cells3_SCT_ATAC <- RunSVD(mesenchymal_cells3_SCT_ATAC)
#SCT assay integration
DefaultAssay(mesenchymal_cells3_SCT_ATAC) <- "SCT"
mesenchymal_cells3_SCT_ATAC <- RunHarmony(mesenchymal_cells3_SCT_ATAC, group.by.vars = "sampleID", 
                                          reduction = "pca", theta = 2, tau = 0, max.iter.harmony = 10, 
                                          max.iter.cluster = 20, epsilon.cluster = 1e-05, 
                                          epsilon.harmony = 1e-04, plot_convergence = TRUE, verbose = TRUE, 
                                          reduction.save = "harmonySCT", assay.use = "SCT", project.dim = FALSE)
#ATAC assay integration
DefaultAssay(mesenchymal_cells3_SCT_ATAC) <- "ATAC"
mesenchymal_cells3_SCT_ATAC <- RunHarmony(mesenchymal_cells3_SCT_ATAC, group.by.vars = "sampleID", 
                                          reduction = "lsi", theta = 2, tau = 0, max.iter.harmony = 10, 
                                          max.iter.cluster = 20, epsilon.cluster = 1e-05, 
                                          epsilon.harmony = 1e-04, plot_convergence = TRUE, verbose = TRUE, 
                                          reduction.save = "harmonyATAC", assay.use = "ATAC", project.dim = FALSE)

mesenchymal3_integrated <- mesenchymal_cells3_SCT_ATAC

#WNN analysis
mesenchymal3_integrated <- FindMultiModalNeighbors(mesenchymal3_integrated, reduction.list = list("harmonySCT", "harmonyATAC"), 
                                                   dims.list = list(1:10, 2:11), k.nn = 50)
mesenchymal3_integrated <- RunUMAP(mesenchymal3_integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                                   reduction.key = "wnnUMAP_")
mesenchymal3_integrated <- FindClusters(mesenchymal3_integrated, graph.name = "wsnn", algorithm = 3, resolution = 0.5)


# ── 5. Subset fibroblast cells ────────────────────────────────────────────────

#fibroblast subset 1 & processing
Idents(mesenchymal3_integrated) <- "seurat_clusters"
fibroblast_cells1 <- subset(mesenchymal3_integrated, idents = c("0","1","2","3","5","7","9","10","11","13"))

#add cell cycle scoring to dataset
#load cell cycle genes from Seurat
cc.genes <- Seurat::cc.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
#convert to mouse names
s.genes.mouse <- homologene(s.genes, inTax = 9606, outTax = 10090)
g2m.genes.mouse <- homologene(g2m.genes, inTax = 9606, outTax = 10090)
s.genes.mouse.vector <- as.vector(s.genes.mouse[[2]])
g2m.genes.mouse.vector <- as.vector(g2m.genes.mouse[[2]])
#run CellCycleScoring function
fibroblast_cells1 <- CellCycleScoring(fibroblast_cells1, s.features = s.genes.mouse.vector, g2m.features = g2m.genes.mouse.vector, set.ident = TRUE)

#RNA processing
DefaultAssay(fibroblast_cells1) <- "RNA"
fibroblast_cells1_SCT <- SCTransform(fibroblast_cells1, vars.to.regress = c("percent.mt","nCount_RNA","S.Score","G2M.Score"))
fibroblast_cells1_SCT <- RunPCA(fibroblast_cells1_SCT)
ElbowPlot(fibroblast_cells1_SCT)
#ATAC processing
DefaultAssay(fibroblast_cells1_SCT) <- "ATAC"
fibroblast_cells1_SCT_ATAC <- FindTopFeatures(fibroblast_cells1_SCT, min.cutoff = "q5")
fibroblast_cells1_SCT_ATAC <- RunTFIDF(fibroblast_cells1_SCT_ATAC)
fibroblast_cells1_SCT_ATAC <- RunSVD(fibroblast_cells1_SCT_ATAC)
#SCT assay integration
DefaultAssay(fibroblast_cells1_SCT_ATAC) <- "SCT"
fibroblast_cells1_SCT_ATAC <- RunHarmony(fibroblast_cells1_SCT_ATAC, group.by.vars = "sampleID", 
                                         reduction = "pca", theta = 2, tau = 0, max.iter.harmony = 10, 
                                         max.iter.cluster = 20, epsilon.cluster = 1e-05, 
                                         epsilon.harmony = 1e-04, plot_convergence = TRUE, verbose = TRUE, 
                                         reduction.save = "harmonySCT", assay.use = "SCT", project.dim = FALSE)
#ATAC assay integration
DefaultAssay(fibroblast_cells1_SCT_ATAC) <- "ATAC"
fibroblast_cells1_SCT_ATAC <- RunHarmony(fibroblast_cells1_SCT_ATAC, group.by.vars = "sampleID", 
                                         reduction = "lsi", theta = 2, tau = 0, max.iter.harmony = 10, 
                                         max.iter.cluster = 20, epsilon.cluster = 1e-05, 
                                         epsilon.harmony = 1e-04, plot_convergence = TRUE, verbose = TRUE, 
                                         reduction.save = "harmonyATAC", assay.use = "ATAC", project.dim = FALSE)

fibroblasts1_integrated <- fibroblast_cells1_SCT_ATAC

#WNN analysis
fibroblasts1_integrated <- FindMultiModalNeighbors(fibroblasts1_integrated, reduction.list = list("harmonySCT", "harmonyATAC"), 
                                                   dims.list = list(1:16, 2:16), k.nn = 50)
fibroblasts1_integrated <- RunUMAP(fibroblasts1_integrated, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
                                   reduction.key = "wnnUMAP_")
fibroblasts1_integrated <- FindClusters(fibroblasts1_integrated, graph.name = "wsnn", algorithm = 3, resolution = 0.4)

#final clustering chosen
fibroblasts_final <- fibroblasts1_integrated

#process RNA assay
DefaultAssay(fibroblasts_final) <- "RNA"
fibroblasts_final <- NormalizeData(fibroblasts_final, scale.factor = 10000)
fibroblasts_final <- FindVariableFeatures(fibroblasts_final, selection.method = "mean.var.plot", nfeatures = 2000)
fibroblasts_final <- ScaleData(fibroblasts_final, features = rownames(fibroblasts_final))


# ── 6. Save objects ───────────────────────────────────────────────────────────

saveRDS(multiome_integrated, "/path/to/multiome_integrated.rds")
saveRDS(fibroblasts_final,   "/path/to/fibroblasts_final.rds")
