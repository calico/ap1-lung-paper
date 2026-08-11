# Mouse lung multiome - part 4 of 4: figures
# Reproduces the main and supplemental figure panels derived from the mouse single-nucleus
# multiome dataset and published scRNA-seq datasets, using the annotated RDS objects 
# provided in the GEO upload together with the objects produced by parts 2 and 3 of this analysis.
# Inputs and their origins are listed in section 1.

# ── 0. Load Packages ──────────────────────────────────────────────────────────

library(Seurat)
library(Signac)
library(Matrix)
library(cowplot)
library(GRaNIE)
library(igraph)
library(ggraph)
library(patchwork)
library(ggplot2)
library(dplyr)
library(tidyr)

# ── 1. Shared settings and objects ────────────────────────────────────────────

#display order used for every fibroblast violin, dot and bar plot below
celltype_order <- c("Activated myofibroblasts","Adult ductal myofibroblasts",
                    "P7 ductal myofibroblasts","Alveolar myofibroblasts",
                    "Matrix fibroblasts","Adventitial fibroblasts",
                    "Dev alveolar fibroblasts","Alveolar fibroblasts")

#annotated objects from the GEO submission
multiome_integrated <- readRDS("/path/to/multiome_integrated.rds")
fibroblasts_final   <- readRDS("/path/to/fibroblasts_final.rds")

#chromVAR and GRaNIE objects produced by part 2
fibroblasts_final_chrom         <- readRDS("/path/to/fibroblasts_final_chrom.rds")
multiome_integrated_downsampled <- readRDS("/path/to/multiome_integrated_downsampled.rds")
grn_object                      <- readRDS("/path/to/grn_object.rds")

#published mouse datasets processed in part 3
fibroblasts_Zepp    <- readRDS("/path/to/fibroblasts_Zepp.rds")
fibroblasts_Strunz  <- readRDS("/path/to/fibroblasts_Strunz.rds")
fibroblasts_Curras  <- readRDS("/path/to/fibroblasts_Curras.rds")
fibroblasts_Narvaez <- readRDS("/path/to/fibroblasts_Narvaez.rds")
fibroblasts_Tsukui  <- readRDS("/path/to/fibroblasts_Tsukui.rds")


# ── 2. Figure 1 ───────────────────────────────────────────────────────────────

#fig 1b
Idents(multiome_integrated) <- "cellType_specific_all"
DimPlot(multiome_integrated, reduction = "wnn.umap", label = FALSE, repel = TRUE, label.size = 2.5, raster = FALSE) + NoLegend()

#fig 1c
DimPlot(multiome_integrated, reduction = "wnn.umap", label = FALSE, repel = FALSE, label.size = 2.5, raster = FALSE, group.by = "condition_broad", shuffle = TRUE) + NoLegend() + ggtitle(NULL)

#fig 1d
cell_types <- c("Myofibroblasts", "Alveolar macrophages", "AT2 cells", "Pericytes")

meta <- multiome_integrated@meta.data %>% mutate(cellType = Idents(multiome_integrated))

df <- meta %>%
  group_by(sampleID, condition_broad) %>%
  mutate(total = n()) %>%
  filter(cellType %in% cell_types) %>%
  count(sampleID, condition_broad, total, cellType, name = "n") %>%
  ungroup() %>%
  tidyr::complete(cellType = cell_types, tidyr::nesting(sampleID, condition_broad, total), fill = list(n = 0)) %>%
  mutate(freq = n / total)

ggplot(df, aes(condition_broad, freq, fill = condition_broad)) +
  geom_boxplot(alpha = .7, outlier.size = 1) +
  geom_jitter(width = .2, size = 1, alpha = .7) +
  facet_wrap(~cellType, scales = "free_y") +
  theme_classic() +
  labs(x = "Condition", y = "Relative frequency of cells") +
  theme(legend.position = "none")


# ── 3. Supplemental Figure 2 ──────────────────────────────────────────────────

#fig s2a
Idents(multiome_integrated) <- "sampleID"
levels(multiome_integrated) <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20")

VlnPlot(multiome_integrated, features = c("nCount_RNA"), ncol = 1, pt.size = 0.001) + NoLegend() + ggtitle(NULL)
VlnPlot(multiome_integrated, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.001) + NoLegend() + ggtitle(NULL)
VlnPlot(multiome_integrated, features = c("percent.mt"), ncol = 1, pt.size = 0.001) + NoLegend() + ggtitle(NULL)
VlnPlot(multiome_integrated, features = c("nCount_ATAC"), ncol = 1, pt.size = 0.001) + NoLegend() + ggtitle(NULL)
VlnPlot(multiome_integrated, features = c("TSS.enrichment"), ncol = 1, pt.size = 0.001) + NoLegend() + ggtitle(NULL)
VlnPlot(multiome_integrated, features = c("nucleosome_signal"), ncol = 1, pt.size = 0.001) + NoLegend() + ggtitle(NULL)

#fig s2b
DensityScatter(multiome_integrated, x = "nFeature_RNA", y = "nCount_RNA", log_x = TRUE, log_y = TRUE, quantiles = TRUE) + ggtitle(NULL)
DensityScatter(multiome_integrated, x = "nFeature_RNA", y = "percent.mt", log_x = TRUE, log_y = FALSE, quantiles = TRUE) + ggtitle(NULL)
DensityScatter(multiome_integrated, x = "nCount_ATAC", y = "TSS.enrichment", log_x = TRUE, log_y = FALSE, quantiles = TRUE) + ggtitle(NULL)
DensityScatter(multiome_integrated, x = "nucleosome_signal", y = "TSS.enrichment", log_x = FALSE, log_y = FALSE, quantiles = TRUE) + ggtitle(NULL)

#fig s2c
DefaultAssay(multiome_integrated) <- "SCT"
FeaturePlot(multiome_integrated, reduction = "wnn.umap", features = "Cdh1", raster = FALSE, order = FALSE) + ggtitle(NULL)
FeaturePlot(multiome_integrated, reduction = "wnn.umap", features = "Pecam1", raster = FALSE, order = FALSE) + ggtitle(NULL)
FeaturePlot(multiome_integrated, reduction = "wnn.umap", features = "Col1a2", raster = FALSE, order = FALSE) + ggtitle(NULL)
FeaturePlot(multiome_integrated, reduction = "wnn.umap", features = "Ptprc", raster = FALSE, order = FALSE) + ggtitle(NULL)

#fig s2d
DefaultAssay(multiome_integrated) <- "SCT"
Idents(multiome_integrated) <- "cellType_specific_all"
levels(multiome_integrated) <- c("Interstitial macrophages","Alveolar macrophages","Dendritic cells","NK cells","B cells","T cells","Aerocytes","General capillary cells","PNEC cells","Club cells","Ciliated cells","Transitional cells","AT2 cells","AT1 cells","Mesothelial cells","Pericytes","SMCs","Myofibroblasts","Fibroblasts")
DotPlot(multiome_integrated, features = c("Col1a2","Pdgfra","Rbfox1","Acta2","Pdgfrb","Upk3b","Cdh1","Ager","Sftpc","Cldn4","Foxj1","Scgb1a1","Calca","Pecam1","Plvap","Car4","Ptprc","Cd3e","Cd19","Gzma","Itgax","Cd68","Marco","Cd163"))

#fig s2e
#pseudobulk gene expression for sample level PCA plot
expr <- GetAssayData(multiome_integrated, assay = "SCT", slot = "data")
samples <- as.character(multiome_integrated$sampleID)

pseudo_bulk <- sapply(unique(samples), function(s) {
  cells <- which(samples == s)
  Matrix::rowMeans(expr[, cells, drop = FALSE])
})

pseudo_bulk_t <- t(pseudo_bulk)
pseudo_bulk_t <- pseudo_bulk_t[, apply(pseudo_bulk_t, 2, var) > 0]

vars <- apply(pseudo_bulk_t, 2, var)
top_genes <- names(sort(vars, decreasing = TRUE))[1:2000]
pseudo_bulk_t <- pseudo_bulk_t[, top_genes]

#run PCA
pca <- prcomp(pseudo_bulk_t, scale. = TRUE)
rownames(pca$x) <- rownames(pseudo_bulk_t)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  sampleID = rownames(pca$x)
)

#plot
ggplot(pca_df, aes(x = PC1, y = PC2, label = sampleID)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.5) +
  theme_minimal()

#fig s2f
ggplot(multiome_integrated@meta.data, aes(x=multiome_integrated$condition, fill=cellType_specific_all)) + geom_bar(position = "fill")
ggplot(multiome_integrated@meta.data, aes(x=multiome_integrated$condition_broad, fill=cellType_specific_all)) + geom_bar(position = "fill")


# ── 4. Figure 2 ───────────────────────────────────────────────────────────────

#fig 2a
Idents(fibroblasts_final) <- "cellType_specific"
DimPlot(fibroblasts_final, reduction = "wnn.umap", label = FALSE, repel = TRUE, label.size = 2.5, raster = FALSE, pt.size = 0.25) + NoLegend()

#fig 2b
#cell type markers
DefaultAssay(fibroblasts_final) <- "SCT"
Idents(fibroblasts_final) <- "cellType_specific"
levels(fibroblasts_final) <- celltype_order
DotPlot(fibroblasts_final, features = c("Col1a2","Col4a1","Fn1","Pdgfra","Npnt","Wnt2","Dcn","Pi16","Rbfox1","Zfp385b","Hhip","Cdh4","Lgr6","Spp1","Cthrc1"))

#fig 2c
Idents(fibroblasts_final) <- "cellType_specific"
DimPlot(fibroblasts_final, reduction = "wnn.umap", label = FALSE, repel = TRUE, label.size = 2.5, raster = FALSE, split.by = "condition_broad", pt.size = 0.25) + NoLegend()

#fig 2d
Idents(fibroblasts_final) <- "cellType_specific"
levels(fibroblasts_final) <- celltype_order
DefaultAssay(fibroblasts_final) <- "RNA"

plots <- VlnPlot(fibroblasts_final, features = c("Curras_fibrosis_sig1","Spp1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 4, colour = "black", shape = 95) + NoLegend() + theme(axis.text.x = element_text(size = 5))
})
combined_plot <- plot_grid(plotlist = plots, ncol = 2, rel_heights = c(1,1), align = 'v', axis = 'lr')
combined_plot

#fig 2e
DefaultAssay(fibroblasts_final) <- "RNA"
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Curras_fibrosis_sig1", raster = FALSE, order = FALSE, min.cutoff = 0.0) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Spp1", raster = FALSE, order = FALSE) + ggtitle(NULL)

#fig 2f
DefaultAssay(fibroblasts_final) <- "RNA"
Idents(fibroblasts_final) <- "cellType_specific"
levels(fibroblasts_final) <- celltype_order

plots <- VlnPlot(fibroblasts_final, features = c("Eln","Rbfox1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 4, colour = "black", shape = 95) + NoLegend() + theme(axis.text.x = element_text(size = 5))
})
combined_plot <- plot_grid(plotlist = plots, ncol = 2, rel_heights = c(1,1), align = 'v', axis = 'lr')
combined_plot

#fig 2g
DefaultAssay(fibroblasts_final) <- "RNA"
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Eln", raster = FALSE, order = FALSE) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Rbfox1", raster = FALSE, order = FALSE) + ggtitle(NULL)


# ── 5. Supplemental Figure 3 ──────────────────────────────────────────────────

#fig s3a
Idents(fibroblasts_final) <- "cellType_specific"
levels(fibroblasts_final) <- celltype_order
ggplot(fibroblasts_final@meta.data, aes(x=fibroblasts_final$condition_broad, fill=cellType_specific)) + geom_bar(position = "fill")

#fig s3b
Idents(fibroblasts_final) <- "cellType_specific"
levels(fibroblasts_final) <- celltype_order
DefaultAssay(fibroblasts_final) <- "RNA"

plots <- VlnPlot(fibroblasts_final, features = c("Tsukui_fibrosis_sig1","Fn1","Pdgfra","Loxl2","Hhip","Fbn2"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 4, colour = "black", shape = 95) + NoLegend() + theme(axis.text.x = element_text(size = 5))
})
combined_plot <- plot_grid(plotlist = plots, ncol = 6, rel_heights = c(1,1,1,1,1,1), align = 'v', axis = 'lr')
combined_plot

#fig s3c
DefaultAssay(fibroblasts_final) <- "RNA"
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Tsukui_fibrosis_sig1", raster = FALSE, order = FALSE, min.cutoff = 0) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Pdgfra", raster = FALSE, order = FALSE) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Hhip", raster = FALSE, order = FALSE) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Loxl2", raster = FALSE, order = FALSE) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Fn1", raster = FALSE, order = FALSE) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Fbn2", raster = FALSE, order = FALSE) + ggtitle(NULL)

#fig s3e
Idents(fibroblasts_final) <- "cellType_specific"
levels(fibroblasts_final) <- celltype_order
DefaultAssay(fibroblasts_final) <- "RNA"

plots <- VlnPlot(fibroblasts_final, features = c("Col3a1","Col4a1","Col6a2","Bgn"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 4, colour = "black", shape = 95) + NoLegend() + theme(axis.text.x = element_text(size = 5))
})
combined_plot <- plot_grid(plotlist = plots, ncol = 4, rel_heights = c(1,1,1,1), align = 'v', axis = 'lr')
combined_plot

#fig s3f
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Col3a1", raster = FALSE, order = FALSE) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Col4a1", raster = FALSE, order = FALSE) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Col6a2", raster = FALSE, order = FALSE) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Bgn", raster = FALSE, order = FALSE) + ggtitle(NULL)

#fig s3g
Idents(fibroblasts_Narvaez) <- "timepoint"
DefaultAssay(fibroblasts_Narvaez) <- "RNA"
plots <- VlnPlot(fibroblasts_Narvaez, features = c("all_mfb_ecm_genes1","matrix_genes_final1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) + NoLegend()
})
combined_plot <- plot_grid(plotlist = plots, ncol = 2, rel_heights = c(1,1), align = 'v', axis = 'lr')
combined_plot

#fig s3h
Idents(fibroblasts_Zepp) <- "orig.ident"
DefaultAssay(fibroblasts_Zepp) <- "RNA"

plots <- VlnPlot(fibroblasts_Zepp, features = c("matrix_genes_final1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
    coord_cartesian(ylim = c(-0.5, 1)) + NoLegend()
})
combined_plot <- plot_grid(plotlist = plots, ncol = 2, rel_heights = c(1,1), align = 'v', axis = 'lr')
combined_plot

#fig s3k
Idents(fibroblasts_Tsukui) <- "condition"
DefaultAssay(fibroblasts_Tsukui) <- "RNA"
plots <- VlnPlot(fibroblasts_Tsukui, features = c("all_mfb_ecm_genes1","matrix_genes_final1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) + NoLegend()
})
combined_plot <- plot_grid(plotlist = plots, ncol = 2, rel_heights = c(1,1), align = 'v', axis = 'lr')
combined_plot


# ── 6. Figure 3 ───────────────────────────────────────────────────────────────

#fig 3a
Idents(fibroblasts_final) <- "cellType_specific"
levels(fibroblasts_final) <- celltype_order

#violin plot
plots <- VlnPlot(fibroblasts_final, features = c("matrix_genes_final1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 4, colour = "black", shape = 95) + NoLegend() + theme(axis.text.x = element_text(size = 5))
})
combined_plot <- plot_grid(plotlist = plots, ncol = 3, rel_heights = c(1), align = 'v', axis = 'lr')
combined_plot

#feature plot
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "matrix_genes_final1", raster = FALSE, order = FALSE, min.cutoff = 0.0) + ggtitle(NULL)

#fig 3b
Idents(fibroblasts_final) <- "cellType_specific"
levels(fibroblasts_final) <- celltype_order

#violin plot
plots <- VlnPlot(fibroblasts_final, features = c("all_mfb_ecm_genes1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 4, colour = "black", shape = 95) + NoLegend() + theme(axis.text.x = element_text(size = 5))
})
combined_plot <- plot_grid(plotlist = plots, ncol = 3, rel_heights = c(1), align = 'v', axis = 'lr')
combined_plot

#feature plot
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "all_mfb_ecm_genes1", raster = FALSE, order = FALSE, min.cutoff = 0.0) + ggtitle(NULL)

#fig 3c
Idents(fibroblasts_Zepp) <- "orig.ident"
DefaultAssay(fibroblasts_Zepp) <- "RNA"

plots <- VlnPlot(fibroblasts_Zepp, features = c("all_mfb_ecm_genes1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) +
    coord_cartesian(ylim = c(-0.5, 1)) + NoLegend()
})
combined_plot <- plot_grid(plotlist = plots, ncol = 2, rel_heights = c(1,1), align = 'v', axis = 'lr')
combined_plot

#fig 3d
DefaultAssay(fibroblasts_Strunz) <- "RNA"
Idents(fibroblasts_Strunz) <- "grouping"
levels(fibroblasts_Strunz) <- c("PBS","d3","d7","d10","d14","d21","d28")

plots <- VlnPlot(fibroblasts_Strunz, features = c("all_mfb_ecm_genes1","matrix_genes_final1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) + NoLegend()
})
combined_plot <- plot_grid(plotlist = plots, ncol = 2, rel_heights = c(1, 1), align = 'v', axis = 'lr')
combined_plot

#fig 3e
Idents(fibroblasts_Curras) <- "early_late"
DefaultAssay(fibroblasts_Curras) <- "RNA"
plots <- VlnPlot(fibroblasts_Curras, features = c("all_mfb_ecm_genes1","matrix_genes_final1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) + NoLegend()
})
combined_plot <- plot_grid(plotlist = plots, ncol = 2, rel_heights = c(1,1), align = 'v', axis = 'lr')
combined_plot


# ── 7. Figure 4 ───────────────────────────────────────────────────────────────

#fig 4a
#add metadata for coloring
Idents(fibroblasts_final) <- "cellType_specific"
fibroblasts_final@meta.data[which(fibroblasts_final@meta.data$cellType_specific %in% "Alveolar myofibroblasts"), "cellType_specific2"] <- "P7 myofibroblasts"
fibroblasts_final@meta.data[which(fibroblasts_final@meta.data$cellType_specific %in% "P7 ductal myofibroblasts"), "cellType_specific2"] <- "P7 myofibroblasts"
fibroblasts_final@meta.data[which(fibroblasts_final@meta.data$cellType_specific %in% "Activated myofibroblasts"), "cellType_specific2"] <- "Activated myofibroblasts"

#plot
Idents(fibroblasts_final) <- "cellType_specific2"
DimPlot(fibroblasts_final, reduction = "wnn.umap", label = F, repel = TRUE, label.size = 2.5, raster = FALSE) + NoLegend()

#fig 4d
DefaultAssay(fibroblasts_final_chrom) <- "chromvar"
Idents(fibroblasts_final_chrom) <- "cellType_specific"
levels(fibroblasts_final_chrom) <- celltype_order

#violin plots
plots <- VlnPlot(fibroblasts_final_chrom, features = c("MA0490.2","MA0840.1"), pt.size = 0.0, combine = FALSE, same.y.lims = FALSE) #junb, creb5
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 4, colour = "black", shape = 95) + NoLegend() + theme(axis.text.x = element_text(size = 10))
})
combined_plot <- plot_grid(plotlist = plots, ncol = 2, rel_heights = c(1,1), align = 'v', axis = 'lr')
combined_plot

#feature plots
FeaturePlot(fibroblasts_final_chrom, reduction = "wnn.umap", features = "MA0490.2", raster = FALSE, label.size = 2, order = F, label = F, min.cutoff = 'q05', max.cutoff = 'q95', repel = TRUE) #junb
FeaturePlot(fibroblasts_final_chrom, reduction = "wnn.umap", features = "MA0840.1", raster = FALSE, label.size = 2, order = F, label = F, min.cutoff = 'q05', max.cutoff = 'q95', repel = TRUE) #creb5

#fig 4e
DefaultAssay(fibroblasts_final) <- "RNA"
Idents(fibroblasts_final) <- "cellType_specific"
levels(fibroblasts_final) <- celltype_order

#violin plot
plots <- VlnPlot(fibroblasts_final, features = c("Junb","Creb5"), pt.size = 0.0, combine = FALSE, same.y.lims = FALSE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 4, colour = "black", shape = 95) + NoLegend() + theme(axis.text.x = element_text(size = 5))
})
combined_plot <- plot_grid(plotlist = plots, ncol = 2, rel_heights = c(1,1), align = 'v', axis = 'lr')
combined_plot

#feature plots
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Junb", raster = FALSE, order = FALSE, min.cutoff = 0) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Creb5", raster = FALSE, order = FALSE, min.cutoff = 0) + ggtitle(NULL)

#fig 4f
plot_tf_network <- function(grn_object, tf_name, top_n_targets = 100, filter_fdr = 0.25) {
  connections <- GRaNIE::getGRNConnections(grn_object)
  #filter to this TF and significant peak-gene links, keep strongest correlations
  tf_network <- connections %>%
    filter(as.character(TF.name) == tf_name,
           peak_gene.p_adj < filter_fdr) %>%
    arrange(desc(abs(peak_gene.r))) %>%
    slice_head(n = top_n_targets)

  if (nrow(tf_network) == 0) return(NULL)

  edges <- tf_network %>%
    select(from = TF.name, to = gene.name, weight = peak_gene.r, fdr = peak_gene.p_adj)
  nodes <- data.frame(
    name = unique(c(edges$from, edges$to)),
    type = ifelse(unique(c(edges$from, edges$to)) == tf_name, "TF", "Gene")
  )
  g <- graph_from_data_frame(edges, directed = TRUE, vertices = nodes)

  ggraph(g, layout = 'star') +
    geom_edge_link(aes(edge_width = abs(weight), edge_color = weight),
                   arrow = arrow(length = unit(3, 'mm')),
                   end_cap = circle(3, 'mm'), alpha = 0.6) +
    geom_node_point(aes(color = type, size = type)) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3) +
    scale_edge_color_gradient2(low = "blue", mid = "gray", high = "red",
                               midpoint = 0, name = "Correlation") +
    scale_edge_width(range = c(0.5, 2), name = "Abs. Corr.") +
    scale_color_manual(values = c("TF" = "#E64B35", "Gene" = "#4DBBD5")) +
    scale_size_manual(values = c("TF" = 8, "Gene" = 2)) +
    labs(title = paste0(tf_name, " Regulatory Network"),
         subtitle = paste0("Top ", nrow(edges), " target genes (FDR < ", filter_fdr, ")")) +
    theme_graph() +
    theme(legend.position = "right")
}

#generate plots for both TFs
p_creb5 <- plot_tf_network(grn_object, "Creb5", top_n_targets = 50)
p_junb  <- plot_tf_network(grn_object, "Junb",  top_n_targets = 35)
p_creb5 | p_junb

#fig 4g
Idents(fibroblasts_final) <- "cellType_specific"
levels(fibroblasts_final) <- celltype_order

plots <- VlnPlot(fibroblasts_final, features = c("junb_genes1","creb5_genes1"), pt.size = 0.0, combine = FALSE, same.y.lims = FALSE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 4, colour = "black", shape = 95) + NoLegend() + theme(axis.text.x = element_text(size = 5))
})
combined_plot <- plot_grid(plotlist = plots, ncol = 2, rel_heights = c(1,1), align = 'v', axis = 'lr')
combined_plot

#fig 4h
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "junb_genes1", raster = FALSE, order = FALSE, min.cutoff = 0) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "creb5_genes1", raster = FALSE, order = FALSE, min.cutoff = 0) + ggtitle(NULL)


# ── 8. Supplemental Figure 4 ──────────────────────────────────────────────────

#fig s4a
DefaultAssay(fibroblasts_final_chrom) <- "chromvar"
FeaturePlot(fibroblasts_final_chrom, reduction = "wnn.umap", features = "MA0489.1", raster = FALSE, label.size = 2, order = F, label = F, min.cutoff = 'q05', max.cutoff = 'q95', repel = TRUE) #jun
FeaturePlot(fibroblasts_final_chrom, reduction = "wnn.umap", features = "MA0491.2", raster = FALSE, label.size = 2, order = F, label = F, min.cutoff = 'q05', max.cutoff = 'q95', repel = TRUE) #jund
FeaturePlot(fibroblasts_final_chrom, reduction = "wnn.umap", features = "MA0476.1", raster = FALSE, label.size = 2, order = F, label = F, min.cutoff = 'q05', max.cutoff = 'q95', repel = TRUE) #fos
FeaturePlot(fibroblasts_final_chrom, reduction = "wnn.umap", features = "MA0477.2", raster = FALSE, label.size = 2, order = F, label = F, min.cutoff = 'q05', max.cutoff = 'q95', repel = TRUE) #fosl1
FeaturePlot(fibroblasts_final_chrom, reduction = "wnn.umap", features = "MA0478.1", raster = FALSE, label.size = 2, order = F, label = F, min.cutoff = 'q05', max.cutoff = 'q95', repel = TRUE) #fosl2

DefaultAssay(fibroblasts_final) <- "RNA"
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Jun", raster = FALSE, order = FALSE, min.cutoff = 0) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Jund", raster = FALSE, order = FALSE, min.cutoff = 0) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Fos", raster = FALSE, order = FALSE, min.cutoff = 0) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Fosl1", raster = FALSE, order = FALSE, min.cutoff = 0) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Fosl2", raster = FALSE, order = FALSE, min.cutoff = 0) + ggtitle(NULL)

#fig s4c
DefaultAssay(fibroblasts_final_chrom) <- "chromvar"
FeaturePlot(fibroblasts_final_chrom, reduction = "wnn.umap", features = "MA1632.1", raster = FALSE, label.size = 2, order = F, label = F, min.cutoff = 'q05', max.cutoff = 'q95', repel = TRUE) #atf2
FeaturePlot(fibroblasts_final_chrom, reduction = "wnn.umap", features = "MA0834.1", raster = FALSE, label.size = 2, order = F, label = F, min.cutoff = 'q05', max.cutoff = 'q95', repel = TRUE) #atf7

DefaultAssay(fibroblasts_final) <- "RNA"
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Atf2", raster = FALSE, order = FALSE, min.cutoff = 0) + ggtitle(NULL)
FeaturePlot(fibroblasts_final, reduction = "wnn.umap", features = "Atf7", raster = FALSE, order = FALSE, min.cutoff = 0) + ggtitle(NULL)

#fig s4d
Idents(multiome_integrated_downsampled) <- "cellType_specific_all"
DefaultAssay(multiome_integrated_downsampled) <- "chromvar"

plots <- VlnPlot(multiome_integrated_downsampled, features = c("MA0490.2","MA0840.1"), pt.size = 0.0, combine = FALSE, same.y.lims = FALSE) #junb, creb5
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 4, colour = "black", shape = 95) + NoLegend() + theme(axis.text.x = element_text(size = 5))
})
combined_plot <- plot_grid(plotlist = plots, ncol = 2, rel_heights = c(1,1), align = 'v', axis = 'lr')
combined_plot

#fig s4g
Idents(fibroblasts_final) <- "cellType_specific"
levels(fibroblasts_final) <- celltype_order

plots <- VlnPlot(fibroblasts_final, features = c("fosl2_genes1","fos_genes1","jun_genes1"), pt.size = 0.0, combine = FALSE, same.y.lims = FALSE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 4, colour = "black", shape = 95) + NoLegend() + theme(axis.text.x = element_text(size = 5))
})
combined_plot <- plot_grid(plotlist = plots, ncol = 3, rel_heights = c(1,1,1), align = 'v', axis = 'lr')
combined_plot

#fig s4h
plot_genomic_locus <- function(grn, tf_name, gene_name) {

  #get the specific TF-gene connection
  connections <- GRaNIE::getGRNConnections(grn) %>%
    filter(as.character(TF.name) == tf_name, as.character(gene.name) == gene_name)

  if (nrow(connections) == 0) {
    stop(paste("No connection found between", tf_name, "and", gene_name))
  }

  #peak and gene annotation
  peaks <- grn@annotation$peaks
  genes <- grn@annotation$genes
  gene_info <- genes %>% filter(as.character(gene.name) == gene_name)

  #parse peak coordinates from peak.ID (format: chr1:start-end)
  peak_info <- peaks %>%
    filter(peak.ID %in% connections$peak.ID) %>%
    mutate(
      chr   = sub(":.*", "", peak.ID),
      start = as.numeric(sub(".*:(\\d+)-\\d+", "\\1", peak.ID)),
      end   = as.numeric(sub(".*:(\\d+)-(\\d+)", "\\2", peak.ID))
    )

  chr       <- as.character(gene_info$gene.chr[1])
  start_pos <- min(gene_info$gene.start, peak_info$start) - 10000
  end_pos   <- max(gene_info$gene.end,   peak_info$end)   + 10000

  #TSS depends on strand (minus-strand genes have TSS at gene.end)
  gene_tss <- if (as.character(gene_info$gene.strand[1]) == "-") {
    gene_info$gene.end[1]
  } else {
    gene_info$gene.start[1]
  }

  ggplot() +
    #gene track
    geom_rect(data = gene_info,
              aes(xmin = gene.start, xmax = gene.end, ymin = 0.8, ymax = 1.2),
              fill = "#4DBBD5", color = "black") +
    geom_text(data = gene_info,
              aes(x = (gene.start + gene.end)/2, y = 1.5, label = gene.name),
              size = 4, fontface = "bold") +
    #peak track
    geom_rect(data = peak_info,
              aes(xmin = start, xmax = end, ymin = 0.2, ymax = 0.6),
              fill = "#E64B35", color = "black", alpha = 0.7) +
    #connection arcs - one per peak, from peak midpoint to gene TSS
    geom_curve(data = peak_info,
               aes(x = (start + end)/2, xend = gene_tss, y = 0.4, yend = 1.0),
               curvature = -0.3, arrow = arrow(length = unit(2, 'mm')),
               color = "darkgray", linewidth = 1) +
    annotate("text", x = start_pos + 5000, y = 1.0,
             label = "Gene", hjust = 0, fontface = "bold") +
    annotate("text", x = start_pos + 5000, y = 0.4,
             label = paste0(tf_name, " Peak"), hjust = 0, fontface = "bold") +
    scale_y_continuous(limits = c(0, 2), expand = c(0, 0)) +
    scale_x_continuous(labels = scales::comma) +
    labs(title = paste0(tf_name, " -> ", gene_name, " Regulatory Connection"),
         subtitle = paste0(chr, ":", format(start_pos, big.mark = ","), "-",
                           format(end_pos, big.mark = ",")),
         x = "Genomic Position (bp)", y = "") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(face = "bold", size = 14))
}

p_creb5_eln <- plot_genomic_locus(grn_object, "Creb5", "Eln")
p_creb5_eln


# ── 9. Figure 6 ───────────────────────────────────────────────────────────────

#fig 6f
DefaultAssay(fibroblasts_final) <- "RNA"
Idents(fibroblasts_final) <- "cellType_specific"
levels(fibroblasts_final) <- celltype_order

plots <- VlnPlot(fibroblasts_final, features = c("ERK_targets1", "Trametinib_up1"), pt.size = 0.0, combine = FALSE, same.y.lims = FALSE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 4, colour = "black", shape = 95) + NoLegend() + theme(axis.text.x = element_text(size = 5))
})
combined_plot <- plot_grid(plotlist = plots, ncol = 1, rel_heights = c(1,1), align = 'v', axis = 'lr')
combined_plot

