# Human lung dataset figures - figure reproduction
# Reproduces the main and supplemental figure panels derived from the human single-cell
# datasets analyzed in human_sc_analysis_01.R
# Inputs: Seurat objects saved in human_sc_analysis_01.R

# ── 0. Load Packages ──────────────────────────────────────────────────────────

library(Seurat)
library(Signac)
library(cowplot)
library(patchwork)
library(ggplot2)

# ── 1. Load objects ───────────────────────────────────────────────────────────

#objects built by the human analysis script
fibroblasts_Habermann    <- readRDS("/path/to/fibroblasts_Habermann.rds")
fibroblast_Adams         <- readRDS("/path/to/fibroblast_Adams.rds")
fibroblasts_human_Tsukui <- readRDS("/path/to/fibroblasts_human_Tsukui.rds")
Valenzi_atac_fibroblasts <- readRDS("/path/to/Valenzi_atac_fibroblasts.rds")


# ── 2. Figure 7 ───────────────────────────────────────────────────────────────

#fig 7a
Idents(fibroblasts_Habermann) <- "Diagnosis"
levels(fibroblasts_Habermann) <- c("Control","IPF","cHP","Unclassifiable ILD","sarcoidosis","NSIP")

plots <- VlnPlot(fibroblasts_Habermann, features = c("matrix_sig_final_human1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) + NoLegend()
})
combined_plot <- plot_grid(plotlist = plots, ncol = 1, rel_heights = c(1), align = 'v', axis = 'lr')
combined_plot

#fig 7c
Idents(fibroblast_Adams) <- "Disease_Identity"

plots <- VlnPlot(fibroblast_Adams, features = c("matrix_sig_final_human1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) + NoLegend()
})
combined_plot <- plot_grid(plotlist = plots, ncol = 1, rel_heights = c(1), align = 'v', axis = 'lr')
combined_plot

#fig 7f
DimPlot(fibroblast_Adams, reduction = "umap", group.by = "Manuscript_Identity", label = FALSE)
DimPlot(fibroblast_Adams, reduction = "umap", group.by = "Disease_Identity", label = FALSE)

#fig 7g
DefaultAssay(fibroblast_Adams) <- "RNA"
FeaturePlot(fibroblast_Adams, features = "CTHRC1", reduction = "umap", label = F, order = F)
FeaturePlot(fibroblast_Adams, features = "matrix_sig_final_human1", reduction = "umap", label = F, order = F, min.cutoff = 0)

#fig 7i
DefaultAssay(Valenzi_atac_fibroblasts) <- "ATAC"
DimPlot(Valenzi_atac_fibroblasts, reduction = "umap_fib", group.by = "disease", label = FALSE)

DefaultAssay(Valenzi_atac_fibroblasts) <- "chromvar"
FeaturePlot(
  Valenzi_atac_fibroblasts, features  = "MA0490.2", reduction = "umap_fib", min.cutoff = "q5", max.cutoff = "q95", order = T) &
  scale_color_distiller(palette = "RdBu", direction = -1) &
  theme(axis.title = element_blank(), axis.text = element_blank())


# ── 3. Supplemental Figure 7 ──────────────────────────────────────────────────

#fig s7a
Idents(fibroblasts_human_Tsukui) <- "condition"
levels(fibroblasts_human_Tsukui) <- c("Ctrl","IPF")

plots <- VlnPlot(fibroblasts_human_Tsukui, features = c("matrix_sig_final_human1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) + NoLegend()
})
combined_plot <- plot_grid(plotlist = plots, ncol = 1, rel_heights = c(1), align = 'v', axis = 'lr')
combined_plot

#fig s7c
Idents(fibroblasts_Habermann) <- "Diagnosis"
levels(fibroblasts_Habermann) <- c("Control","IPF","cHP","Unclassifiable ILD","sarcoidosis","NSIP")

plots <- VlnPlot(fibroblasts_Habermann, features = c("all_mfb_ecm_genes_human1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) + NoLegend()
})
combined_plot <- plot_grid(plotlist = plots, ncol = 1, rel_heights = c(1), align = 'v', axis = 'lr')
combined_plot

#fig s7e
Idents(fibroblasts_human_Tsukui) <- "condition"
levels(fibroblasts_human_Tsukui) <- c("Ctrl","IPF")

plots <- VlnPlot(fibroblasts_human_Tsukui, features = c("all_mfb_ecm_genes_human1"), pt.size = 0.0, combine = FALSE, same.y.lims = TRUE)
plots <- lapply(plots, function(p) {
  p + stat_summary(fun.y = median, geom = 'point', size = 10, colour = "black", shape = 95) + NoLegend()
})
combined_plot <- plot_grid(plotlist = plots, ncol = 1, rel_heights = c(1), align = 'v', axis = 'lr')
combined_plot
