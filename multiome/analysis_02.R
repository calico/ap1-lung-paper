# Mouse lung multiome - part 2 of 4: analysis of multiome RDS objects
# chromVAR motif activity, differential expression and accessibility, module scoring,
# pseudobulking, and GRaNIE gene regulatory network inference.
#
# Inputs: the annotated multiome_integrated and fibroblasts_final objects deposited with
#         the GEO submission
#         HOCOMOCO v12 PWMScan TFBS bed files and translationTable_mm10.csv for GRaNIE
# Outputs: fibroblasts_final_chrom, multiome_integrated_downsampled, marker and differential
#          accessibility tables, pseudobulk matrices, and grn_object

# ── 0. Load Packages ──────────────────────────────────────────────────────────

library(Seurat)
library(Signac)
library(GenomeInfoDb)
library(GenomicRanges)
library(BSgenome.Mmusculus.UCSC.mm10)
library(org.Mm.eg.db)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(chromVAR)
library(GRaNIE)
library(dplyr)

# ── 1. Load the annotated objects ─────────────────────────────────────────────

#annotated objects from the GEO submission
multiome_integrated <- readRDS("/path/to/multiome_integrated.rds")
fibroblasts_final   <- readRDS("/path/to/fibroblasts_final.rds")


# ── 2. chromVAR analysis ──────────────────────────────────────────────────────

#adding chromVAR analysis to fibroblasts_final
DefaultAssay(fibroblasts_final) <- "ATAC"
valid_chroms <- seqlevels(BSgenome.Mmusculus.UCSC.mm10)
valid_peaks <- which(seqnames(fibroblasts_final[["ATAC"]]@ranges) %in% valid_chroms)

#subset the counts matrix and the peak ranges
fibroblasts_final_chrom <- fibroblasts_final
fibroblasts_final_chrom[["ATAC"]] <- subset(fibroblasts_final[["ATAC"]], features = rownames(fibroblasts_final[["ATAC"]])[valid_peaks])

#get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(x = JASPAR2020, opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE))

#add motif information
fibroblasts_final_chrom <- AddMotifs(object = fibroblasts_final_chrom, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)

#run chromVAR
fibroblasts_final_chrom <- RunChromVAR(object = fibroblasts_final_chrom, genome = BSgenome.Mmusculus.UCSC.mm10)


#adding chromVAR analysis to multiome_integrated
#downsample multiome_integrated for chromvar to run efficiently
downsampled_n_cells <- 50000 

if (ncol(multiome_integrated) > downsampled_n_cells) {
  set.seed(123)
  cells_to_keep <- sample(Cells(multiome_integrated), size = downsampled_n_cells, replace = FALSE)
  multiome_integrated_downsampled <- subset(multiome_integrated, cells = cells_to_keep)
  print(paste("Seurat object downsampled to", ncol(multiome_integrated_downsampled), "cells."))
} else {
  print("Seurat object already has fewer or equal cells than the downsampling target.")
  multiome_integrated_downsampled <- multiome_integrated
}

DefaultAssay(multiome_integrated_downsampled) <- "ATAC"
#add motif information
multiome_integrated_downsampled <- AddMotifs(object = multiome_integrated_downsampled, genome = BSgenome.Mmusculus.UCSC.mm10, pfm = pfm)

valid_peaks <- which(seqnames(multiome_integrated_downsampled[["ATAC"]]@ranges) %in% valid_chroms)

#subset the counts matrix and the peak ranges
multiome_integrated_downsampled[["ATAC"]] <- subset(multiome_integrated_downsampled[["ATAC"]], features = rownames(multiome_integrated_downsampled[["ATAC"]])[valid_peaks])

#run chromVAR
multiome_integrated_downsampled <- RunChromVAR(object = multiome_integrated_downsampled, genome = BSgenome.Mmusculus.UCSC.mm10)


# ── 3. Differential expression/accessibility analysis ─────────────────────────

#identify differentially expressed genes between developmental and activated myofibroblasts
DefaultAssay(fibroblasts_final) <- "RNA"
Idents(fibroblasts_final) <- "cellType_specific"
markers <- FindMarkers(object = fibroblasts_final, ident.1 = "Activated myofibroblasts", ident.2 = c("Alveolar myofibroblasts","P7 ductal myofibroblasts"), test.use = 'wilcox', min.pct = 0.01)

#identify differentially expressed genes between activated myofibroblasts and the rest of the fibroblast dataset
DefaultAssay(fibroblasts_final) <- "RNA"
Idents(fibroblasts_final) <- "cellType_specific"
markers <- FindMarkers(object = fibroblasts_final, ident.1 = "Activated myofibroblasts", test.use = 'wilcox', min.pct = 0.01)

#identify differentially expressed peaks among all fibroblasts
DefaultAssay(fibroblasts_final) <- "ATAC"
Idents(fibroblasts_final) <- "cellType_specific"
atac_markers <- FindAllMarkers(object = fibroblasts_final, only.pos = TRUE, test.use = 'LR', min.pct = 0.05, latent.vars = 'nCount_ATAC')

#identify differentially expressed peaks between developmental and activated myofibroblasts
DefaultAssay(fibroblasts_final) <- "ATAC"
Idents(fibroblasts_final) <- "cellType_specific"
da_peaks <- FindMarkers(object = fibroblasts_final, ident.1 = c("Alveolar myofibroblasts","P7 ductal myofibroblasts"), ident.2 = "Activated myofibroblasts", only.pos = FALSE, test.use = 'LR', min.pct = 0.05, latent.vars = 'nCount_ATAC')


# ── 4. Gene expression signatures ─────────────────────────────────────────────

#ecm signatures
#fibrosis signature developed from Tsukui et al.
Tsukui_fibrosis_sig <- c("Actn1",	"Adamts2",	"Arpc1b",	"Atp6v0e",	"Bax",	"Calm1",	"Ccdc80",	"Ccnd1",	"Ccnd2",	"Cd44",	"Cd63",	"Cd74",	"Cdk6",	"Cdkn1a",	"Chd4",	"Cilp",	"Col12a1",	"Col14a1",	"Col15a1",	"Col16a1",	"Col18a1",	"Col1a1",	"Col1a2",	"Col3a1",	"Col4a1",	"Col4a2",	"Col5a1",	"Col5a2",	"Col6a2",	"Col6a3",	"Csrp2",	"Ctsk",	"Cxcl12",	"Dag1",	"Dbi",	"Ddr2",	"Dpt",	"Efemp2",	"Eln",	"Ext1",	"Fbn1",	"Fn1",	"Fst",	"Fzd1",	"Gas6",	"Gnai2",	"Gng12",	"Grb10",	"Gucy1a1",	"H2-D1",	"H2-K1",	"H2-Q7",	"Hif1a",	"Hsp90b1",	"Hspg2",	"Igf1",	"Igfbp5",	"Il1r1",	"Iqgap1",	"Isg15",	"Itga1",	"Itgb1",	"Itgb5",	"Jchain",	"Kras",	"Loxl2",	"Ltbp1",	"Ltbp2",	"Mdk",	"Mef2c",	"Mfap2",	"Mfap5",	"Mmp14",	"Mmp2",	"Mmp23",	"Myh11",	"Myh9",	"Myl6",	"Ncor1",	"Ndnf",	"Nid1",	"Nrp2",	"Nrtn",	"P4hb",	"Pdgfrb",	"Pfn1",	"Pik3r1",	"Plat",	"Plod2",	"Postn",	"Pten",	"Rock1",	"Runx1",	"Serpinh1",	"Sh3pxd2a",	"Slit2",	"Sparc",	"Spp1",	"Spred1",	"Stat1",	"Sulf1",	"Tcf4",	"Tgfbr2",	"Thbs2",	"Timp1",	"Tnc",	"Tpm1",	"Vcan")

#fibrosis signature developed from Curras-Alonso et al.
Curras_fibrosis_sig <- c("Bgn",	"Cd44",	"Col3a1",	"Col4a2",	"Col5a1",	"Col5a2",	"Col1a1",	"Col1a2",	"Ctsd",	"Ctsk",	"Ctss",	"Eln",	"Fbn1",	"Fn1",	"Ltbp2",	"Sparc",	"Spp1",	"Sdc3",	"Timp1",	"Tnc",	"Emilin1",	"Il4ra",	"Tnfrsf1a",	"Hspg2",	"Cdkn1a",	"Igf1",	"Osmr",	"Bcl3",	"H2-D1",	"H2-K1",	"H2-Q7",	"Isg15",	"Tgfbr2",	"Sphk1",	"Arpc1b",	"B2m",	"Casp4",	"Fat1",	"Arpc2",	"C3",	"Tnfrsf12a",	"Ecscr",	"Runx1",	"Dbi",	"Tyrobp")

#developmental ECM signature
all_mfb_ecm_genes <- c("Adam12", "Adam19", "Adamts10", "Adamts17", "Adamts6", "Adamts9", "Aspn", "Bmp2", "Cask", "Col14a1", "Col24a1", "Col25a1", "Col27a1", "Dst", "Eln", "Eng", "Fbln1", "Fbln5", "Fbn2", "Hmcn1", "Hpse2", "Htra1", "Itga1", "Itga9", "Itgav", "Itgb1", "Lrp12", "Ltbp2", "Mfap2", "Ntn4", "P3h2", "P4ha3", "Ppib", "Reck", "Sdc2", "Sparc", "Tgfb2", "Tgfbi")

#fibrotic ECM signature
matrix_genes_final <- c("Adamts12","Adamts15","Bgn","Bmp1","Col16a1","Col1a1","Col1a2","Col23a1","Col5a2","Col5a3","Col7a1","Ctsd","Fbn1","Fgf2","Fn1","Has2","Hspg2","Itga5","Lama5","Ndnf","Plec","Serpine1","Spp1","Thbs1","Timp1","Tnc")

DefaultAssay(fibroblasts_final) <- "RNA"
fibroblasts_final <- AddModuleScore(fibroblasts_final, features = list(Tsukui_fibrosis_sig), name = "Tsukui_fibrosis_sig")
fibroblasts_final <- AddModuleScore(fibroblasts_final, features = list(Curras_fibrosis_sig), name = "Curras_fibrosis_sig")
fibroblasts_final <- AddModuleScore(fibroblasts_final, features = list(all_mfb_ecm_genes), name = "all_mfb_ecm_genes")
fibroblasts_final <- AddModuleScore(fibroblasts_final, features = list(matrix_genes_final), name = "matrix_genes_final")

#grn signatures
creb5_genes <- c(
  "Mark1", "Ggta1", "Bcas1", "Bcas1os2", "Zfp704", "Rbis", "Gm10710", "Dchs2", "Nes", "Bach2",
  "Ror1", "Kcnq4", "Tex46", "4933427I22Rik", "Gm13032", "Slc45a1", "Dancr", "Rasl11b", "Tbx3os1",
  "Gm16063", "Tbx3", "Eln", "Stx1a", "Vps37d", "Plxna4os1", "Ifitm10", "Lsp1", "Ankrd10", "Adcy7",
  "Brd7", "Arhgap42", "Gm16833", "Plekhg1", "Samd5", "Ros1", "Unc5b", "Myocd", "Gm12295",
  "Cacna1g", "Hoxb8", "Hoxb6", "Hoxb4", "Hoxb2", "Skap1", "Meox1", "Mxra7", "C1qtnf1",
  "F730043M19Rik", "Vash1", "Net1", "Gm8765", "Tgfbi", "Smad5", "Atp6ap1l", "Lhfpl2", "Tnpo1",
  "Map1b", "Wnt5a", "Npy4r", "Tnfrsf19", "Defb42", "Fam167a", "Dock5", "Farp1", "Ghr", "Hoxc9",
  "Cd200", "Filip1l", "Pde10a", "Qpct", "Megf10", "Tjp2", "Pdlim1", "Entpd1", "Aldh18a1",
  "Neurl1a", "Sh3pxd2a", "Gm19557", "Rdh10", "Olfml2b", "Ccdc3", "Platr4", "Syt11", "Zfp462",
  "Met", "Pdzrn3", "P4ha3", "Zfp827", "Ubash3b", "2900052N01Rik", "Smad6", "Ulk4", "4930432O09Rik",
  "Fgf18", "Hoxb1", "Rundc3a", "Sox11", "Egln3", "Slc24a4", "Kif26a", "Sema4d", "Gm38397", "Arl11",
  "Dleu2", "Kcnrg", "Klf10", "Ptk2", "Fndc1", "Erlin1"
)

junb_genes <- c(
  "Ckmt1", "S100a6", "Ptpn12", "Rhoh", "Anxa4", "Chl1", "Alox5", "Ptpn6", "Iqgap1", 
  "Itgax", "Tspan4", "Crispld2", "Zcchc14", "Icam1", "Amigo3", "Gm29825", "Myd88", 
  "Vsir", "Pwp2", "Cbarp", "Midn", "Slc41a2", "Eif4enif1", "Trim16", "Vmp1", "Nags", 
  "Cd300lf", "Sfxn1", "Ankrd55", "Egr3", "Slc39a14", "Gm16311", "C9", "Deptor", 
  "Maff", "Vdr", "Prr13", "Gm21859", "Txndc11", "H2-Eb1", "Lrg1", "Sema6b", "C3", 
  "Ltbp1", "Atp8b1", "Npas4", "Rela", "Neat1"
)

fosl2_genes <- c(
  "Cryga", "Prlh", "Col6a3", "Fasl", "Zfp804a", "Slc20a1", "Kcns1", "Wfdc5", "A730032A03Rik",
  "Ttc39b", "C1qa", "Npffr2", "Cxcl5", "Cxcl1", "Cxcl2", "Ereg", "Gm13822", "Ddx54", "Tmem130",
  "Klrb1a", "Slc12a4", "Dpep3", "Trim43a", "Stx11", "Gm49339", "Gas7", "Mir22hg", "Myo1c", "Slfn2",
  "Slfn4", "Gip", "Hoxb9", "Gm34639", "Ankrd55", "Gdf2", "Gch1", "Gm16311", "C9", "Ext1", "Atf4",
  "Gm4719", "Hbegf", "E230025N22Rik", "Fosl1", "Ehbp1l1", "Jak2", "Cd274"
)

jun_genes <- c(
  "Ptprn", "Il1rn", "Shc4", "Rrbp1", "Fndc3b", "Fosl2", "Chrna9", "Cxcl5",
  "Cxcl1", "Cxcl2", "Ereg", "Ncor2", "Ctsc", "Ccnd1", "Srp68", "Atp13a3",
  "Gsk3b", "Smim3", "Cd74", "Tnc"
)

fos_genes <- c(
  "6720483E21Rik", "S100a6", "Anxa4", "Raet1e", "Lgals3", 
  "Egr3", "Gm21859", "Txndc11", "C3"
)

DefaultAssay(fibroblasts_final) <- "RNA"
fibroblasts_final <- AddModuleScore(fibroblasts_final, features = list(creb5_genes), name = "creb5_genes")
fibroblasts_final <- AddModuleScore(fibroblasts_final, features = list(junb_genes), name = "junb_genes")
fibroblasts_final <- AddModuleScore(fibroblasts_final, features = list(fosl2_genes), name = "fosl2_genes")
fibroblasts_final <- AddModuleScore(fibroblasts_final, features = list(jun_genes), name = "jun_genes")
fibroblasts_final <- AddModuleScore(fibroblasts_final, features = list(fos_genes), name = "fos_genes")

#mapk related gene signatures
ERK_targets <- c("Ercc1", "Gramd1b", "Vmp1", "Epn1", "Abca7", "Pkm", "Fosl2", "Actb", "Kif22", "Mfsd11", "Hook2", "Ikbkb", "Hspa8", "Sfpq", "Nr4a1", "Rbck1", "Prdx5", "Nup214", "Cd68", "Sat1", "Ldlr", "Iqcn", "Dnajb1", "Tpt1", "Slc38a2", "Cars2", "Hnrnpa1", "Flot1", "Man2c1", "Igf1r", "Samd1", "Rpl13a", "Slc20a1", "Nr4a2", "Ing1", "Eef1a1", "Eif4a2", "Ier2", "Rpl8", "Mfsd12", "Srsf2", "S100a11", "Zc3h12a", "B2m", "Tpm4", "Midn", "Ttc21a", "Mat2a", "Slc50a1", "Krt8", "Junb", "Trmt112", "Esrra", "Rbm4", "Snapc5", "Sgf29", "Rcc1", "Anxa2", "Actg1", "Pgp", "Plscr1", "S100a16", "Ppia", "Cd55", "Hspa1b", "Exosc6", "Patl2", "Gpx1", "Gas5", "Rgl2", "Neat1", "Malat1")
Trametinib_up <- c("Ccn5", "Adgrd1", "Adamtsl2", "Nkain4", "Ogn", "Wscd2", "Krt14", "Bmp3", "Gm14341", "Cemip", "4833415N18Rik", "C1qtnf3", "Asb12", "Map3k7cl", "Fgf9", "Stmn4", "Diras2", "Cntn2", "Prss35", "Ptn", "Xirp2", "Fgf16", "Gm36099", "Dpp4", "Serpinb9b", "Myom2", "Klhl30", "Ccdc141", "Fras1", "Ryr2", "Upk1b", "Cmah", "Eln", "Insc", "A330058E17Rik", "Gal3st2", "Myh1", "Chrm2", "Gm42604", "Ntf3", "Epha3", "Gm30003", "Rgs6", "Serpinb1a", "Aspn", "Wnt10b", "Adamts5", "Vegfd", "BB123696", "Dnah6", "Erbb4", "Arsi", "Slc1a6", "Serpina3h", "Plppr4", "Npy6r", "Ramp1", "Gm20560", "Atp2b4", "Serpinb6b", "Ptprz1", "Hr", "Gadd45b", "Rspo2", "Dbh", "Gm30476", "Sh3rf2", "Serpinb6c", "Eppk1", "Klhl38", "Tbxa2r", "Rasl11b", "Hmcn1", "Clec3b", "Lox", "Sh3bgr", "Nr4a3", "Gm56960", "Ppm1e", "Elobl", "Gm41724", "Gm20631", "Myo18b", "Serpina3h", "Mustn1", "Xkr4", "Upk3b", "Tph2", "Nog", "Wnt16", "Smpd3", "Gadd45g", "Gm45470", "Mir7225", "Sertm2", "Tppp", "Krt80", "Cbr2", "Tmem252", "Dpt", "Lrrc15", "Tgfb2", "Mir6903", "Rspo1", "Prelp", "Chrdl1", "Fam163b", "Slc24a3", "Gm29491", "Sema3b", "Hey2", "Snta1", "Adap1", "Gm35546", "Tmem196", "Myh2", "Prkaa2", "Gm10824", "Kcnj6", "Crabp2", "Nat8l", "Gm35853", "Cd248", "Krt7", "Fyb2", "Rtl3", "Pgm5", "Cnn1", "Kcnn2", "Avpr1a", "Tlr5", "Cthrc1", "Ror1", "Zfp185", "Mrgprf", "Alkal1", "Gm41505", "Igfbp2", "Dusp8", "Itga7", "Podn", "Gm15851", "Mcf2l", "Fmod", "Fam124a", "Bco1", "Eno3", "Gabrb2", "Adcy1", "Rflnb", "Tek", "Tcp11x2", "Id4", "Gm14066", "Muc16", "Pdlim3", "Pkia", "A730020E08Rik", "Prss23os", "Itgb1bp2", "Gm20744", "Prss23", "Lrrn4cl", "Abca6", "Plin4", "Col19a1", "Gm49291", "Cdh8", "Ryr3", "Fbln5", "Sync", "Igfbp5", "Susd2", "Ackr3", "Nppb", "Hspb7", "Gm19938", "Trim29", "Mmp17", "Plxna4os3", "Actg2", "Ism1", "Gm56664", "Kctd16", "A330033J07Rik", "Parm1", "Chrd", "Plxna4", "Gm57353", "Smim10l2a", "Siglecg", "Sh2d4a", "Dipk2a", "Anxa8", "Ccdc68", "Srpk3", "Lims2", "Tgm2", "Fzd6", "Acta2")

DefaultAssay(fibroblasts_final) <- "RNA"
fibroblasts_final <- AddModuleScore(fibroblasts_final, features = list(ERK_targets), name = "ERK_targets")
fibroblasts_final <- AddModuleScore(fibroblasts_final, features = list(Trametinib_up), name = "Trametinib_up")


# ── 5. Pseudobulking multiome dataset ─────────────────────────────────────────

#pseudobulking RNA expression and gene signature module scores
#choose pseudobulking variables
group_vars   <- c("sampleID", "cellType_specific")  #grouping variables
min_cells    <- 30  #minimum cells per group
module_scores <- c(
  "Curras_fibrosis_sig1",
  "Tsukui_fibrosis_sig1",
  "all_mfb_ecm_genes1",
  "matrix_genes_final1",
  "creb5_genes1","junb_genes1","fos_genes1","jun_genes1","fosl2_genes1" #grn programs
)

DefaultAssay(fibroblasts_final) <- "RNA"

#create combined group label and filter small groups
fibroblasts_final$pb_group <- paste(
  fibroblasts_final@meta.data[[group_vars[1]]],
  fibroblasts_final@meta.data[[group_vars[2]]],
  sep = "_"
)

group_counts <- table(fibroblasts_final$pb_group)
keep_groups  <- names(group_counts[group_counts >= min_cells])

#subset object
obj_filtered <- subset(fibroblasts_final, pb_group %in% keep_groups)

#aggregate gene expression
pb_rna <- AggregateExpression(obj_filtered, assays = "RNA", group.by = "pb_group", return.seurat = FALSE)$RNA

#pseudobulk module scores
module_mat <- obj_filtered@meta.data %>%
  select(all_of(c("pb_group", module_scores))) %>%
  group_by(pb_group) %>%
  summarise(across(all_of(module_scores), mean, na.rm = TRUE), .groups = "drop")

#create metadata
group_meta <- obj_filtered@meta.data %>%
  select(pb_group, all_of(group_vars), condition_broad) %>%
  distinct() %>%
  left_join(
    obj_filtered@meta.data %>% count(pb_group, name = "n_cells"),
    by = "pb_group"
  ) %>%
  arrange(pb_group) %>%
  left_join(module_mat, by = "pb_group")

#align RNA columns to match group_meta row order
pb_rna <- pb_rna[, group_meta$pb_group]


# ── 6. GRN analysis ───────────────────────────────────────────────────────────

#running GRaNIE
output_dir      <- "/path/to/dir"
genome_assembly <- "mm10"
organism_name   <- "Mus musculus"

#pseudobulk RNA + ATAC counts by cell type
rna_matrix_pb  <- AggregateExpression(fibroblasts_final, assays = "RNA",  slot = "counts",
                                      group.by = "cellType_specific", return.seurat = FALSE)$RNA
atac_matrix_pb <- AggregateExpression(fibroblasts_final, assays = "ATAC", slot = "counts",
                                      group.by = "cellType_specific", return.seurat = FALSE)$ATAC
if(!is.matrix(rna_matrix_pb))  rna_matrix_pb  <- as.matrix(rna_matrix_pb)
if(!is.matrix(atac_matrix_pb)) atac_matrix_pb <- as.matrix(atac_matrix_pb)

#map RNA rownames: gene symbol to ensembl id
original_symbols <- rownames(rna_matrix_pb)
ensembl_ids <- mapIds(org.Mm.eg.db, keys = original_symbols,
                      column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
mapped_idx          <- !is.na(ensembl_ids)
rna_matrix_pb_mapped<- rna_matrix_pb[mapped_idx, , drop = FALSE]
ensembl_ids_mapped  <- ensembl_ids[mapped_idx]

if (any(duplicated(ensembl_ids_mapped))) {
  # Sum counts when multiple symbols map to the same Ensembl ID
  rna_matrix_pb_ensembl <- rowsum(as.matrix(rna_matrix_pb_mapped),
                                  group = ensembl_ids_mapped, reorder = TRUE, na.rm = TRUE)
} else {
  rownames(rna_matrix_pb_mapped) <- ensembl_ids_mapped
  rna_matrix_pb_ensembl <- rna_matrix_pb_mapped
}

#build data frames for addData
idColumn_RNA_name   <- "EnsemblID"
rna_df_for_granie   <- data.frame(rownames(rna_matrix_pb_ensembl), rna_matrix_pb_ensembl,
                                  check.names = FALSE)
colnames(rna_df_for_granie)[1] <- idColumn_RNA_name
rownames(rna_df_for_granie)    <- NULL

idColumn_peaks_name <- "peakID"
atac_df_for_granie  <- data.frame(rownames(atac_matrix_pb), atac_matrix_pb, check.names = FALSE)
colnames(atac_df_for_granie)[1] <- idColumn_peaks_name
rownames(atac_df_for_granie)    <- NULL

cluster_ids        <- colnames(rna_matrix_pb_ensembl)
sampleMetadata_df  <- data.frame(SampleID = cluster_ids, row.names = cluster_ids)

#initialize GRN and add data
objectMetadata.l <- list(project_name   = "GRaNIE_Mouse_Multiome",
                         organism       = organism_name,
                         genomeAssembly = genome_assembly)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
grn_object <- initializeGRN(objectMetadata = objectMetadata.l,
                            outputFolder   = output_dir,
                            genomeAssembly = genome_assembly)

grn_object <- addData(grn_object,
                      counts_peaks = atac_df_for_granie, normalization_peaks = "DESeq2_sizeFactors",
                      idColumn_peaks = idColumn_peaks_name,
                      counts_rna   = rna_df_for_granie,  normalization_rna   = "limma_quantile",
                      idColumn_RNA = idColumn_RNA_name,
                      sampleMetadata = sampleMetadata_df, forceRerun = TRUE)

#TFBS and connections
grn_object <- addTFBS(grn_object,
                      motifFolder = "/path/to/PWMScan_HOCOMOCOv12/H12INVIVO/pwmscan_filt",
                      TFs = "all", filesTFBSPattern = "_TFBS", fileEnding = ".bed.gz",
                      translationTable = "translationTable_mm10.csv", forceRerun = TRUE)
grn_object <- overlapPeaksAndTFBS(grn_object, nCores = 1, forceRerun = TRUE)
grn_object <- addConnections_TF_peak(grn_object, corMethod = "pearson", forceRerun = TRUE)
grn_object <- addConnections_peak_gene(grn_object, promoterRange = 250000,
                                       corMethod = "pearson", forceRerun = TRUE)

#filter and finalize
grn_object <- filterGRNAndConnectGenes(grn_object,
                                       TF_peak.fdr.threshold = 0.3, peak_gene.fdr.threshold = 0.3,
                                       peak_gene.fdr.method = "BH",
                                       gene.types = c("protein_coding", "lincRNA"),
                                       allowMissingTFs = FALSE, allowMissingGenes = FALSE, forceRerun = TRUE)
grn_object <- add_TF_gene_correlation(grn_object, corMethod = "pearson",
                                      nCores = 1, forceRerun = TRUE)

GRN_connections.all <- getGRNConnections(grn_object, type = "all.filtered",
                                         include_TF_gene_correlations = TRUE,
                                         include_geneMetadata = TRUE)

# ── 7. Save objects ───────────────────────────────────────────────────────────

saveRDS(fibroblasts_final_chrom,         "/path/to/fibroblasts_final_chrom.rds")
saveRDS(multiome_integrated_downsampled, "/path/to/multiome_integrated_downsampled.rds")
saveRDS(grn_object,                      "/path/to/grn_object.rds")
