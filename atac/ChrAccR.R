# ATAC-seq analysis with ChrAccR
# End-to-end ATAC-seq processing (peak calling, normalization, differential accessibility)
# for EGFP- vs AFOS-transduced mouse pulmonary fibroblasts
# Inputs: sampleAnnotation.tsv  sample metadata
#         BAM_files/            aligned BAM files named in the bamFilename column produced from atac_process.sh

# ── 0. Load Packages ──────────────────────────────────────────────────────────

library(ChrAccR)
library(muRtools)
library(ggplot2)
library(chromVAR)        
library(chromVARmotifs) 

# ── 1. Set configs ────────────────────────────────────────────────────────────

theme_set(muRtools::theme_nogrid())
sampleAnnot     <- file.path("/path/to/directory", "sampleAnnotation.tsv")
bamDir          <- file.path("/path/to/directory", "BAM_files")
sampleAnnotATAC <- read.table(sampleAnnot, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#add a column that ChrAccR can use to find the correct bam file for each sample
sampleAnnotATAC[, "bamFilenameFull"] <- file.path(bamDir, sampleAnnotATAC[, "bamFilename"])

setConfigElement("annotationColumns", c("condition"))
setConfigElement("colorSchemes", c(getConfigElement("colorSchemes"),
  list(
    "condition" = c("EGFP" = "#78F542", "AFOS" = "#F5426C")
  )
))
setConfigElement("filteringCovgCount", 1L)
setConfigElement("filteringCovgReqSamples", 0.25)
setConfigElement("filteringSexChroms", TRUE)
setConfigElement("normalizationMethod", "quantile")
diffCompNames <- c(
  "EGFP vs AFOS [condition]"
)
setConfigElement("differentialCompNames", diffCompNames)
setConfigElement("regionTypes", c("promoter", ".peaks.cons"))
setConfigElement("doPeakCalling", TRUE)
setConfigElement("annotationPeakGroupColumn", "condition")
setConfigElement("annotationPeakGroupAgreePerc", 0.5)

# ── 2. Run ChrAccR ────────────────────────────────────────────────────────────

run_atac("/path/to/output", "bamFilenameFull", sampleAnnotATAC,
         genome = "mm10", sampleIdCol = "sampleID")
