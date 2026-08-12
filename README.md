# Code for *AP-1 specifies developmental versus fibrotic extracellular matrix transcriptional programs in the lung*

Analysis code accompanying **AM Kaiser, A Luo, K Sivasubramanian, G Dueñas, B Martin-McNulty, W Cedron-Craft, JM Patino, AJ Chang, J Riegler, & JG Ruby**, *AP-1 specifies developmental versus fibrotic extracellular matrix transcriptional programs in the lung*.

This repository contains the R scripts used to process and analyze the single-cell multiome, bulk RNA-seq, ATAC-seq, and published human datasets reported in the paper, and to generate the figures.

--------------

## Instructions

### Clone git repository

```
git clone https://github.com/<GITHUB_USERNAME>/ap1-lung-paper.git
```

### Data and file paths

The scripts are not organized as an R package; each is a standalone script meant to be run section by section in an R session (or with `Rscript`). Nothing is installed from this repository — the packages the scripts depend on are installed from CRAN, Bioconductor, and GitHub. The user needs to download raw and processed data from online repositories to run the analyses in this repository. Inputs for each script are provided in data deposited at the GEO and are listed at the top of each script. The scripts use placeholder paths (e.g. `/path/to/...`). To reproduce an analysis, set these to the location of the corresponding data on your system. 

Raw and processed data are deposited at:

- This study: {GEO accession number}
- Published datasets reanalyzed here: Habermann et al: GSE135893, Adams et al: GSE136831, Tsukui et al: GSE132771, Valenzi et al: GSE214085, Zepp et al: GSE149563, Strunz et al: GSE141259, Curras-Alonso et al: GSE211713, Narvaez del Pilar et al: GSE180822, LGRC: GSE47460 — see each dataset's original publication for terms of use.

### R environments

The analyses were run in two R environments. `package_versions.txt` lists which environment each script was run in, and the version of every package used.

--------------

## Repository structure

- **`multiome/`** — mouse lung single-nuclei multiome (RNA + ATAC)
  - `initial_processing_01.R` — load 10x outputs, QC, integrate, write the base Seurat objects
  - `analysis_02.R` — downstream analysis, chromatin scoring, GRN inference
  - `external_sc_dataset_analysis_03.R` — reprocessing of published external datasets for comparison
  - `figures_04.R` — reads the objects above and produces the multiome figure panels
- **`human/`** — reanalysis of published human lung datasets
  - `human_sc_analysis_01.R`, `human_sc_figures_02.R` — analysis of published single-cell RNA- and ATAC-seq datasets and associated figures
  - `human_microarray_analysis_03.R` - whole-lung microarray dataset analysis
- **`bulk_rnaseq/`** — bulk RNA-seq analyses:
  - `AFOS.R`, `KD_OE.R`, `inhibitor_screen.R`
- **`atac/`** — bulk ATAC-seq accessibility analysis (`ChrAccR.R`).

### Dependencies between scripts

Most analyses are contained within their experimental folder. However, a few analyses depend on objects produced elsewhere in the repo: `bulk_rnaseq/KD_OE.R` reads the fibroblasts_fina RDS object in the GEO repository, and `bulk_rnaseq/inhibitor_screen.R` uses the AFOS sample table and counts matrix.
