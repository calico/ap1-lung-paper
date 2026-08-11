# Analysis of the published whole-lung microarray dataset GSE47460 (LGRC)
# ssGSEA scoring of the fibrotic and developmental (AFOS) ECM signatures across
# disease groups, and their association with lung function.
# Input: GSE47460 is downloaded directly from GEO in this script.

# ── 0. Load Packages ──────────────────────────────────────────────────────

library(GEOquery)
library(Biobase)
library(limma)
library(GSVA)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

# ── 1. GSE47460 microarray dataset analysis ────────────────────────────

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
