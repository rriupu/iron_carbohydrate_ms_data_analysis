library(DEP)
library(dplyr)
library(ggplot2)

# https://www.bioconductor.org/packages/devel/bioc/vignettes/DEP/inst/doc/DEP.html

## ---------- ##
## Proteomics ##
## ---------- ##

output_dir = "~/Documents/proteomics/proteomics/"
dir.create(output_dir, recursive = T, showWarnings = F)

# --------------- #
# Read input file #
# --------------- #
path = "~/Documents/proteomics/proteomics/proteinGroups.txt"

df = read.table(path, sep = "\t", header = T)

# ---------------------------------------------------------------- #
# Filtering and preparing the file to build a SummarizedExperiment #
# ---------------------------------------------------------------- #

df = df %>% filter(Reverse != "+", Potential.contaminant != "+")

df_unique = make_unique(df, "Protein.IDs", "Protein.IDs", delim = ";")

# name column should have gene symbol, ID column should have UniProt
# some entries have zz|Y-FGCZContXXXXX| (X are numbers)
parts = strsplit(df_unique$name, "|", fixed = T)
names_col = sapply(
  parts,
  function(x) {
    if (x[1] == "zz") {
      return(x[2])
    }
    if (x[1] == "sp") {
      return(x[3])
    }
  }
)
names_col = gsub("_HUMAN", "", names_col)
IDs_col = sapply(parts, "[[", 2)

df_unique$name = names_col
df_unique$ID = IDs_col

# --------------------------------------------------- #
# Prepare other objects to build SummarizedExperiment #
# --------------------------------------------------- #

LFQ_columns = grep("LFQ.", colnames(df_unique))

colnames_of_interest = colnames(df_unique)[LFQ_columns][!grepl("Monofer", colnames(df_unique)[LFQ_columns])]

LFQ_columns = which(colnames(df_unique) %in% colnames_of_interest)

conditions = sapply(strsplit(colnames_of_interest, "_"), "[[", 7)
donors = sapply(strsplit(colnames_of_interest, "_"), "[[", 8)
donors = as.numeric(gsub("R", "", donors))

experimental_design = data.frame(
    label = colnames_of_interest,
    condition = conditions,
    replicate = donors)

# -------------------------- #
# Build SummarizedExperiment #
# -------------------------- #

df_se = make_se(df_unique, LFQ_columns, experimental_design)

# ----------------------------------------------------------------- #
# Remove proteins for which there are no observations in any sample #
# ----------------------------------------------------------------- #

exploration_plots_dir = file.path(output_dir, "data_exploration")
dir.create(exploration_plots_dir, recursive = T, showWarnings = F)

freq_plot = plot_frequency(df_se)
ggsave(file.path(exploration_plots_dir, "protein_frequency_before_filtering.pdf"), freq_plot)

df_filt = filter_missval(df_se, thr = 0)
freq_plot_after_filtering = plot_frequency(df_filt)
ggsave(file.path(exploration_plots_dir, "protein_frequency_after_filtering.pdf"), freq_plot_after_filtering)

proteins_per_sample_plot = plot_numbers(df_filt)
ggsave(file.path(exploration_plots_dir, "proteins_per_sample_and_condition.pdf"), proteins_per_sample_plot)

protein_coverage_plot = plot_coverage(df_filt)
ggsave(file.path(exploration_plots_dir, "protein_coverage.pdf"), protein_coverage_plot)

# ------------------ #
# Data normalization #
# ------------------ #

df_norm = normalize_vsn(df_filt)

# ---------------- #
# Data exploration #
# ---------------- #

pca_plot = plot_pca(df_norm, x = 1, y = 2, n = 500, point_size = 4)
ggsave(file.path(exploration_plots_dir, "pca_no_correction.pdf"), pca_plot)

# See if donor effect can be removed
df_norm_nobatch = df_norm
assay(df_norm_nobatch) = limma::removeBatchEffect(
  assay(df_norm_nobatch),
  df_norm_nobatch$replicate)

pca_plot_nobatch = plot_pca(df_norm_nobatch, x = 1, y = 2, n = 500, point_size = 4)
ggsave(file.path(exploration_plots_dir, "pca_corrected.pdf"), pca_plot_nobatch)

mean_sd_plot = meanSdPlot(df_norm)
pdf(file.path(exploration_plots_dir, "mean_sd_plot.pdf"))
print(mean_sd_plot)
dev.off()
normalization_plot = plot_normalization(df_filt, df_norm)
ggsave(file.path(exploration_plots_dir, "normalization.pdf"), normalization_plot)

missval_plot = plot_missval(df_filt)
pdf(file.path(exploration_plots_dir, "missval.pdf"))
print(missval_plot)
dev.off()
miss_val_log2_plot = plot_detect(df_filt)
ggsave(file.path(exploration_plots_dir, "missval_log2_distribution.pdf"), miss_val_log2_plot)

# --------------- #
# Data imputation #
# --------------- #

dea_dir = file.path(output_dir, "dea")
dir.create(dea_dir, recursive = T, showWarnings = F)

df_imp_man = impute(df_norm, fun = "man", shift = 1.4, scale = 0.4)
imputation_plot = plot_imputation(df_norm, df_imp_man)
ggsave(file.path(dea_dir, "imputation.pdf"), imputation_plot)

# --------------------- #
# Differential analysis #
# --------------------- #

df_diff = test_diff(
  df_imp_man,
  type = "manual",
  test = c("Ferinject_vs_Venofer", "Venofer_vs_UT", "Ferinject_vs_UT"),
  design_formula = formula(~ 0 + condition + replicate))

dep = add_rejections(df_diff, alpha = 0.05, lfc = log2(1.5))

# ---------------- #
# Plotting results #
# ---------------- #

corplot = plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")
pdf(file.path(exploration_plots_dir, "correlation_plot.pdf"))
print(corplot)
dev.off()

htmp_centered = plot_heatmap(
  dep,
  type = "centered",
  kmeans = TRUE,
  k = 4,
  col_limit = 4,
  show_row_names = T,
  indicate = c("condition", "replicate"))

pdf(file.path(exploration_plots_dir, "heatmap_centered.pdf"))
print(htmp_centered)
dev.off()

htmp_contrast = plot_heatmap(
  dep,
  type = "contrast",
  kmeans = TRUE,
  k = 3,
  col_limit = 10,
  show_row_names = TRUE)

pdf(file.path(dea_dir, "heatmap_contrast.pdf"))
print(htmp_contrast)
dev.off()

volcano_fer_veno = plot_volcano(
  dep,
  contrast = "Ferinject_vs_Venofer",
  label_size = 2,
  add_names = TRUE) + 
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")

ggsave(file.path(dea_dir, "volcano_ferinject_vs_venofer.pdf"), volcano_fer_veno)

volcano_veno_ut = plot_volcano(
  dep,
  contrast = "Venofer_vs_UT",
  label_size = 2,
  add_names = TRUE) + 
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")

ggsave(file.path(dea_dir, "volcano_venofer_vs_untreated.pdf"), volcano_veno_ut)

volcano_fer_ut = plot_volcano(
  dep,
  contrast = "Ferinject_vs_UT",
  label_size = 2,
  add_names = TRUE) + 
    geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed") + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed")
    
ggsave(file.path(dea_dir, "volcano_ferinject_vs_untreated.pdf"), volcano_fer_ut)

diff_exp_prots = as.data.frame(rowData(dep)) %>% filter(significant) %>% pull(name)

single_gene_plots = plot_single(dep, proteins = diff_exp_prots)
ggsave(file.path(dea_dir, "single_gene_plots.pdf"), single_gene_plots)

condition_overlap_plot = plot_cond(dep)
ggsave(file.path(dea_dir, "condition_overlap.pdf"), condition_overlap_plot)

df_results = get_results(dep)

write.table(
  df_results,
  file.path(dea_dir, "results.tsv"),
  col.names = T,
  row.names = F,
  sep = "\t",
  quote = F)
  
veno_ut_pvals = rowData(dep) %>%
  as.data.frame() %>%
  pull(Venofer_vs_UT_p.val)

veno_ut_pvals_hist = ggplot(data.frame(pvals = veno_ut_pvals), aes(x = pvals)) + geom_histogram(bins = 100)
ggsave(file.path(dea_dir, "venofer_vs_ut_pval_hist.pdf"), veno_ut_pvals_hist)

fer_ut_pvals = rowData(dep) %>%
  as.data.frame() %>%
  pull(Ferinject_vs_UT_p.val)

fer_ut_pvals_hist = ggplot(data.frame(pvals = fer_ut_pvals), aes(x = pvals)) + geom_histogram(bins = 100)
ggsave(file.path(dea_dir, "ferinject_vs_ut_pval_hist.pdf"), fer_ut_pvals_hist)

fer_veno_pvals = rowData(dep) %>%
  as.data.frame() %>%
  pull(Ferinject_vs_Venofer_p.val)

fer_veno_pvals_hist = ggplot(data.frame(pvals = fer_veno_pvals), aes(x = pvals)) + geom_histogram(bins = 100)
ggsave(file.path(dea_dir, "ferinject_vs_venofer_pval_hist.pdf"), fer_veno_pvals_hist)

RAB_genes = as.data.frame(rowData(dep)) %>% filter(grepl("RAB", rowData(dep)$name)) %>% pull(name)
rab_genes = plot_single(dep, proteins = RAB_genes)
ggsave(file.path(dea_dir, "RAB_genes.pdf"), rab_genes)

lamp1 = plot_single(dep, proteins = "LAMP1")
ggsave(file.path(dea_dir, "lamp1.pdf"), lamp1)