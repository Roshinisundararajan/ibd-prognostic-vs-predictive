# =====================================================================
# 00_setup.R — run this first. Installs are the slowest part of Day 1,
# so start it and read the runbook while it works.
#
# The list is deliberately short. pROC, PRROC and sva were in the first
# draft and have been removed: the pipeline now computes AUROC, PR-AUC
# and the ROC curve itself (see _helpers.R), and it never runs ComBat.
# Fewer Bioconductor compiles, fewer ways for Day 1 to be eaten alive.
# =====================================================================

cran <- c("glmnet",      # elastic net
          "ranger",      # random forest
          "data.table",  # fast IO
          "matrixStats", # fast column statistics inside CV folds
          "ggplot2",     # figures
          "logistf")     # Firth-penalised logistic regression, for the
                         # interaction test at ~60 events

bioc <- c("GEOquery",    # GEO download
          "Biobase",     # ExpressionSet accessors
          "singscore")   # rank-based single-sample scoring (optional:
                         # script 06 falls back to a built-in equivalent)

new <- setdiff(cran, rownames(installed.packages()))
if (length(new)) install.packages(new, repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
new <- setdiff(bioc, rownames(installed.packages()))
if (length(new)) BiocManager::install(new, ask = FALSE, update = FALSE)

for (d in c("out", "data", "data/signatures"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

cat("\n--- availability check ---\n")
for (p in c(cran, bioc))
  cat(sprintf("  %-14s %s\n", p, requireNamespace(p, quietly = TRUE)))
cat("\nsingscore FALSE is survivable; logistf FALSE is not — script 06\n")
cat("will fall back to unpenalised glm and the odds ratios will blow up.\n\n")
sessionInfo()
