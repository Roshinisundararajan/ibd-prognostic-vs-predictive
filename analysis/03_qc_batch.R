# =====================================================================
# 03_qc_batch.R
#
# Quick QC. Budget 15 minutes, produce two figures, move on.
#
# Batch correction policy: do NOT ComBat the discovery cohort. One
# platform, one trial — correcting within it removes signal and buys
# nothing. The commented block at the bottom is for a cross-cohort
# PICTURE only. Never fit a model on ComBat output that mixes cohorts,
# and never pass the outcome into `mod`, which leaks it into the data.
#
# IN   out/02_cohort_<GSE>.rds
# OUT  out/03_pca.png, out/03_qc_boxplot.png, out/03_qc_notes.txt
# =====================================================================

suppressPackageStartupMessages({ library(ggplot2); library(matrixStats) })
source("_helpers.R")
set.seed(20261017)

GSE_ID <- "GSE206285"
co <- readRDS(file.path("out", paste0("02_cohort_", GSE_ID, ".rds")))
ex <- co$expr
notes <- character(0)
note <- function(...) { s <- paste0(...); cat(s, "\n"); notes <<- c(notes, s) }

## ---- per-sample intensity distributions ------------------------------
png("out/03_qc_boxplot.png", width = 1400, height = 700)
idx <- sample(ncol(ex), min(60, ncol(ex)))
boxplot(ex[, idx], outline = FALSE, las = 2, cex.axis = 0.6,
        main = sprintf("Per-sample intensity, %d random samples", length(idx)),
        ylab = "log2 expression")
dev.off()

## ---- array-level outliers, flagged not dropped -----------------------
# Correlation of each array with the cohort median profile. Arrays far
# below the rest are flagged in the report; they are not removed, because
# dropping samples on a QC metric this crude costs more than it saves.
medp <- rowMedians(ex)
r    <- apply(ex, 2, function(z) cor(z, medp, method = "spearman"))
z    <- (r - median(r)) / mad(r)
note(sprintf("median array-vs-cohort correlation: %.3f", median(r)))
note(sprintf("arrays with robust z < -3.5: %d  (%s)", sum(z < -3.5),
             paste(head(names(which(z < -3.5)), 10), collapse = " ")))

## ---- PCA on the most variable genes ----------------------------------
top <- head(order(rowVars(ex), decreasing = TRUE), 2000)
pca <- prcomp(t(ex[top, ]), scale. = TRUE)
pv  <- 100 * pca$sdev^2 / sum(pca$sdev^2)

df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                 arm = co$arm,
                 outcome = factor(co$y, levels = c(0, 1),
                                  labels = c("no remission", "remission")))
p <- ggplot(df, aes(PC1, PC2, colour = outcome, shape = arm)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = sprintf("PC1 (%.1f%%)", pv[1]), y = sprintf("PC2 (%.1f%%)", pv[2]),
       title = "Baseline mucosal transcriptome, discovery cohort") +
  theme_minimal(base_size = 13)
ggsave("out/03_pca.png", p, width = 7, height = 5, dpi = 200)

note(sprintf("PC1 variance explained: %.1f%%   PC2: %.1f%%", pv[1], pv[2]))

## ---- two sanity checks worth putting in the report -------------------
# 1. The outcome should NOT separate on PC1/PC2. If it does, suspect a
#    label or batch artefact before celebrating a finding.
note(sprintf("PC1 ~ outcome: p = %.3g",
             t.test(PC1 ~ outcome, data = df)$p.value))
note(sprintf("PC2 ~ outcome: p = %.3g",
             t.test(PC2 ~ outcome, data = df)$p.value))

# 2. The arm should not separate either — these are baseline, pre-treatment
#    biopsies from a randomised trial, so an arm difference here would mean
#    randomisation or processing went wrong. A null result is the pass.
if (length(unique(na.omit(df$arm))) == 2)
  note(sprintf("PC1 ~ arm (randomisation check, expect null): p = %.3g",
               t.test(PC1 ~ arm, data = df)$p.value))

# 3. Technical processing codes, if script 01 recovered them.
if (!is.null(co$clin) && "b_code" %in% names(co$clin))
  note("b_code present in clinical block; test association with PC1 if needed")

writeLines(notes, "out/03_qc_notes.txt")
cat("\nwrote out/03_pca.png, out/03_qc_boxplot.png, out/03_qc_notes.txt\n")
cat("NEXT: 04_models_cv.R\n")

# ---- ComBat, for a cross-cohort PICTURE only. Do not model on this. ----
# library(sva)
# common    <- intersect(rownames(ex_disc), rownames(ex_repl))
# combined  <- cbind(ex_disc[common, ], ex_repl[common, ])
# batch     <- c(rep("discovery", ncol(ex_disc)), rep("replication", ncol(ex_repl)))
# corrected <- ComBat(dat = combined, batch = batch, mod = NULL)   # no outcome
