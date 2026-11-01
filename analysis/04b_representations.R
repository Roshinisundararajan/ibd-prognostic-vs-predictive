install.packages("remotes")
remotes::install_github("dviraran/xCell")

# =====================================================================
# 04b_representations.R — which compression of the transcriptome carries
# the signal? Runs BETWEEN 03 and 04.
#
# Twenty thousand genes against sixty-odd events is hopeless, so every
# model needs the genes compressed somehow. This script compares five
# ways of doing that under one cross-validation scheme.
#
# THE SELECTION TRAP, READ THIS BEFORE YOU USE THE OUTPUT
# -------------------------------------------------------
# Picking the best of five here and then reporting its AUROC as the
# headline is cheating, quietly: you chose it using the same ~63 events
# you then report it on. At five options the inflation is not small, and
# it is the same error that `enet_leaky` in script 04 exists to expose.
#
# So PATHWAY SCORES ARE PRE-SPECIFIED AS PRIMARY, before you look at
# anything. Fewest features, most interpretable, least able to overfit at
# this event count. The other four are secondary and are reported as a
# grid with intervals, not ranked into a winner. If you want a different
# primary, change PRIMARY below and write down the date you changed it,
# BEFORE running this.
#
# IN   out/02_cohort_<GSE>.rds
#      data/genesets/*.gmt          (optional) pathway collection
#      out/06_signature_scores.csv  (optional) published signature scores
# OUT  out/04b_representation_grid.csv
#      out/04b_representation_preds.csv
# =====================================================================

suppressPackageStartupMessages({
  library(glmnet); library(data.table); library(matrixStats)
})
source("_helpers.R")
set.seed(20261017)

GSE_ID    <- "GSE206285"
PRIMARY   <- "pathway"        # pre-specified. Do not change after looking.
N_FOLDS   <- 5
N_REPEATS <- 3
TOP_K     <- 2000
ALPHA     <- 0.5
N_PC      <- 20

co   <- readRDS(file.path("out", paste0("02_cohort_", GSE_ID, ".rds")))
y    <- co$y
ids  <- co$ids
exf  <- t(co$expr)                                   # variance-filtered, samples x genes
exa  <- if (!is.null(co$expr_all)) co$expr_all else co$expr   # unfiltered, genes x samples

cat(sprintf("%d patients, %d events (%.1f%%)\n", length(y), sum(y), 100 * mean(y)))

## =====================================================================
## Build the representations
## =====================================================================
# Everything here is OUTCOME-BLIND, so it can be computed once on all
# samples. Rank-based scores in particular use only a patient's own genes
# and nothing about other patients' labels. The steps that DO see the
# outcome — the gene filter, the PCA basis, scaling, elastic-net selection
# — are refit per fold further down.
reps <- list()

## 1. clinical only ------------------------------------------------------
clin_mm <- build_clin_mm(co$clin)
if (ncol(clin_mm)) reps$clinical <- clin_mm

## 2. top-variance genes -------------------------------------------------
reps$topvar <- exf

## 3. principal components ----------------------------------------------
# Placeholder only; the basis is refit inside each fold below. Carrying
# the matrix here keeps the loop uniform.
reps$pca <- exf

## 4. pathway scores — THE PRE-SPECIFIED PRIMARY -------------------------
# Needs a gene-set collection. Put Hallmark or Reactome .gmt files in
# data/genesets/. Without them this falls back to the published signature
# files, which is a much smaller collection — say which you used.
read_gmt <- function(path) {
  ln <- strsplit(readLines(path, warn = FALSE), "\t")
  ln <- ln[lengths(ln) > 3]
  setNames(lapply(ln, function(z) toupper(trimws(z[-c(1, 2)]))),
           vapply(ln, function(z) z[1], character(1)))
}
rank_score <- function(ex, up, dn = character(0)) {
  R  <- apply(ex, 2, function(z) rank(z, ties.method = "average") / length(z))
  su <- if (length(up)) colMeans(R[up, , drop = FALSE]) else rep(0.5, ncol(ex))
  if (length(dn)) su - colMeans(R[dn, , drop = FALSE]) else su - 0.5
}

gmts <- list.files("data/genesets", pattern = "\\.gmt$", full.names = TRUE)
sets <- list()
if (length(gmts)) {
  for (g in gmts) sets <- c(sets, read_gmt(g))
  cat(sprintf("gene sets: %d from %d GMT file(s)\n", length(sets), length(gmts)))
} else if (length(list.files("data/signatures", pattern = "\\.txt$"))) {
  for (f in list.files("data/signatures", pattern = "\\.txt$", full.names = TRUE))
    sets[[sub("\\.txt$", "", basename(f))]] <- toupper(trimws(readLines(f, warn = FALSE)))
  cat(sprintf("no GMT files — falling back to %d published signature lists.\n", length(sets)))
  cat("State this in the methods; it is a much smaller collection than Hallmark.\n")
} else {
  cat("*** no gene sets and no signature files: the PRIMARY representation\n")
  cat("    cannot be built. Put a .gmt in data/genesets/ and rerun.\n")
}

if (length(sets)) {
  sets <- lapply(sets, function(g) intersect(g, rownames(exa)))
  sets <- sets[lengths(sets) >= 5]
  cat(sprintf("usable gene sets (>= 5 genes present): %d\n", length(sets)))
  if (length(sets)) {
    use_ss <- requireNamespace("singscore", quietly = TRUE)
    if (use_ss) {
      rk <- singscore::rankGenes(exa)
      P  <- vapply(sets, function(g)
             as.numeric(singscore::simpleScore(rk, upSet = g)$TotalScore),
             numeric(ncol(exa)))
    } else {
      cat("singscore not installed — using the built-in rank score.\n")
      P <- vapply(sets, function(g) rank_score(exa, g), numeric(ncol(exa)))
    }
    rownames(P) <- colnames(exa)
    reps$pathway <- P
  }
}

## 5. cell-type deconvolution -------------------------------------------
# A biopsy's response is substantially a cell-composition story, so this
# is the representation most likely to beat the others on biology rather
# than on variance. xCell MUST see the unfiltered matrix.
if (requireNamespace("xCell", quietly = TRUE)) {
  X <- xCell::xCellAnalysis(as.data.frame(exa))
  reps$deconv <- t(X)
  cat(sprintf("xCell: %d cell types\n", nrow(X)))
} else {
  cat("xCell not installed — the deconvolution representation is skipped.\n")
  cat("  BiocManager::install('xCell')  or  remotes::install_github('dviraran/xCell')\n")
}

## 6. published signature scores, if script 06 has already run ----------
sig_path <- "out/06_signature_scores.csv"
if (file.exists(sig_path)) {
  S <- fread(sig_path, data.table = FALSE)
  rownames(S) <- S$sample
  S <- S[ids, setdiff(names(S), c("sample", "y", "arm")), drop = FALSE]
  if (ncol(S)) reps$signatures <- as.matrix(S)
}

for (nm in names(reps)) stopifnot(nrow(reps[[nm]]) == length(y))
cat("\nrepresentations built: ", paste(sprintf("%s (%d features)", names(reps),
    vapply(reps, ncol, integer(1))), collapse = " · "), "\n", sep = "")
if (!PRIMARY %in% names(reps))
  cat("\n*** the pre-specified primary ('", PRIMARY, "') could not be built.\n",
      "    Do not silently promote another one. Fix the inputs, or record the\n",
      "    change of primary with its date before you look at the grid.\n", sep = "")

## =====================================================================
## One cross-validation scheme, applied to every representation
## =====================================================================
strat_folds <- function(y, k) {
  f <- integer(length(y))
  for (cl in unique(y)) {
    idx <- which(y == cl)
    f[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  f
}

preds <- list(); t0 <- Sys.time()
for (rep_i in seq_len(N_REPEATS)) {
  folds <- strat_folds(y, N_FOLDS)
  for (k in seq_len(N_FOLDS)) {
    tr <- which(folds != k); te <- which(folds == k)
    ytr <- y[tr]
    row <- data.table(rep = rep_i, fold = k, sample = ids[te], y = y[te])

    for (nm in names(reps)) {
      M <- reps[[nm]]

      # impute, then scale, on training rows only
      med <- impute_fit(M[tr, , drop = FALSE])
      Mtr <- impute_apply(M[tr, , drop = FALSE], med)
      Mte <- impute_apply(M[te, , drop = FALSE], med)

      if (nm %in% c("topvar", "pca")) {
        tst <- abs(col_t(Mtr, ytr))                      # in-fold gene filter
        sel <- names(sort(tst, decreasing = TRUE))[seq_len(min(TOP_K, length(tst)))]
        Mtr <- Mtr[, sel, drop = FALSE]; Mte <- Mte[, sel, drop = FALSE]
      }

      mu <- colMeans(Mtr); sdv <- matrixStats::colSds(Mtr); sdv[sdv <= 0] <- 1
      Mtr <- scale(Mtr, mu, sdv); Mte <- scale(Mte, mu, sdv)

      if (nm == "pca") {                                 # basis refit in-fold
        pc  <- prcomp(Mtr, center = FALSE, scale. = FALSE)
        np  <- min(N_PC, ncol(pc$rotation))
        Mtr <- pc$x[, seq_len(np), drop = FALSE]
        Mte <- Mte %*% pc$rotation[, seq_len(np), drop = FALSE]
      }

      if (ncol(Mtr) < 2) { row[[nm]] <- NA_real_; next }
      cv <- cv.glmnet(Mtr, ytr, family = "binomial", alpha = ALPHA, nfolds = 5)
      row[[nm]] <- as.numeric(predict(cv, Mte, s = "lambda.min", type = "response"))
    }
    preds[[length(preds) + 1]] <- row
  }
  cat(sprintf("repeat %d done (%s)\n", rep_i,
              format(round(difftime(Sys.time(), t0, units = "mins"), 1))))
}
pred_dt <- rbindlist(preds)
fwrite(pred_dt, "out/04b_representation_preds.csv")

## =====================================================================
## The grid — reported as a grid, not ranked into a winner
## =====================================================================
agg <- pred_dt[, c(list(y = y[1]),
                   lapply(.SD, function(z) mean(z, na.rm = TRUE))),
               by = sample, .SDcols = names(reps)]

B <- 2000; n <- nrow(agg)
idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = n)
grid <- rbindlist(lapply(names(reps), function(nm) {
  b <- vapply(seq_len(B), function(j) {
         i <- idx[, j]; auc_fast(agg$y[i], agg[[nm]][i]) }, numeric(1))
  data.table(representation = nm,
             role     = ifelse(nm == PRIMARY, "PRIMARY (pre-specified)", "secondary"),
             # the dimension actually handed to the model, after the
             # in-fold filter or projection, not the raw input width
             features = if (nm == "pca") min(N_PC, ncol(reps[[nm]]))
                        else if (nm == "topvar") min(TOP_K, ncol(reps[[nm]]))
                        else ncol(reps[[nm]]),
             AUROC = auc_fast(agg$y, agg[[nm]]),
             lo = quantile(b, .025, na.rm = TRUE), hi = quantile(b, .975, na.rm = TRUE),
             PR_AUC = pr_auc(agg$y, agg[[nm]]),
             Brier = mean((agg[[nm]] - agg$y)^2, na.rm = TRUE))
}))
setorder(grid, -AUROC)
fwrite(grid, "out/04b_representation_grid.csv")

cat("\n")
print(as.data.frame(grid), row.names = FALSE)
cat("\nReport this whole grid. The row marked PRIMARY is the one whose number\n")
cat("goes in the abstract; the others are secondary and are described as such.\n")
cat("If a secondary row beats it, that is an interesting observation to discuss,\n")
cat("NOT a result to promote — it was selected on the same events it is scored on.\n")
cat("\nwrote out/04b_representation_grid.csv\n")
cat("NEXT: 04_models_cv.R\n")
