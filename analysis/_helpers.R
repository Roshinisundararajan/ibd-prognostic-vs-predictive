# =====================================================================
# _helpers.R — small functions, sourced by scripts 02-07.
# No side effects, no file writing, no globals.
# =====================================================================

## ---- fast row/column statistics -------------------------------------
row_var <- function(X) { m <- rowMeans(X); rowSums((X - m)^2) / (ncol(X) - 1) }
row_sd  <- function(X) sqrt(row_var(X))

## ---- Mann-Whitney AUC ------------------------------------------------
# Rank-based, so it costs nothing. Used inside the bootstrap, where
# pROC::roc() would be 1000x too slow.
auc_fast <- function(y, p) {
  ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(p)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

## ---- two-sample pooled-variance t statistic, one per COLUMN ---------
# X is samples x features. This replaces apply(X, 2, ...) which is far
# too slow to call once per cross-validation fold.
#
# Note the denominator: the pooled WITHIN-GROUP sd. Using the overall sd
# of the feature (a tempting shortcut) penalises exactly the genes that
# separate the groups, which is backwards.
col_t <- function(X, y) {
  g1 <- y == 1; n1 <- sum(g1); n0 <- sum(!g1)
  if (n1 < 3 || n0 < 3) return(rep(NA_real_, ncol(X)))
  X1 <- X[g1, , drop = FALSE]; X0 <- X[!g1, , drop = FALSE]
  m1 <- colMeans(X1); m0 <- colMeans(X0)
  v1 <- matrixStats::colVars(X1); v0 <- matrixStats::colVars(X0)
  se <- sqrt(((n1 - 1) * v1 + (n0 - 1) * v0) / (n1 + n0 - 2)) *
        sqrt(1 / n1 + 1 / n0)
  se[!is.na(se) & se <= 0] <- NA_real_
  (m1 - m0) / se
}

## ---- clinical design matrix, MISSING VALUES PRESERVED ---------------
# model.matrix() silently drops rows containing NA, which de-synchronises
# the clinical block from the expression matrix and the outcome vector.
# That failure is silent and produces plausible-looking nonsense, so the
# dummy coding is done by hand here and NAs are carried through to be
# imputed inside the training fold.
build_clin_mm <- function(clin) {
  out <- list()
  for (j in names(clin)) {
    v <- clin[[j]]
    if (is.numeric(v) || is.logical(v)) {
      out[[j]] <- as.numeric(v)
    } else {
      v <- factor(as.character(v))
      if (nlevels(v) < 2) next
      for (l in levels(v)[-1])
        out[[paste0(j, "_", make.names(l))]] <- as.integer(v == l)  # NA stays NA
    }
  }
  if (!length(out)) return(matrix(numeric(0), nrow = nrow(clin), ncol = 0))
  M <- do.call(cbind, out)
  rownames(M) <- rownames(clin)
  keep <- apply(M, 2, function(z) length(unique(z[!is.na(z)])) > 1)
  M[, keep, drop = FALSE]
}

## ---- in-fold imputation ---------------------------------------------
# Medians are learned on TRAINING ROWS ONLY and then applied to both.
impute_fit <- function(M) {
  med <- apply(M, 2, median, na.rm = TRUE)
  med[!is.finite(med)] <- 0            # column entirely missing in this fold
  med
}
impute_apply <- function(M, med) {
  for (j in seq_len(ncol(M))) {
    k <- is.na(M[, j])
    if (any(k)) M[k, j] <- med[j]
  }
  M
}

## ---- binning that survives tied predictions -------------------------
# quantile() on a prediction vector with many ties returns duplicate
# break points, and cut() then errors out. This is a real failure mode
# for the clinical-only model, which has few distinct values.
safe_bins <- function(p, nb = 10) {
  br <- unique(quantile(p, seq(0, 1, length.out = nb + 1), na.rm = TRUE))
  if (length(br) < 3) return(rep(1L, length(p)))
  cut(p, breaks = br, include.lowest = TRUE, labels = FALSE)
}

## ---- expected calibration error -------------------------------------
ece <- function(y, p, nb = 10) {
  ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  b  <- safe_bins(p, nb)
  w  <- as.numeric(table(b)) / length(b)
  sum(abs(tapply(y, b, mean) - tapply(p, b, mean)) * w)
}

## ---- precision-recall AUC -------------------------------------------
# Written out so the pipeline does not need PRROC. Equivalent to the
# "integral" variant: trapezoid under precision plotted against recall.
pr_auc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  if (sum(y) == 0 || sum(y) == length(y)) return(NA_real_)
  o <- order(p, decreasing = TRUE); y <- y[o]
  tp <- cumsum(y); fp <- cumsum(1 - y)
  rec  <- tp / sum(y)
  prec <- tp / (tp + fp)
  sum(diff(c(0, rec)) * (c(prec[1], prec[-length(prec)]) + prec) / 2)
}

## ---- ROC curve coordinates ------------------------------------------
roc_points <- function(y, p) {
  ok <- !is.na(y) & !is.na(p); y <- y[ok]; p <- p[ok]
  o <- order(p, decreasing = TRUE); y <- y[o]
  data.frame(fpr = c(0, cumsum(1 - y) / sum(1 - y)),
             tpr = c(0, cumsum(y) / sum(y)))
}

## ---- net benefit (decision curve) -----------------------------------
net_benefit <- function(y, p, th) {
  tp <- sum(p >= th & y == 1); fp <- sum(p >= th & y == 0); n <- length(y)
  tp / n - (fp / n) * (th / (1 - th))
}

## ---- exclusion log (CONSORT-style) ----------------------------------
excl_new <- function() data.frame(step = character(0), n_before = integer(0),
                                  n_excluded = integer(0), n_after = integer(0),
                                  reason = character(0))
excl_add <- function(tbl, step, n_before, n_after, reason)
  rbind(tbl, data.frame(step = step, n_before = n_before,
                        n_excluded = n_before - n_after, n_after = n_after,
                        reason = reason))

## ---- "key: value" characteristics parsing ---------------------------
# GEO hides the useful fields inside characteristics_ch1.* as free text.
# Turning them into named columns here is what removes the "open the CSV
# and guess which column is the outcome" step from script 02.
parse_characteristics <- function(pd) {
  ch <- grep("^characteristics_ch1", colnames(pd), value = TRUE)
  if (!length(ch)) return(data.frame(gsm = rownames(pd)))
  rows <- do.call(rbind, lapply(ch, function(cc)
    data.frame(gsm = rownames(pd), raw = as.character(pd[[cc]]),
               stringsAsFactors = FALSE)))
  rows <- rows[!is.na(rows$raw) & grepl(":", rows$raw), ]
  rows$key <- make.names(tolower(trimws(sub(":.*$", "", rows$raw))))
  rows$val <- trimws(sub("^[^:]*:\\s*", "", rows$raw))
  keys <- unique(rows$key)
  wide <- data.frame(gsm = rownames(pd), stringsAsFactors = FALSE)
  for (k in keys) {
    sub_k <- rows[rows$key == k, ]
    wide[[k]] <- sub_k$val[match(wide$gsm, sub_k$gsm)]
  }
  # GEO writes the literal string "NA" rather than a missing value.
  for (k in keys) wide[[k]][wide[[k]] %in% c("NA", "na", "")] <- NA
  wide
}
