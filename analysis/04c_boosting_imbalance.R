# =====================================================================
# 04c_boosting_imbalance.R
#
# Closes two promises in proposal section 7 and one in section 8:
#   (a) gradient boosting (gbm) under the SAME leak-free CV as script 04
#   (b) class-imbalance handling, compared head to head:
#         none | class weights | random oversampling INSIDE the training
#         fold | random oversampling BEFORE the split (the leak demo)
#   (c) F1, sensitivity, specificity, PPV at a threshold
#
# Feature set is identical to enet_both in script 04: arm + clinical +
# top-2000 genes by in-fold t-statistic, everything fitted on training
# rows only. Only the learner and the imbalance handling change, so any
# difference is attributable to them.
#
# Thresholds for F1: two are reported.
#   prev   = training-fold prevalence (fixed in advance, honest)
#   youden = Youden-optimal on the pooled out-of-fold predictions
#            (post hoc and optimistic; label it as such)
#
# IN   out/02_cohort_GSE206285.rds
# OUT  out/04c_predictions.csv, out/04c_metrics.csv,
#      out/04c_comparison.csv, out/04c_leak_demo.csv
# Runtime: about 10-20 min on 2 cores (gbm is the slow part).
# =====================================================================

suppressPackageStartupMessages({
  library(glmnet); library(gbm); library(data.table); library(matrixStats)
})
source("_helpers.R")
set.seed(20261019)

GSE_ID    <- "GSE206285"
N_FOLDS   <- 5
N_REPEATS <- 3
TOP_K     <- 2000
ALPHA     <- 0.5
B         <- 2000

co  <- readRDS(file.path("out", paste0("02_cohort_", GSE_ID, ".rds")))
ex  <- t(co$expr); y <- co$y; ids <- co$ids
arm_mm  <- matrix(as.integer(co$arm == "treated"), ncol = 1,
                  dimnames = list(ids, "arm_treated"))
clin_mm <- build_clin_mm(co$clin)
stopifnot(nrow(clin_mm) == nrow(ex), identical(rownames(ex), ids))
cat(sprintf("samples %d | events %d (%.1f%%)\n", length(y), sum(y), 100 * mean(y)))

strat_folds <- function(y, k) {
  f <- integer(length(y))
  for (cl in unique(y)) { i <- which(y == cl); f[i] <- sample(rep(seq_len(k), length.out = length(i))) }
  f
}

# build the in-fold design: arm + imputed/scaled clinical + filtered/scaled genes
make_design <- function(tr, te) {
  med <- impute_fit(clin_mm[tr, , drop = FALSE])
  Ctr <- impute_apply(clin_mm[tr, , drop = FALSE], med)
  Cte <- impute_apply(clin_mm[te, , drop = FALSE], med)
  cm <- colMeans(Ctr); cs <- colSds(Ctr); cs[cs <= 0] <- 1
  Ctr <- scale(Ctr, cm, cs); Cte <- scale(Cte, cm, cs)
  ts  <- abs(col_t(ex[tr, , drop = FALSE], y[tr]))
  sel <- names(sort(ts, decreasing = TRUE))[seq_len(TOP_K)]
  mu  <- colMeans(ex[tr, sel]); sd <- colSds(ex[tr, sel]); sd[sd <= 0] <- 1
  list(tr = cbind(arm_mm[tr, , drop = FALSE], Ctr, scale(ex[tr, sel], mu, sd)),
       te = cbind(arm_mm[te, , drop = FALSE], Cte, scale(ex[te, sel], mu, sd)))
}

class_w <- function(yy) ifelse(yy == 1, sum(yy == 0) / sum(yy == 1), 1)

oversample <- function(idx, yy) {          # duplicate minority rows to balance
  pos <- idx[yy[idx] == 1]; neg <- idx[yy[idx] == 0]
  c(neg, pos, sample(pos, length(neg) - length(pos), replace = TRUE))
}

fit_enet <- function(X, yy, w = NULL) {
  cv <- cv.glmnet(X, yy, family = "binomial", alpha = ALPHA, nfolds = 5, weights = w)
  function(Z) as.numeric(predict(cv, Z, s = "lambda.min", type = "response"))
}

fit_gbm <- function(X, yy, w = NULL) {
  d <- data.frame(y_ = yy, X, check.names = FALSE)
  m <- gbm(y_ ~ ., data = d, distribution = "bernoulli", weights = w,
           n.trees = 1000, interaction.depth = 2, shrinkage = 0.01,
           bag.fraction = 0.5, n.minobsinnode = 10, cv.folds = 3,
           n.cores = 1, verbose = FALSE)
  nt <- suppressWarnings(gbm.perf(m, method = "cv", plot.it = FALSE))
  function(Z) as.numeric(predict(m, data.frame(Z, check.names = FALSE),
                                 n.trees = nt, type = "response"))
}

## ---- honest CV -------------------------------------------------------
preds <- list(); thr_prev <- list(); t0 <- Sys.time()
for (r in seq_len(N_REPEATS)) {
  folds <- strat_folds(y, N_FOLDS)
  for (k in seq_len(N_FOLDS)) {
    tr <- which(folds != k); te <- which(folds == k); ytr <- y[tr]
    D  <- make_design(tr, te)
    w  <- class_w(ytr)
    os <- match(oversample(tr, y), tr)               # positions within tr
    p <- list(
      enet_none     = fit_enet(D$tr, ytr)(D$te),
      enet_weights  = fit_enet(D$tr, ytr, w)(D$te),
      enet_oversamp = fit_enet(D$tr[os, ], ytr[os])(D$te),
      gbm_none      = fit_gbm(D$tr, ytr)(D$te),
      gbm_weights   = fit_gbm(D$tr, ytr, w)(D$te))
    preds[[length(preds) + 1]] <- data.table(rep = r, fold = k, sample = ids[te],
                                             y = y[te], arm = co$arm[te], as.data.table(p))
    thr_prev[[length(thr_prev) + 1]] <- data.table(rep = r, fold = k, sample = ids[te],
                                                   thr = mean(ytr))
    cat(sprintf("repeat %d fold %d (%.1f min)\n", r, k,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
P <- rbindlist(preds); fwrite(P, "out/04c_predictions.csv")
models <- setdiff(names(P), c("rep", "fold", "sample", "y", "arm"))

## ---- patient-level aggregation and metrics ---------------------------
# Weighted/oversampled models are not calibrated to the true prevalence by
# design, so the "prev" threshold is applied to rank-percentile instead:
# flag the top (prevalence x 100)% of predicted risk. Same rule for all models.
A <- P[, lapply(.SD, mean), by = .(sample, y), .SDcols = models]
prev <- mean(A$y)

thr_metrics <- function(yy, p, cut) {
  pred <- as.integer(p >= cut)
  tp <- sum(pred == 1 & yy == 1); fp <- sum(pred == 1 & yy == 0)
  fn <- sum(pred == 0 & yy == 1); tn <- sum(pred == 0 & yy == 0)
  sens <- tp / (tp + fn); spec <- tn / (tn + fp); ppv <- if (tp + fp) tp / (tp + fp) else NA
  c(sens = sens, spec = spec, ppv = ppv,
    f1 = if (is.na(ppv) || ppv + sens == 0) 0 else 2 * ppv * sens / (ppv + sens))
}
youden_cut <- function(yy, p) {
  cs <- sort(unique(p)); j <- sapply(cs, function(c) { m <- thr_metrics(yy, p, c); m[1] + m[2] - 1 })
  cs[which.max(j)]
}

met <- rbindlist(lapply(models, function(m) {
  p <- A[[m]]
  bs <- replicate(B, { i <- sample.int(nrow(A), replace = TRUE); auc_fast(A$y[i], p[i]) })
  top <- thr_metrics(A$y, p, quantile(p, 1 - prev, type = 1))
  yo  <- thr_metrics(A$y, p, youden_cut(A$y, p))
  data.table(model = m, AUROC = auc_fast(A$y, p),
             AUROC_lo = quantile(bs, .025), AUROC_hi = quantile(bs, .975),
             PR_AUC = pr_auc(A$y, p), Brier = mean((p - A$y)^2),
             mean_pred = mean(p), prevalence = prev,
             F1_topprev = top["f1"], sens_topprev = top["sens"], spec_topprev = top["spec"], PPV_topprev = top["ppv"],
             F1_youden_posthoc = yo["f1"], sens_youden = yo["sens"], spec_youden = yo["spec"], PPV_youden = yo["ppv"])
}))
fwrite(met, "out/04c_metrics.csv")

## ---- paired bootstrap vs the reference (enet_none = enet_both in 04) ---
ref <- "enet_none"
cmp <- rbindlist(lapply(setdiff(models, ref), function(m) {
  d <- replicate(B, { i <- sample.int(nrow(A), replace = TRUE)
    auc_fast(A$y[i], A[[m]][i]) - auc_fast(A$y[i], A[[ref]][i]) })
  data.table(model = m, vs = ref, dAUROC = auc_fast(A$y, A[[m]]) - auc_fast(A$y, A[[ref]]),
             lo = quantile(d, .025), hi = quantile(d, .975),
             p_boot = 2 * min(mean(d <= 0), mean(d >= 0)))
}))
fwrite(cmp, "out/04c_comparison.csv")

## ---- the leak: oversample the WHOLE cohort, then split ----------------
# Copies of the same minority patient land in both training and test folds,
# so the model is scored partly on patients it has memorised.
leak <- rbindlist(lapply(seq_len(N_REPEATS), function(r) {
  big <- oversample(seq_along(y), y)               # row indices, duplicates included
  yb  <- y[big]; fb <- strat_folds(yb, N_FOLDS); out <- numeric(length(big))
  for (k in seq_len(N_FOLDS)) {
    trb <- which(fb != k); teb <- which(fb == k)
    D <- make_design(big[trb], big[teb])
    out[teb] <- fit_enet(D$tr, yb[trb])(D$te)
  }
  data.table(rep = r, AUROC_leaky = auc_fast(yb, out))
}))
leak[, AUROC_honest_enet_oversamp := met[model == "enet_oversamp", AUROC]]
fwrite(leak, "out/04c_leak_demo.csv")

cat("\n==== metrics (patient level, mean over repeats) ====\n")
print(as.data.frame(met[, .(model, AUROC, AUROC_lo, AUROC_hi, PR_AUC, Brier, mean_pred,
                            F1_topprev, F1_youden_posthoc)]), digits = 3, row.names = FALSE)
cat("\n==== paired bootstrap vs enet_none ====\n")
print(as.data.frame(cmp), digits = 3, row.names = FALSE)
cat("\n==== pre-split oversampling leak ====\n")
print(as.data.frame(leak), digits = 3, row.names = FALSE)
cat(sprintf("\nprevalence %.3f, so a no-skill model has PR-AUC ~%.3f and F1 ~%.3f at the top-prevalence cut\n",
            prev, prev, prev))
