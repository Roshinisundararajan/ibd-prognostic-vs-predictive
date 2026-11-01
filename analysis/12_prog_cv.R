# =====================================================================
# 12_prog_cv.R  —  fully out-of-sample check of the script-11 models
#
# Script 11 used Harrell's bootstrap: the logistic model is refitted in
# every resample, but PC1 was computed once on all 550 patients and the
# calibration / decision curves were drawn from in-sample predictions.
# PC1 never sees the outcome, so this is not outcome leakage, but the
# full PIPELINE was not re-learned. This script re-learns everything
# inside each training fold and evaluates ONLY held-out predictions:
#
#   per outer fold:  median imputation -> scaling -> top-2000 variance
#                    genes -> PCA basis -> PC1 sign (vs Mayo) -> PC1
#                    scaling -> model fit        (all on training rows)
#   held-out fold:   projected / scaled with the training parameters
#
# Two clinical baselines, each with and without the biopsy score:
#   PRIMARY   (proposal benchmark)  7 routine clinical variables + arm,
#              ridge logistic (8-9 parameters exceed the ~5 that 63
#              events support unpenalised, so ridge is used for BOTH
#              the baseline and the +PC1 model)
#   SECONDARY (sample-size-limited) arm + Mayo + prior anti-TNF +
#              calprotectin, plain logistic = script 11's M0 / M1
#   EXPLORATORY (post hoc)          secondary + OXPHOS instead of PC1
#
# IN   out/02_cohort_GSE206285.rds, out/06b_hallmark_scores.csv
# OUT  out/12_oof_predictions.csv, out/12_performance.csv,
#      out/12_added_value.csv, out/12_calibration.png,
#      out/12_decision_curve.png
# Runtime: ~3-6 min.
# =====================================================================

suppressPackageStartupMessages({
  library(data.table); library(matrixStats); library(glmnet)
})
source("_helpers.R")
set.seed(20261021)

K <- 5; R <- 10; B <- 2000; TOPV <- 2000

co  <- readRDS("out/02_cohort_GSE206285.rds")
y   <- co$y; ids <- co$ids; n <- length(y)
W   <- as.integer(co$arm == "treated")
X   <- t(co$expr)                                   # patients x genes
stopifnot(identical(rownames(X), ids), nrow(X) == n)

clin7 <- c("Mayo_score", "Disease_duration", "Neutrophils_LP",
           "Neutrophils_Epithelium", "Anti_TNF_history", "CRP_log2",
           "Calprotectin_log2")
clin4 <- c("Mayo_score", "Anti_TNF_history", "Calprotectin_log2")
stopifnot(all(clin7 %in% names(co$clin)))
C7 <- as.matrix(co$clin[, clin7]); C4 <- as.matrix(co$clin[, clin4])
cat("missing values per clinical variable (imputed with TRAINING-fold medians):\n")
print(colSums(is.na(C7)))

H  <- fread("out/06b_hallmark_scores.csv", data.table = FALSE)
OX <- H$OXIDATIVE_PHOSPHORYLATION[match(ids, H$sample)]
stopifnot(!anyNA(OX))

folds_arm_y <- function() {                      # stratify on arm x outcome
  g <- paste(W, y); f <- integer(n)
  for (s in unique(g)) { i <- which(g == s); f[i] <- sample(rep(seq_len(K), length.out = length(i))) }
  f
}
prep <- function(A, tr, te) {                    # impute + scale, train-only
  med <- impute_fit(A[tr, , drop = FALSE])
  Atr <- impute_apply(A[tr, , drop = FALSE], med); Ate <- impute_apply(A[te, , drop = FALSE], med)
  mu <- colMeans(Atr); s <- colSds(Atr); s[!is.finite(s) | s == 0] <- 1
  list(tr = scale(Atr, mu, s), te = scale(Ate, mu, s))
}
zpair <- function(a, b) { m <- mean(a); s <- sd(a); if (!is.finite(s) || s == 0) s <- 1
                          list(tr = (a - m) / s, te = (b - m) / s) }
fit_glm <- function(Xtr, ytr, Xte) {
  d <- data.frame(y_ = ytr, Xtr)
  m <- suppressWarnings(glm(y_ ~ ., data = d, family = binomial))
  as.numeric(predict(m, newdata = data.frame(Xte), type = "response"))
}
fit_ridge <- function(Xtr, ytr, Xte) {
  fid <- integer(length(ytr))                    # stratified inner folds
  for (cl in 0:1) { i <- which(ytr == cl); fid[i] <- sample(rep(1:5, length.out = length(i))) }
  m <- cv.glmnet(as.matrix(Xtr), ytr, family = "binomial", alpha = 0,
                 foldid = fid, type.measure = "deviance")
  as.numeric(predict(m, as.matrix(Xte), s = "lambda.min", type = "response"))
}

out <- list(); t0 <- Sys.time()
for (r in seq_len(R)) {
  f <- folds_arm_y()
  for (k in seq_len(K)) {
    tr <- which(f != k); te <- which(f == k)
    P7 <- prep(C7, tr, te); P4 <- prep(C4, tr, te)
    # biopsy score: every step learned on training patients only
    top <- head(order(colVars(X[tr, ]), decreasing = TRUE), TOPV)
    pca <- prcomp(X[tr, top], center = TRUE, scale. = TRUE, rank. = 1)
    ptr <- pca$x[, 1]; pte <- predict(pca, X[te, top, drop = FALSE])[, 1]
    mt  <- co$clin$Mayo_score[tr]
    if (cor(ptr, mt, method = "spearman", use = "complete.obs") < 0) { ptr <- -ptr; pte <- -pte }
    PC <- zpair(ptr, pte); OXz <- zpair(OX[tr], OX[te])

    A7tr <- cbind(arm = W[tr], P7$tr); A7te <- cbind(arm = W[te], P7$te)
    A4tr <- cbind(arm = W[tr], P4$tr); A4te <- cbind(arm = W[te], P4$te)
    out[[length(out) + 1]] <- data.table(rep = r, fold = k, sample = ids[te], y = y[te], arm = W[te],
      P_clin7       = fit_ridge(A7tr, y[tr], A7te),
      P_clin7_pc1   = fit_ridge(cbind(A7tr, pc1 = PC$tr), y[tr], cbind(A7te, pc1 = PC$te)),
      S_clin4       = fit_glm(A4tr, y[tr], A4te),
      S_clin4_pc1   = fit_glm(cbind(A4tr, pc1 = PC$tr), y[tr], cbind(A4te, pc1 = PC$te)),
      E_clin4_oxphos_posthoc = fit_glm(cbind(A4tr, ox = OXz$tr), y[tr], cbind(A4te, ox = OXz$te)))
  }
  cat(sprintf("repeat %d of %d (%.1f min)\n", r, R, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
O <- rbindlist(out); fwrite(O, "out/12_oof_predictions.csv")
models <- c("P_clin7", "P_clin7_pc1", "S_clin4", "S_clin4_pc1", "E_clin4_oxphos_posthoc")
A <- O[, lapply(.SD, mean), by = .(sample, y), .SDcols = models]   # one row per patient

metrics <- function(yy, p) {
  lp <- qlogis(pmin(pmax(p, 1e-8), 1 - 1e-8))
  c(AUROC = auc_fast(yy, p), Brier = mean((p - yy)^2),
    slope = unname(coef(glm(yy ~ lp, family = binomial))[2]),
    CITL  = unname(coef(glm(yy ~ offset(lp), family = binomial))[1]))
}
idx <- matrix(sample.int(nrow(A), nrow(A) * B, replace = TRUE), nrow = nrow(A))
perf <- rbindlist(lapply(models, function(m) {
  bs <- apply(idx, 2, function(i) auc_fast(A$y[i], A[[m]][i]))
  v  <- metrics(A$y, A[[m]])
  data.table(model = m, AUROC = v["AUROC"], AUROC_lo = quantile(bs, .025),
             AUROC_hi = quantile(bs, .975), Brier = v["Brier"],
             cal_slope = v["slope"], CITL = v["CITL"])
}))
fwrite(perf, "out/12_performance.csv")

pairs <- list(c("P_clin7", "P_clin7_pc1"), c("S_clin4", "S_clin4_pc1"),
              c("S_clin4", "E_clin4_oxphos_posthoc"))
added <- rbindlist(lapply(pairs, function(pq) {
  a <- A[[pq[1]]]; b <- A[[pq[2]]]
  dA <- apply(idx, 2, function(i) auc_fast(A$y[i], b[i]) - auc_fast(A$y[i], a[i]))
  dB <- apply(idx, 2, function(i) mean((b[i] - A$y[i])^2) - mean((a[i] - A$y[i])^2))
  data.table(baseline = pq[1], model = pq[2],
             dAUROC = auc_fast(A$y, b) - auc_fast(A$y, a),
             dAUROC_lo = quantile(dA, .025), dAUROC_hi = quantile(dA, .975),
             dBrier = mean((b - A$y)^2) - mean((a - A$y)^2),
             dBrier_lo = quantile(dB, .025), dBrier_hi = quantile(dB, .975))
}))
fwrite(added, "out/12_added_value.csv")

png("out/12_calibration.png", 1600, 800, res = 150); par(mfrow = c(1, 2))
for (pq in pairs[1:2]) {
  plot(NA, xlim = c(0, .45), ylim = c(0, .45), xlab = "predicted (held-out)",
       ylab = "observed proportion", main = paste(pq, collapse = " vs "))
  abline(0, 1, lty = 2)
  for (j in 1:2) { p <- A[[pq[j]]]; o <- order(p); lo <- loess(A$y ~ p, span = .75, degree = 1)
                   lines(p[o], predict(lo)[o], lwd = 2, col = c("steelblue", "firebrick")[j]) }
  legend("topleft", pq, col = c("steelblue", "firebrick"), lwd = 2, bty = "n")
}
dev.off()

th <- seq(0.02, 0.40, by = 0.01); prev <- mean(A$y)
nb <- sapply(models[1:4], function(m) sapply(th, function(t) net_benefit(A$y, A[[m]], t)))
png("out/12_decision_curve.png", 1100, 800, res = 150)
cols <- c("steelblue", "firebrick", "darkgreen", "darkorange", "grey40", "black")
matplot(th, cbind(nb, prev - (1 - prev) * th / (1 - th), 0), type = "l", lwd = 2,
        lty = c(1, 1, 2, 2, 2, 3), col = cols, ylim = c(-0.02, 0.11),
        xlab = "threshold probability of remission", ylab = "net benefit",
        main = "Decision curve, held-out predictions")
legend("topright", c(models[1:4], "treat all", "treat none"), col = cols,
       lty = c(1, 1, 2, 2, 2, 3), lwd = 2, bty = "n")
dev.off()

cat("\n==== held-out performance (", R, "x", K, "-fold CV, everything re-learned in-fold) ====\n")
print(as.data.frame(perf), digits = 3, row.names = FALSE)
cat("\n==== does the biopsy add? (paired patient bootstrap; dBrier < 0 = better) ====\n")
print(as.data.frame(added), digits = 3, row.names = FALSE)
cat("\nwrote out/12_*\n")
