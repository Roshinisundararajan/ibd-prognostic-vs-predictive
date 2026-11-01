# =====================================================================
# 11_prognostic_model.R
#
# A small, pre-specified PROGNOSTIC model for week-8 clinical remission,
# built the way the prognosis-research guidance asks
# (prognosisresearch.com/guidance-prognostic-models; Riley BMJ 2020;
#  Steyerberg "Clinical Prediction Models"; reported with TRIPOD+AI).
#
# Question it answers:
#   "Given what we know at baseline, how likely is THIS patient to reach
#    remission - and does the biopsy transcriptome add anything to the
#    routine clinical picture?"
#
# Design decisions, fixed BEFORE running (write them in Methods):
#   - Sample size (Riley criterion 1, shrinkage >= 0.9): 63 events support
#     only ~4-6 predictor parameters. So no 2000-gene model here.
#   - M0 (clinical): arm + Mayo score + prior anti-TNF + faecal calprotectin
#   - M1 (clinical + biopsy): M0 + ONE transcriptomic score = PC1 of the
#     2000 most variable genes. PC1 is unsupervised (never sees the
#     outcome) and script 03 showed it is the mucosal inflammation axis.
#   - M2 (sensitivity only, post hoc): M0 + OXPHOS Hallmark score, the
#     pathway that was most prognostic in 06b. Chosen AFTER seeing data,
#     so it is optimistic. Report it as exploratory.
#   - Plain logistic regression (no class-imbalance correction: 04c and
#     the guidance both show it damages calibration).
#   - Internal validation: Harrell bootstrap optimism correction, B = 500.
#   - Uniform shrinkage by the optimism-corrected calibration slope,
#     intercept re-estimated -> the final, reportable equation.
#   - Missing values (CRP/calprotectin): median imputation, outcome-blind.
#
# IN   out/02_cohort_GSE206285.rds  (+ out/06b_hallmark_scores.csv for M2)
# OUT  out/11_performance.csv     apparent vs optimism-corrected metrics
#      out/11_added_value.csv     M1 vs M0 (and M2 vs M0)
#      out/11_final_model.csv     shrunk coefficients = the calculator
#      out/11_predictions.csv     per-patient risks
#      out/11_calibration.png, out/11_decision_curve.png,
#      out/11_instability.png
# Runtime: ~2-4 minutes.
# =====================================================================

suppressPackageStartupMessages({ library(data.table); library(matrixStats) })
source("_helpers.R")
set.seed(20261020)
B <- 500

co <- readRDS("out/02_cohort_GSE206285.rds")
y  <- co$y; n <- length(y)

## ---- 1. sample size: how many predictors can 63 events support? -------
p_events <- mean(y)
lnL0  <- sum(y) * log(p_events) + (n - sum(y)) * log(1 - p_events)
maxR2 <- 1 - exp(2 * lnL0 / n)
cat(sprintf("n = %d, events = %d (%.1f%%), max Cox-Snell R2 = %.3f\n",
            n, sum(y), 100 * p_events, maxR2))
for (R2 in c(0.05, 0.075, 0.10))
  cat(sprintf("  if R2cs = %.3f -> at most %.1f predictor parameters (Riley, S = 0.9)\n",
              R2, n * (0.9 - 1) * log(1 - R2 / 0.9)))

## ---- 2. predictors ----------------------------------------------------
med <- function(v) { v[is.na(v)] <- median(v, na.rm = TRUE); v }
ex  <- co$expr
top <- head(order(rowVars(ex), decreasing = TRUE), 2000)
pc1 <- prcomp(t(ex[top, ]), scale. = TRUE)$x[, 1]
# orient PC1 so that higher = more inflamed (higher Mayo)
if (cor(pc1, co$clin$Mayo_score, method = "spearman") < 0) pc1 <- -pc1

D <- data.frame(
  y        = y,
  arm      = as.integer(co$arm == "treated"),
  mayo     = co$clin$Mayo_score,
  antiTNF  = as.integer(co$clin$Anti_TNF_history),
  calpro   = med(co$clin$Calprotectin_log2),
  pc1      = as.numeric(scale(pc1)))
H <- tryCatch({
  h <- fread("out/06b_hallmark_scores.csv", data.table = FALSE)
  as.numeric(scale(h$OXIDATIVE_PHOSPHORYLATION[match(co$ids, h$sample)]))
}, error = function(e) NULL)
if (!is.null(H)) D$oxphos <- H
cat(sprintf("missing imputed: calprotectin %d\n", sum(is.na(co$clin$Calprotectin_log2))))

forms <- list(
  M0_clinical        = y ~ arm + mayo + antiTNF + calpro,
  M1_clinical_biopsy = y ~ arm + mayo + antiTNF + calpro + pc1)
if (!is.null(H)) forms$M2_oxphos_posthoc <- y ~ arm + mayo + antiTNF + calpro + oxphos

## ---- 3. metrics -------------------------------------------------------
metrics <- function(yy, p) {
  lp <- qlogis(pmin(pmax(p, 1e-8), 1 - 1e-8))
  c(AUROC = auc_fast(yy, p),
    Brier = mean((p - yy)^2),
    slope = unname(coef(glm(yy ~ lp, family = binomial))[2]),
    CITL  = unname(coef(glm(yy ~ offset(lp), family = binomial))[1]))
}

## ---- 4. Harrell bootstrap optimism correction ---------------------------
validate <- function(f) {
  fit  <- glm(f, data = D, family = binomial)
  app  <- metrics(y, fitted(fit))
  boot <- t(replicate(B, {
    i  <- sample.int(n, replace = TRUE)
    fb <- suppressWarnings(glm(f, data = D[i, ], family = binomial))
    metrics(y[i], fitted(fb)) -
      metrics(y, predict(fb, newdata = D, type = "response"))
  }))
  opt <- colMeans(boot, na.rm = TRUE)
  list(fit = fit, apparent = app, optimism = opt, corrected = app - opt)
}

V <- lapply(forms, validate)
perf <- rbindlist(lapply(names(V), function(m) {
  v <- V[[m]]
  data.table(model = m, metric = names(v$apparent),
             apparent = v$apparent, optimism = v$optimism, corrected = v$corrected)
}))
# corrected = what to expect in new patients from the same population
# (ideal: calibration slope 1, calibration-in-the-large CITL 0)
fwrite(perf, "out/11_performance.csv")

## ---- 5. added value of the biopsy -------------------------------------
lr <- function(a, b) anova(V[[a]]$fit, V[[b]]$fit, test = "LRT")[2, "Pr(>Chi)"]
dboot <- function(a, b) {               # bootstrap CI of apparent dAUROC, paired
  quantile(replicate(1000, { i <- sample.int(n, replace = TRUE)
  auc_fast(y[i], fitted(V[[b]]$fit)[i]) - auc_fast(y[i], fitted(V[[a]]$fit)[i]) }),
  c(.025, .975))
}
added <- rbindlist(lapply(setdiff(names(V), "M0_clinical"), function(m) {
  ci <- dboot("M0_clinical", m)
  data.table(model = m, vs = "M0_clinical",
             dAUROC_apparent  = V[[m]]$apparent["AUROC"] - V$M0_clinical$apparent["AUROC"],
             dAUROC_lo = ci[1], dAUROC_hi = ci[2],
             dAUROC_corrected = V[[m]]$corrected["AUROC"] - V$M0_clinical$corrected["AUROC"],
             LR_test_p = lr("M0_clinical", m))
}))
fwrite(added, "out/11_added_value.csv")

## ---- 6. final model: uniform shrinkage, intercept re-estimated ----------
shrink <- function(v) {
  s  <- min(1, v$corrected["slope"])
  b  <- coef(v$fit); b[-1] <- b[-1] * s
  X  <- model.matrix(v$fit)
  lp <- as.numeric(X[, -1, drop = FALSE] %*% b[-1])
  b[1] <- coef(glm(y ~ offset(lp), family = binomial))[1]
  list(coef = b, shrinkage = s, p = plogis(b[1] + lp))
}
S <- lapply(V, shrink)
final <- rbindlist(lapply(names(S), function(m)
  data.table(model = m, term = names(S[[m]]$coef),
             coef_original = coef(V[[m]]$fit), coef_shrunk = S[[m]]$coef,
             OR_shrunk = exp(S[[m]]$coef), shrinkage = S[[m]]$shrinkage)))
fwrite(final, "out/11_final_model.csv")
fwrite(data.table(sample = co$ids, y = y, arm = co$arm,
                  as.data.table(lapply(S, `[[`, "p"))), "out/11_predictions.csv")

## ---- 7. calibration plot (flexible, loess) --------------------------------
png("out/11_calibration.png", 1500, 750, res = 150); par(mfrow = c(1, 2))
for (m in names(S)[1:2]) {
  p <- S[[m]]$p; o <- order(p)
  lo <- loess(y ~ p, span = 0.75, degree = 1)
  plot(p[o], predict(lo)[o], type = "l", lwd = 2, xlim = c(0, 0.45), ylim = c(0, 0.45),
       xlab = "predicted probability of remission", ylab = "observed proportion",
       main = m); abline(0, 1, lty = 2)
  rug(p[y == 1], side = 3, col = "darkgreen"); rug(p[y == 0], side = 1, col = "grey50")
}
dev.off()

## ---- 8. decision curve -------------------------------------------------
th <- seq(0.02, 0.40, by = 0.01)
nb <- sapply(names(S), function(m) sapply(th, function(t) net_benefit(y, S[[m]]$p, t)))
nb_all <- p_events - (1 - p_events) * th / (1 - th)
png("out/11_decision_curve.png", 1100, 800, res = 150)
matplot(th, cbind(nb, nb_all, 0), type = "l", lty = c(1, 1, 1, 2, 3), lwd = 2,
        col = c("steelblue", "firebrick", "darkorange", "grey40", "black")[seq_len(ncol(nb) + 2)],
        ylim = c(-0.02, max(nb, nb_all, na.rm = TRUE) + 0.01),
        xlab = "threshold probability of remission", ylab = "net benefit",
        main = "Decision curve")
legend("topright", c(names(S), "treat all", "treat none"), bty = "n", lwd = 2,
       col = c("steelblue", "firebrick", "darkorange", "grey40", "black")[seq_len(ncol(nb) + 2)],
       lty = c(rep(1, ncol(nb)), 2, 3))
dev.off()

## ---- 9. instability plot (Riley 2023): how much would a patient's risk ----
##      change if the model were rebuilt on a different sample of 550?
inst <- function(f, Bi = 200) {
  p0 <- fitted(glm(f, data = D, family = binomial))
  pb <- replicate(Bi, { i <- sample.int(n, replace = TRUE)
  predict(suppressWarnings(glm(f, data = D[i, ], family = binomial)), D, type = "response") })
  list(p0 = p0, pb = pb, MAPE = mean(abs(pb - p0)))
}
I <- lapply(forms[1:2], inst)
png("out/11_instability.png", 1500, 750, res = 150); par(mfrow = c(1, 2))
for (m in names(I)) {
  matplot(I[[m]]$p0, I[[m]]$pb[, 1:100], pch = ".", col = rgb(0, 0, 0, .15),
          xlab = "risk from the final model", ylab = "risk from bootstrap models",
          main = sprintf("%s  (MAPE %.3f)", m, I[[m]]$MAPE), xlim = c(0, .5), ylim = c(0, .5))
  abline(0, 1, col = "red")
}
dev.off()

## ---- 10. print ---------------------------------------------------------
cat("\n==== internal validation (bootstrap optimism-corrected, B = 500) ====\n")
print(dcast(perf, model ~ metric, value.var = "corrected"), digits = 3)
cat("(apparent values and optimism are in out/11_performance.csv)\n")
cat("\n==== does the biopsy add to the clinical model? ====\n")
print(as.data.frame(added), digits = 3, row.names = FALSE)
cat("\n==== final shrunk model (this is the calculator) ====\n")
print(as.data.frame(final[, .(model, term, OR_shrunk, shrinkage)]), digits = 3, row.names = FALSE)
cat(sprintf("\nprediction instability (mean absolute change in risk): M0 %.3f, M1 %.3f\n",
            I[[1]]$MAPE, I[[2]]$MAPE))
cat("\nwrote out/11_* (performance, added_value, final_model, predictions, 3 figures)\n")
