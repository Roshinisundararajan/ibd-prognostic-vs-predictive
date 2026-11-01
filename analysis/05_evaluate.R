# =====================================================================
# 05_evaluate.R
#
# Discrimination, calibration and clinical usefulness, with intervals
# that mean what they say.
#
# Two things here are easy to get wrong and both were wrong in the first
# draft of this pipeline:
#
#  1. With N_REPEATS > 1 every patient appears several times in
#     04_cv_predictions.csv. Bootstrapping those ROWS treats one patient
#     as several independent people and returns intervals that are far too
#     narrow. Predictions are averaged to one per patient first, and the
#     bootstrap resamples PATIENTS.
#  2. Comparing two models with DeLong's test on those duplicated rows has
#     the same problem. The comparison here is a paired bootstrap on the
#     same resampled patients, which is what the duplicated-row test was
#     pretending to be.
#
# IN   out/04_cv_predictions.csv
# OUT  out/05_metrics.csv, out/05_model_comparison.csv
#      out/05_roc.png, out/05_calibration.png, out/05_decision_curve.png
# =====================================================================

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
source("_helpers.R")
set.seed(20261017)

B <- 2000
p <- fread("out/04_cv_predictions.csv")

LADDER <- c("arm_only", "clinical", "arm_sig", "enet_expr", "enet_both",
            "enet_both_cal", "rf", "enet_leaky")
models <- intersect(LADDER, names(p))
models <- models[vapply(models, function(m) !all(is.na(p[[m]])), logical(1))]
cat("models evaluated:", paste(models, collapse = ", "), "\n")
cat("repeats in file  :", length(unique(p$rep)), "\n")

## ---- one row per patient ---------------------------------------------
agg <- p[, c(list(y = y[1], arm = arm[1]),
             lapply(.SD, function(z) mean(z, na.rm = TRUE))),
         by = sample, .SDcols = models]
cat(sprintf("rows %d -> patients %d\n", nrow(p), nrow(agg)))
prev <- mean(agg$y)
cat(sprintf("prevalence %.3f (%d events)\n\n", prev, sum(agg$y)))

## ---- one set of bootstrap patient draws, shared by every model -------
# Sharing the draws is what makes the model-vs-model comparison paired.
n   <- nrow(agg)
idx <- matrix(sample.int(n, n * B, replace = TRUE), nrow = n)
bA  <- vapply(seq_len(B), function(b) {
  i <- idx[, b]
  vapply(models, function(m) auc_fast(agg$y[i], agg[[m]][i]), numeric(1))
}, numeric(length(models)))
rownames(bA) <- models

## ---- metric table -----------------------------------------------------
res <- rbindlist(lapply(models, function(m) {
  s <- agg[[m]]; y <- agg$y
  ci <- quantile(bA[m, ], c(0.025, 0.975), na.rm = TRUE)
  
  lp  <- qlogis(pmin(pmax(s, 1e-6), 1 - 1e-6))
  cal <- suppressWarnings(glm(y ~ lp, family = binomial))
  
  # Youden threshold. Chosen on these same predictions, so everything
  # downstream of it is optimistic — label it as post hoc in the report.
  rp  <- roc_points(y, s)
  th  <- sort(unique(s), decreasing = TRUE)[which.max(rp$tpr[-1] - rp$fpr[-1])]
  b   <- s >= th
  tp <- sum(b & y == 1); fp <- sum(b & y == 0)
  fn <- sum(!b & y == 1); tn <- sum(!b & y == 0)
  
  data.table(model = m,
             AUROC = auc_fast(y, s), AUROC_lo = ci[[1]], AUROC_hi = ci[[2]],
             PR_AUC = pr_auc(y, s), prevalence = prev,
             Brier = mean((s - y)^2),
             cal_intercept = coef(cal)[[1]], cal_slope = coef(cal)[[2]],
             OE_ratio = mean(y) / mean(s), ECE = ece(y, s),
             threshold_posthoc = th,
             sens = tp / (tp + fn), spec = tn / (tn + fp),
             PPV = tp / (tp + fp), NPV = tn / (tn + fn))
}))
fwrite(res, "out/05_metrics.csv")
print(res[, .(model, AUROC, AUROC_lo, AUROC_hi, PR_AUC, Brier, cal_slope, ECE)])

## ---- paired comparison against the clinical baseline -----------------
if ("clinical" %in% models) {
  cmp <- rbindlist(lapply(setdiff(models, "clinical"), function(m) {
    d <- bA[m, ] - bA["clinical", ]
    data.table(model = m,
               dAUROC = auc_fast(agg$y, agg[[m]]) - auc_fast(agg$y, agg$clinical),
               lo = quantile(d, 0.025, na.rm = TRUE),
               hi = quantile(d, 0.975, na.rm = TRUE),
               p_boot = 2 * min(mean(d <= 0, na.rm = TRUE),
                                mean(d >= 0, na.rm = TRUE)))
  }))
  fwrite(cmp, "out/05_model_comparison.csv")
  cat("\nchange in AUROC relative to clinical variables alone",
      "\n(paired patient-level bootstrap):\n")
  print(cmp, row.names = FALSE)
}

## ---- what the one deliberate shortcut was worth ----------------------
if (all(c("enet_leaky", "enet_expr") %in% models)) {
  d <- bA["enet_leaky", ] - bA["enet_expr", ]
  cat(sprintf(paste0("\noptimism from selecting genes outside the folds: ",
                     "%+.3f AUROC (95%% %+.3f to %+.3f)\n"),
              mean(d, na.rm = TRUE), quantile(d, .025, na.rm = TRUE),
              quantile(d, .975, na.rm = TRUE)))
  cat("That gap is the number to quote when explaining why in-fold refitting\n")
  cat("matters: it measures a mistake, not a result. Expect it to be positive\n")
  cat("where there is real signal to overfit, and near zero where there is not.\n")
}

## ---- figures ----------------------------------------------------------
honest <- setdiff(models, "enet_leaky")

roc_df <- rbindlist(lapply(honest, function(m) {
  d <- roc_points(agg$y, agg[[m]])
  d$model <- sprintf("%s (%.3f)", m, auc_fast(agg$y, agg[[m]])); d
}))
g1 <- ggplot(roc_df, aes(fpr, tpr, colour = model)) +
  geom_abline(linetype = 2, colour = "grey60") + geom_line(linewidth = 0.7) +
  coord_equal() +
  labs(title = "Out-of-fold ROC, discovery cohort",
       subtitle = "predictions averaged per patient across cross-validation repeats",
       x = "1 - specificity", y = "sensitivity", colour = "model (AUROC)") +
  theme_minimal(base_size = 12)
ggsave("out/05_roc.png", g1, width = 7.5, height = 5.5, dpi = 200)

cal_df <- rbindlist(lapply(honest, function(m) {
  b <- safe_bins(agg[[m]], 10)
  data.table(model = m,
             predicted = as.numeric(tapply(agg[[m]], b, mean)),
             observed  = as.numeric(tapply(agg$y,    b, mean)),
             n         = as.numeric(table(b)))
}))
g2 <- ggplot(cal_df, aes(predicted, observed, colour = model)) +
  geom_abline(linetype = 2, colour = "grey60") +
  geom_line(linewidth = 0.6) + geom_point(aes(size = n), alpha = 0.8) +
  geom_hline(yintercept = prev, linetype = 3, colour = "grey40") +
  labs(title = "Calibration, deciles of predicted risk",
       subtitle = sprintf("dotted line = cohort remission rate (%.1f%%)", 100 * prev),
       x = "predicted probability", y = "observed fraction") +
  theme_minimal(base_size = 12)
ggsave("out/05_calibration.png", g2, width = 7.5, height = 5.5, dpi = 200)

ths <- seq(0.01, 0.6, by = 0.01)
dc  <- rbindlist(c(
  lapply(honest, function(m) data.table(model = m, threshold = ths,
                                        nb = vapply(ths, function(t) net_benefit(agg$y, agg[[m]], t), numeric(1)))),
  list(data.table(model = "treat all",  threshold = ths,
                  nb = prev - (1 - prev) * ths / (1 - ths)),
       data.table(model = "treat none", threshold = ths, nb = 0))))
g3 <- ggplot(dc, aes(threshold, nb, colour = model)) + geom_line(linewidth = 0.7) +
  coord_cartesian(ylim = c(-0.05, max(dc$nb, na.rm = TRUE) * 1.15)) +
  labs(title = "Decision curve", x = "threshold probability", y = "net benefit",
       subtitle = "a model is useful only where it sits above both reference lines") +
  theme_minimal(base_size = 12)
ggsave("out/05_decision_curve.png", g3, width = 7.5, height = 5.5, dpi = 200)

cat("\nwrote out/05_metrics.csv, out/05_model_comparison.csv and three figures\n")
cat("NEXT: 06_signatures_interaction.R\n")
