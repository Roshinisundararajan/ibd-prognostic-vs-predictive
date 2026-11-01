# =====================================================================
# 07_replicate.R   OPTIONAL — cut this without regret if Day 2 is tight.
#
# Lock the model on the discovery cohort, apply it unchanged to the
# replication cohort. No matrix merging, no ComBat: expression is
# converted to WITHIN-SAMPLE RANKS, which are platform-independent, so a
# model fitted on one array platform can be applied to another without
# any correction step at all. That is the whole trick.
#
# Prerequisites, in order:
#   1. edit GSE_ID in 01_download_annotate.R to the replication series, run it
#   2. confirm the SPEC entry for that series in 02_build_cohort.R against
#      out/01_pheno_<GSE>.csv, edit GSE_ID there, run it
#   3. run this
#
# IN   out/02_cohort_<DISC>.rds, out/02_cohort_<REPL>.rds
# OUT  out/07_replication.csv, out/07_locked_model.rds
# =====================================================================

suppressPackageStartupMessages({
  library(glmnet); library(data.table)
})
source("_helpers.R")
set.seed(20261017)

DISC <- "GSE206285"
REPL <- "GSE207022"

f_d <- file.path("out", paste0("02_cohort_", DISC, ".rds"))
f_r <- file.path("out", paste0("02_cohort_", REPL, ".rds"))
if (!file.exists(f_r))
  stop(f_r, " does not exist. Build the replication cohort first — see the\n",
       "  three prerequisites at the top of this script. There is no shortcut.")

disc <- readRDS(f_d); repl <- readRDS(f_r)

common <- intersect(rownames(disc$expr), rownames(repl$expr))
cat(sprintf("%d genes shared across the two platforms\n", length(common)))
if (length(common) < 2000)
  cat("*** fewer than 2000 shared genes. Check that both cohorts used the\n",
      "    same symbol vocabulary before reading anything into the result.\n", sep = "")

## ---- rank transform, within each sample, within each cohort ----------
rank_within <- function(m) apply(m, 2, function(z) rank(z, ties.method = "average") / length(z))
Xd <- t(rank_within(disc$expr[common, , drop = FALSE]))
Xr <- t(rank_within(repl$expr[common, , drop = FALSE]))

## ---- lock the model on the full discovery cohort ---------------------
# Fitted once, on everything, and then frozen. Nothing from the
# replication cohort is allowed to influence it, including the choice of
# lambda.
fit <- cv.glmnet(Xd, disc$y, family = "binomial", alpha = 0.5, nfolds = 5)
saveRDS(list(fit = fit, genes = common, built = Sys.time()), "out/07_locked_model.rds")

cf   <- coef(fit, s = "lambda.min")
kept <- setdiff(rownames(cf)[which(as.numeric(cf) != 0)], "(Intercept)")
cat(sprintf("locked model uses %d genes\n", length(kept)))

## ---- apply it ---------------------------------------------------------
pr <- as.numeric(predict(fit, Xr, s = "lambda.min", type = "response"))
a  <- auc_fast(repl$y, pr)

B <- 2000; n <- length(pr)
bs <- replicate(B, { i <- sample.int(n, n, replace = TRUE)
                     auc_fast(repl$y[i], pr[i]) })
ci <- quantile(bs, c(0.025, 0.975), na.rm = TRUE)

cat(sprintf("\nreplication AUROC = %.3f (95%% %.3f to %.3f), %d patients, %d events\n",
            a, ci[[1]], ci[[2]], length(repl$y), sum(repl$y)))
cat(sprintf("replication Brier = %.4f   prevalence = %.3f\n",
            mean((pr - repl$y)^2), mean(repl$y)))
cat(sprintf("calibration slope = %.3f\n",
            coef(suppressWarnings(glm(repl$y ~ qlogis(pmin(pmax(pr, 1e-6), 1 - 1e-6)),
                                      family = binomial)))[[2]]))
cat("\nA slope well below 1 means the model transfers in RANKING but not in\n")
cat("absolute risk — common, and worth saying plainly rather than hiding.\n")

fwrite(data.table(sample = repl$ids, y = repl$y,
                  arm = if (is.null(repl$arm)) NA else repl$arm,
                  predicted = pr), "out/07_replication.csv")
cat("\nwrote out/07_replication.csv\n")
