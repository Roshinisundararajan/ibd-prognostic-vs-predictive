# =====================================================================
# 04_models_cv.R
#
# Repeated stratified cross-validation over a model ladder.
#
# Every learned step — clinical imputation, the univariate gene filter,
# scaling, elastic-net selection, the recalibration map — is refit inside
# the training fold. The test fold is touched only for prediction.
#
# One deliberate exception: `enet_leaky` selects its genes using the whole
# dataset, test folds included. It is there to be WRONG, so the report can
# quantify how much optimism that one shortcut buys. Do not read it as a
# result.
#
# IN   out/02_cohort_<GSE>.rds
# OUT  out/04_cv_predictions.csv   out-of-fold predictions, one row per
#                                  patient per repeat per model
#      out/04_stability.csv        selection frequency per gene
#
# Runtime: roughly 10-30 min at 550 samples. Set N_REPEATS <- 1 if short.
# =====================================================================

suppressPackageStartupMessages({
  library(glmnet); library(ranger); library(data.table); library(matrixStats)
})
source("_helpers.R")
set.seed(20261017)

GSE_ID    <- "GSE206285"
N_FOLDS   <- 5
N_REPEATS <- 3        # drop to 1 if short on time; say which you used
TOP_K     <- 2000     # univariate filter width, refit per fold
ALPHA     <- 0.5      # elastic net mixing
NTHREAD   <- max(1, parallel::detectCores() - 1)

co   <- readRDS(file.path("out", paste0("02_cohort_", GSE_ID, ".rds")))
ex   <- t(co$expr)                       # samples x genes
y    <- co$y
arm  <- co$arm
ids  <- co$ids
stopifnot(nrow(ex) == length(y), identical(rownames(ex), ids))

## ---- clinical design matrix ------------------------------------------
# build_clin_mm keeps NAs. model.matrix() would have dropped those rows and
# silently shifted the clinical block out of step with y.
clin_mm <- build_clin_mm(co$clin)
stopifnot(nrow(clin_mm) == nrow(ex))     # the assertion that catches it
has_clin <- ncol(clin_mm) > 0
cat(sprintf("samples %d | events %d (%.1f%%) | genes %d | clinical columns %d\n",
            nrow(ex), sum(y), 100 * mean(y), ncol(ex), ncol(clin_mm)))
if (!has_clin)
  cat("\n*** no clinical covariates. The clinical-only rung will be skipped,\n",
      "    which removes the comparison the project is built around.\n", sep = "")

## ---- randomised arm, in every rung -----------------------------------
# The ladder in the proposal starts from arm alone and adds one layer at
# a time, so arm belongs in each rung rather than only in the later ones.
arm_mm <- matrix(as.integer(arm == "treated"), ncol = 1,
                 dimnames = list(ids, "arm_treated"))
has_arm <- length(unique(na.omit(arm))) == 2
if (!has_arm) cat("single-arm cohort: the arm rung is dropped\n")

## ---- published signature scores, if script 06 has already run --------
# Rung 3 of the ladder. Script 06 depends only on 02, so run it BEFORE
# this script if you want this rung — otherwise it is skipped and the
# rest of the ladder is unaffected.
sig_mm <- NULL
sig_path <- "out/06_signature_scores.csv"
if (file.exists(sig_path)) {
  S <- data.table::fread(sig_path, data.table = FALSE)
  rownames(S) <- S$sample
  S <- S[ids, setdiff(names(S), c("sample", "y", "arm")), drop = FALSE]
  if (ncol(S)) {
    sig_mm <- as.matrix(S)
    cat(sprintf("signature scores: %d loaded from %s\n", ncol(sig_mm), sig_path))
  }
} else {
  cat("no out/06_signature_scores.csv yet — the arm+signature rung is skipped.\n")
  cat("  Run 06_signatures_interaction.R first if you want it.\n")
}

## ---- the deliberate leak, computed once on everything ----------------
sel_leaky <- names(sort(abs(col_t(ex, y)), decreasing = TRUE))[seq_len(min(TOP_K, ncol(ex)))]

## ---- stratified folds -------------------------------------------------
strat_folds <- function(y, k) {
  f <- integer(length(y))
  for (cl in unique(y)) {
    idx <- which(y == cl)
    f[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  f
}

preds <- list(); stab <- list(); t0 <- Sys.time()

for (rep_i in seq_len(N_REPEATS)) {
  folds <- strat_folds(y, N_FOLDS)
  for (k in seq_len(N_FOLDS)) {
    tr <- which(folds != k); te <- which(folds == k)
    ytr <- y[tr]; yte <- y[te]
    
    ## ---------- everything from here is fit on TRAINING ROWS ONLY -----
    # clinical: impute with training medians
    if (has_clin) {
      med <- impute_fit(clin_mm[tr, , drop = FALSE])
      Ctr <- impute_apply(clin_mm[tr, , drop = FALSE], med)
      Cte <- impute_apply(clin_mm[te, , drop = FALSE], med)
      cm  <- colMeans(Ctr); cs <- row_sd(t(Ctr)); cs[cs <= 0] <- 1
      Ctr <- scale(Ctr, cm, cs); Cte <- scale(Cte, cm, cs)
    } else {
      Ctr <- matrix(numeric(0), nrow = length(tr), ncol = 0)
      Cte <- matrix(numeric(0), nrow = length(te), ncol = 0)
    }
    
    # univariate gene filter, honest version
    tstat <- abs(col_t(ex[tr, , drop = FALSE], ytr))
    sel   <- names(sort(tstat, decreasing = TRUE))[seq_len(min(TOP_K, length(tstat)))]
    
    scale_pair <- function(cols) {
      mu  <- colMeans(ex[tr, cols, drop = FALSE])
      sdv <- matrixStats::colSds(ex[tr, cols, drop = FALSE]); sdv[sdv <= 0] <- 1
      list(tr = scale(ex[tr, cols, drop = FALSE], mu, sdv),
           te = scale(ex[te, cols, drop = FALSE], mu, sdv))
    }
    S  <- scale_pair(sel)
    SL <- scale_pair(sel_leaky)
    
    Atr <- arm_mm[tr, , drop = FALSE]; Ate <- arm_mm[te, , drop = FALSE]
    
    fp <- data.table(rep = rep_i, fold = k, sample = ids[te], y = yte, arm = arm[te])
    
    # A small helper: logistic regression on a design matrix, so the three
    # low-dimensional rungs read the same way.
    glm_rung <- function(Xtr, Xte) {
      if (is.null(Xtr) || ncol(Xtr) == 0) return(rep(NA_real_, nrow(Xte)))
      d <- as.data.frame(Xtr); d$outcome_ <- ytr
      m <- suppressWarnings(glm(outcome_ ~ ., data = d, family = binomial))
      as.numeric(predict(m, newdata = as.data.frame(Xte), type = "response"))
    }
    
    ## rung 0 — arm alone, the intercept-plus-treatment-effect model
    fp$arm_only <- if (has_arm) glm_rung(Atr, Ate) else NA_real_
    
    ## rung 1 — arm plus clinical, the bar everything else has to clear
    if (has_clin) {
      fp$clinical <- glm_rung(if (has_arm) cbind(Atr, Ctr) else Ctr,
                              if (has_arm) cbind(Ate, Cte) else Cte)
    } else fp$clinical <- NA_real_
    
    ## rung 1b — arm plus published signature scores
    if (!is.null(sig_mm)) {
      smed <- impute_fit(sig_mm[tr, , drop = FALSE])
      Str  <- impute_apply(sig_mm[tr, , drop = FALSE], smed)
      Ste  <- impute_apply(sig_mm[te, , drop = FALSE], smed)
      smu  <- colMeans(Str); ssd <- matrixStats::colSds(Str); ssd[ssd <= 0] <- 1
      Str  <- scale(Str, smu, ssd); Ste <- scale(Ste, smu, ssd)
      fp$arm_sig <- glm_rung(if (has_arm) cbind(Atr, Str) else Str,
                             if (has_arm) cbind(Ate, Ste) else Ste)
    } else fp$arm_sig <- NA_real_
    
    ## rung 2 — elastic net on expression alone
    cv2 <- cv.glmnet(S$tr, ytr, family = "binomial", alpha = ALPHA, nfolds = 5)
    fp$enet_expr <- as.numeric(predict(cv2, S$te, s = "lambda.min", type = "response"))
    
    ## rung 3 — elastic net on clinical + expression
    XtrC <- if (has_arm) cbind(Atr, Ctr, S$tr) else cbind(Ctr, S$tr)
    XteC <- if (has_arm) cbind(Ate, Cte, S$te) else cbind(Cte, S$te)
    cv3  <- cv.glmnet(XtrC, ytr, family = "binomial", alpha = ALPHA,
                      nfolds = 5, keep = TRUE)
    fp$enet_both <- as.numeric(predict(cv3, XteC, s = "lambda.min", type = "response"))
    
    ## rung 3c — the same model, recalibrated
    # cv.glmnet(keep = TRUE) returns prevalidated (inner out-of-fold) linear
    # predictors for free, so the two-parameter Platt map is fitted without
    # ever seeing the test fold. Two parameters is about all ~60 events will
    # support; isotonic regression needs far more and was dropped for that
    # reason, not forgotten.
    li  <- which(cv3$lambda == cv3$lambda.min)[1]
    pre <- cv3$fit.preval[, li]
    lte <- as.numeric(predict(cv3, XteC, s = "lambda.min", type = "link"))
    if (all(is.finite(pre)) && length(unique(pre)) > 2) {
      pl <- suppressWarnings(glm(ytr ~ pre, family = binomial))
      fp$enet_both_cal <- as.numeric(predict(pl, newdata = data.frame(pre = lte),
                                             type = "response"))
    } else fp$enet_both_cal <- fp$enet_both
    
    ## rung 4 — random forest, class-weighted
    w  <- ifelse(ytr == 1, sum(ytr == 0) / sum(ytr == 1), 1)
    rf <- ranger(x = as.data.frame(XtrC), y = factor(ytr), probability = TRUE,
                 num.trees = 500, case.weights = w, min.node.size = 5,
                 num.threads = NTHREAD, seed = 20261017 + 100 * rep_i + k)
    fp$rf <- predict(rf, as.data.frame(XteC), num.threads = NTHREAD)$predictions[, "1"]
    
    ## rung 5 — the deliberate leak, for the optimism estimate only
    cvL <- cv.glmnet(SL$tr, ytr, family = "binomial", alpha = ALPHA, nfolds = 5)
    fp$enet_leaky <- as.numeric(predict(cvL, SL$te, s = "lambda.min", type = "response"))
    
    preds[[length(preds) + 1]] <- fp
    
    cf   <- coef(cv3, s = "lambda.min")
    kept <- rownames(cf)[which(as.numeric(cf) != 0)]
    stab[[length(stab) + 1]] <- setdiff(kept, c("(Intercept)", colnames(Ctr)))
    
    cat(sprintf("repeat %d fold %d  (%s elapsed)\n", rep_i, k,
                format(round(difftime(Sys.time(), t0, units = "mins"), 1))))
  }
}

pred_dt <- rbindlist(preds)
fwrite(pred_dt, "out/04_cv_predictions.csv")

n_fits   <- N_FOLDS * N_REPEATS
stab_tab <- sort(table(unlist(stab)), decreasing = TRUE) / n_fits
fwrite(data.table(gene = names(stab_tab), selection_frequency = as.numeric(stab_tab)),
       "out/04_stability.csv")

cat(sprintf("\ntotal time %s\n", format(round(difftime(Sys.time(), t0, units = "mins"), 1))))
cat("wrote out/04_cv_predictions.csv and out/04_stability.csv\n")
cat(sprintf("genes selected in >=80%% of folds: %d\n", sum(stab_tab >= 0.8)))
cat("NEXT: 05_evaluate.R\n")
