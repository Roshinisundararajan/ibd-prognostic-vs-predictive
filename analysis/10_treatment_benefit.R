# =====================================================================
# 10_treatment_benefit.R  —  causal / heterogeneous-treatment-effect ML
#
# Scripts 04-05 ask "WHO WILL BE IN REMISSION?" (a prognostic question).
# Script 06 tests one signature at a time for an arm x signature
# interaction. Neither asks the clinical question directly:
#
#     "For THIS patient, how much MORE likely is remission on
#      ustekinumab than on placebo?"        tau(x) = P(Y=1|T,x) - P(Y=1|P,x)
#
# That quantity is the individual treatment benefit, and it is what a
# biologic-selection tool has to estimate. Here it is learned with
# meta-learners, out of fold, and then CHECKED against the randomisation:
#
#   learners (each on two representations: clinical; clinical + 50 Hallmark)
#     T_enet  : separate elastic-net models in each arm, tau = mu1 - mu0
#     T_rf    : the same with probability random forests
#     S_enet  : one elastic net with arm x feature interactions
#     X_enet  : X-learner (Kunzel 2019), randomisation probability known
#
#   evaluation (all on out-of-fold tau, patients never seen by the model)
#     1. benefit by tertile of predicted benefit: observed remission rate
#        on drug minus on placebo, inside each tertile. If the model works,
#        the top tertile shows a bigger real difference than the bottom.
#     2. calibration-for-benefit: slope of the doubly robust (AIPW) effect
#        score on predicted tau. 0 = no usable heterogeneity, 1 = calibrated.
#     3. RATE / AUTOC (Yadlowsky 2021): does treating the patients ranked
#        highest by tau give more benefit than treating at random?
#        AUTOC > 0 with a CI excluding 0 = the ranking carries real
#        treatment-effect information. Bootstrap over patients.
#
# Power caveat: the placebo arm has 14 events. Report intervals, not labels.
#
# IN   out/02_cohort_GSE206285.rds, out/06b_hallmark_scores.csv
# OUT  out/10_tau_oof.csv, out/10_benefit_tertiles.csv, out/10_rate.csv,
#      out/10_toc.png
# Runtime: ~5-10 min.
# =====================================================================

suppressPackageStartupMessages({
  library(glmnet); library(ranger); library(data.table); library(matrixStats)
})
source("_helpers.R")
set.seed(20261019)

N_FOLDS <- 5; N_REPEATS <- 3; ALPHA <- 0.5; B <- 2000

co  <- readRDS("out/02_cohort_GSE206285.rds")
y   <- co$y; ids <- co$ids; W <- as.integer(co$arm == "treated")
H   <- fread("out/06b_hallmark_scores.csv", data.table = FALSE)
H   <- as.matrix(H[match(ids, H$sample), setdiff(names(H), c("sample", "y", "arm"))])
clin_mm <- build_clin_mm(co$clin)
stopifnot(nrow(H) == length(y), nrow(clin_mm) == length(y))
cat(sprintf("n %d | treated %d (%d events) | placebo %d (%d events)\n",
            length(y), sum(W), sum(y[W == 1]), sum(W == 0), sum(y[W == 0])))

reps <- list(clinical = function(tr, te) NULL,           # Hallmark block off
             clin_hallmark = function(tr, te) list(tr = H[tr, ], te = H[te, ]))

design <- function(tr, te, rep) {
  med <- impute_fit(clin_mm[tr, , drop = FALSE])
  Ctr <- impute_apply(clin_mm[tr, , drop = FALSE], med)
  Cte <- impute_apply(clin_mm[te, , drop = FALSE], med)
  cm <- colMeans(Ctr); cs <- colSds(Ctr); cs[cs <= 0] <- 1
  Xtr <- scale(Ctr, cm, cs); Xte <- scale(Cte, cm, cs)
  extra <- reps[[rep]](tr, te)
  if (!is.null(extra)) {
    mu <- colMeans(extra$tr); sd <- colSds(extra$tr); sd[sd <= 0] <- 1
    Xtr <- cbind(Xtr, scale(extra$tr, mu, sd)); Xte <- cbind(Xte, scale(extra$te, mu, sd))
  }
  list(tr = Xtr, te = Xte)
}

# elastic net that degrades to the arm mean when events are too few to tune
enet <- function(X, yy, fam = "binomial", pf = NULL) {
  if (is.null(pf)) pf <- rep(1, ncol(X))
  fit <- tryCatch(cv.glmnet(X, yy, family = fam, alpha = ALPHA, nfolds = 5,
                            penalty.factor = pf), error = function(e) NULL)
  if (is.null(fit)) { m <- mean(yy); return(function(Z) rep(m, nrow(Z))) }
  function(Z) as.numeric(predict(fit, Z, s = "lambda.min",
                                 type = if (fam == "binomial") "response" else "link"))
}
rfp <- function(X, yy) {
  m <- ranger(x = X, y = factor(yy, levels = 0:1), probability = TRUE,
              num.trees = 1000, min.node.size = 10, num.threads = 1)
  function(Z) predict(m, Z)$predictions[, "1"]
}

learn <- function(Xtr, ytr, Wtr, Xte) {
  t1 <- Wtr == 1; t0 <- Wtr == 0; e <- mean(Wtr)
  m1 <- enet(Xtr[t1, ], ytr[t1]); m0 <- enet(Xtr[t0, ], ytr[t0])
  mu1 <- m1(Xte); mu0 <- m0(Xte)
  r1 <- rfp(Xtr[t1, ], ytr[t1]); r0 <- rfp(Xtr[t0, ], ytr[t0])
  # S-learner: [X, W, (W - e) * X]; arm left unpenalised
  Sx  <- function(X, w) cbind(X, arm = w, (w - e) * X)
  sm  <- enet(Sx(Xtr, Wtr), ytr, pf = c(rep(1, ncol(Xtr)), 0, rep(1, ncol(Xtr))))
  # X-learner: imputed individual effects, regressed within each arm
  D1 <- ytr[t1] - m0(Xtr[t1, , drop = FALSE])
  D0 <- m1(Xtr[t0, , drop = FALSE]) - ytr[t0]
  g1 <- enet(Xtr[t1, ], D1, fam = "gaussian"); g0 <- enet(Xtr[t0, ], D0, fam = "gaussian")
  list(mu1 = mu1, mu0 = mu0,
       T_enet = mu1 - mu0,
       T_rf   = r1(Xte) - r0(Xte),
       S_enet = sm(Sx(Xte, 1)) - sm(Sx(Xte, 0)),
       X_enet = e * g0(Xte) + (1 - e) * g1(Xte))
}

strat2 <- function(k) {                  # stratify on arm x outcome
  g <- paste(W, y); f <- integer(length(y))
  for (s in unique(g)) { i <- which(g == s); f[i] <- sample(rep(seq_len(k), length.out = length(i))) }
  f
}

out <- list(); t0 <- Sys.time()
for (r in seq_len(N_REPEATS)) {
  f <- strat2(N_FOLDS)
  for (k in seq_len(N_FOLDS)) {
    tr <- which(f != k); te <- which(f == k)
    for (rp in names(reps)) {
      D <- design(tr, te, rp)
      L <- learn(D$tr, y[tr], W[tr], D$te)
      out[[length(out) + 1]] <- data.table(rep = r, fold = k, representation = rp,
        sample = ids[te], y = y[te], W = W[te], as.data.table(L))
    }
    cat(sprintf("repeat %d fold %d (%.1f min)\n", r, k,
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}
O <- rbindlist(out); fwrite(O, "out/10_tau_oof.csv")
learners <- c("T_enet", "T_rf", "S_enet", "X_enet")
A <- O[, lapply(.SD, mean), by = .(representation, sample, y, W),
       .SDcols = c("mu1", "mu0", learners)]

## ---- doubly robust (AIPW) effect score per patient ---------------------
# e is KNOWN (randomised), so Gamma is unbiased for tau(x) whatever mu does.
e_hat <- mean(W)
A[, Gamma := (mu1 - mu0) + W * (y - mu1) / e_hat - (1 - W) * (y - mu0) / (1 - e_hat)]
ATE <- A[representation == "clinical", mean(Gamma)]
cat(sprintf("\naverage treatment effect (AIPW): %+.3f  (raw: %+.3f)\n",
            ATE, mean(y[W == 1]) - mean(y[W == 0])))

## ---- RATE / AUTOC ---------------------------------------------------------
toc_curve <- function(tau, G, qs = seq(0.05, 1, by = 0.05)) {
  o <- order(tau, decreasing = TRUE); Gs <- G[o]; n <- length(G)
  sapply(qs, function(q) mean(Gs[seq_len(max(1, floor(q * n)))])) - mean(G)
}
autoc <- function(tau, G) {                      # average over every cut
  o <- order(tau, decreasing = TRUE); cs <- cumsum(G[o]) / seq_along(G)
  mean(cs - mean(G))
}

res <- list(); terc <- list(); curves <- list()
for (rp in names(reps)) for (l in learners) {
  d <- A[representation == rp]; tau <- d[[l]]
  est <- autoc(tau, d$Gamma)
  bs  <- replicate(B, { i <- sample.int(nrow(d), replace = TRUE); autoc(tau[i], d$Gamma[i]) })
  slope <- coef(lm(d$Gamma ~ tau))[2]
  sb  <- replicate(B, { i <- sample.int(nrow(d), replace = TRUE)
                        coef(lm(d$Gamma[i] ~ tau[i]))[2] })
  ip  <- summary(glm(y ~ W * scale(tau), family = binomial, data = d))$coefficients
  res[[length(res) + 1]] <- data.table(representation = rp, learner = l,
    sd_tau = sd(tau), AUTOC = est, AUTOC_lo = quantile(bs, .025), AUTOC_hi = quantile(bs, .975),
    p_AUTOC = 2 * min(mean(bs <= 0), mean(bs >= 0)),
    benefit_slope = slope, slope_lo = quantile(sb, .025, na.rm = TRUE),
    slope_hi = quantile(sb, .975, na.rm = TRUE),
    interaction_p = ip[nrow(ip), 4])
  tt <- cut(rank(tau, ties.method = "first"), 3, labels = c("low", "mid", "high"))
  terc[[length(terc) + 1]] <- d[, .(n = .N,
      pred_benefit = mean(get(l)),
      remit_drug = mean(y[W == 1]), n_drug = sum(W == 1),
      remit_placebo = mean(y[W == 0]), n_placebo = sum(W == 0),
      events_placebo = sum(y[W == 0]),
      observed_benefit = mean(y[W == 1]) - mean(y[W == 0]),
      se = sqrt(mean(y[W == 1]) * (1 - mean(y[W == 1])) / sum(W == 1) +
                mean(y[W == 0]) * (1 - mean(y[W == 0])) / sum(W == 0))),
      by = .(tertile = tt)][order(tertile)][, `:=`(representation = rp, learner = l)]
  curves[[paste(rp, l)]] <- toc_curve(tau, d$Gamma)
}
R  <- rbindlist(res); TT <- rbindlist(terc)
TT[, `:=`(lo = observed_benefit - 1.96 * se, hi = observed_benefit + 1.96 * se)]
fwrite(R, "out/10_rate.csv"); fwrite(TT, "out/10_benefit_tertiles.csv")

png("out/10_toc.png", 1600, 800, res = 150); par(mfrow = c(1, 2))
qs <- seq(0.05, 1, by = 0.05); cols <- c("black", "darkorange", "steelblue", "firebrick")
for (rp in names(reps)) {
  m <- sapply(learners, function(l) curves[[paste(rp, l)]])
  matplot(qs, m, type = "l", lty = 1, lwd = 2, col = cols, ylim = range(c(m, -0.1, 0.15)),
          xlab = "fraction treated (highest predicted benefit first)",
          ylab = "extra benefit vs treating at random", main = rp)
  abline(h = 0, lty = 2); legend("topright", learners, col = cols, lwd = 2, bty = "n")
}
dev.off()

cat("\n==== does the predicted benefit ranking carry real treatment-effect information? ====\n")
print(as.data.frame(R), digits = 3, row.names = FALSE)
cat("\n==== observed drug-minus-placebo remission by tertile of predicted benefit ====\n")
print(as.data.frame(TT[, .(representation, learner, tertile, n, pred_benefit,
      remit_drug, remit_placebo, events_placebo, observed_benefit, lo, hi)]),
      digits = 3, row.names = FALSE)
cat("\nwrote out/10_tau_oof.csv, out/10_rate.csv, out/10_benefit_tertiles.csv, out/10_toc.png\n")
