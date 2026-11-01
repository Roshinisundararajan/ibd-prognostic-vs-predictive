# =====================================================================
# 08_crossdrug_final.R — objective (iii), final run.
#
# Changes from the first pass:
#   - 10,000 permutations, not 2,000. At 2,000 the p-value moved between
#     0.046 and 0.059 across runs, which is Monte Carlo noise sitting on
#     the decision boundary.
#   - a Hallmark pathway-level version beside the gene-level one. This is
#     the enrichment half of objective (iii): 50 pathway scores instead of
#     17,480 genes, so the multiple-testing space is small enough that a
#     real effect could survive it.
#   - scan day carried through, because in this cohort it is confounded
#     with both drug and response.
#
# Statistic: ordinary t. limma is unavailable here (no Bioconductor), so
# the moderated version must be re-run on the laptop.
#
# IN   out/02_cohort_GSE73661.rds, data/genesets/hallmark.gmt
# OUT  out/08_crossdrug_summary.csv, out/08_gene_stats_{un,}adjusted.csv,
#      out/08_hallmark_stats.csv, out/08_confounding_within_drug.csv
# =====================================================================
suppressPackageStartupMessages({ library(matrixStats); library(data.table) })
set.seed(20261020)
N_PERM <- 10000

co <- readRDS("out/02_cohort_GSE73661.rds")
keep <- co$drug %in% c("IFX", "VDZ")
ex <- co$expr[, keep, drop = FALSE]; ex_all <- co$expr_all[, keep, drop = FALSE]
y  <- co$y[keep]; drug <- co$drug[keep]; sday <- co$scan_day[keep]

col_t <- function(X, y) {
  g1 <- y == 1; n1 <- sum(g1); n0 <- sum(!g1)
  X1 <- X[g1, , drop = FALSE]; X0 <- X[!g1, , drop = FALSE]
  m1 <- colMeans(X1); m0 <- colMeans(X0)
  v1 <- colVars(X1); v0 <- colVars(X0)
  se <- sqrt(((n1 - 1) * v1 + (n0 - 1) * v0) / (n1 + n0 - 2)) * sqrt(1 / n1 + 1 / n0)
  se[!is.na(se) & se <= 0] <- NA_real_
  (m1 - m0) / se
}
resid_on <- function(X, f) {
  f <- droplevels(factor(f)); if (nlevels(f) < 2) return(X)
  M <- model.matrix(~ f); X - M %*% qr.solve(crossprod(M), crossprod(M, X))
}

## ---- confounding, within each drug ------------------------------------
conf <- rbindlist(lapply(c("IFX", "VDZ"), function(d) {
  i <- drug == d
  data.table(drug = d, n = sum(i), responders = sum(y[i]),
             scan_days = length(unique(sday[i])),
             fisher_p = fisher.test(table(sday[i], y[i]), simulate.p.value = TRUE,
                                    B = 20000)$p.value,
             r2_response_on_scan_day = summary(lm(y[i] ~ factor(sday[i])))$r.squared,
             residual_df = sum(i) - length(unique(sday[i])) - 1)
}))
fwrite(conf, "out/08_confounding_within_drug.csv"); print(conf)

## ---- one engine for gene level and pathway level ----------------------
run_pair <- function(M, adjust, level) {       # M: features x samples
  parts <- lapply(c("IFX", "VDZ"), function(d) {
    i <- which(drug == d); X <- t(M[, i, drop = FALSE])
    if (adjust) X <- resid_on(X, sday[i])
    list(X = X, y = y[i], nd = length(unique(sday[i])))
  })
  t1 <- col_t(parts[[1]]$X, parts[[1]]$y); t2 <- col_t(parts[[2]]$X, parts[[2]]$y)
  rho  <- cor(t1, t2, method = "spearman", use = "complete.obs")
  perm <- replicate(N_PERM, cor(col_t(parts[[1]]$X, sample(parts[[1]]$y)),
                                col_t(parts[[2]]$X, sample(parts[[2]]$y)),
                                method = "spearman", use = "complete.obs"))
  df1 <- length(parts[[1]]$y) - 2 - adjust * (parts[[1]]$nd - 1)
  df2 <- length(parts[[2]]$y) - 2 - adjust * (parts[[2]]$nd - 1)
  q1 <- p.adjust(2 * pt(-abs(t1), df1), "BH"); q2 <- p.adjust(2 * pt(-abs(t2), df2), "BH")
  list(row = data.table(level = level,
         model = if (adjust) "adjusted for scan day" else "unadjusted",
         features = nrow(M), n_IFX = length(parts[[1]]$y), n_VDZ = length(parts[[2]]$y),
         events_IFX = sum(parts[[1]]$y), events_VDZ = sum(parts[[2]]$y),
         spearman_rho = rho,
         null_lo = unname(quantile(perm, .025)), null_hi = unname(quantile(perm, .975)),
         n_perm = N_PERM,
         p_permutation = (1 + sum(abs(perm) >= abs(rho))) / (N_PERM + 1),
         sig_IFX_q10 = sum(q1 < .1, na.rm = TRUE), sig_VDZ_q10 = sum(q2 < .1, na.rm = TRUE),
         shared_q10 = sum(q1 < .1 & q2 < .1 & sign(t1) == sign(t2), na.rm = TRUE)),
       stats = data.table(feature = rownames(M), t_IFX = t1, q_IFX = q1,
                          t_VDZ = t2, q_VDZ = q2))
}

## ---- gene level --------------------------------------------------------
g_un <- run_pair(ex, FALSE, "gene"); g_ad <- run_pair(ex, TRUE, "gene")
fwrite(g_un$stats, "out/08_gene_stats_unadjusted.csv")
fwrite(g_ad$stats, "out/08_gene_stats_adjusted.csv")
hit <- g_un$stats$q_IFX < .1 & !is.na(g_un$stats$q_IFX)
cat(sprintf("\nIFX genes q<0.10 unadjusted: %d | still q<0.10 after adjustment: %d\n",
            sum(hit), sum(g_ad$stats$q_IFX[hit] < .1, na.rm = TRUE)))

## ---- Hallmark pathway level -------------------------------------------
read_gmt <- function(path) {
  ln <- strsplit(readLines(path), "\t")
  setNames(lapply(ln, function(x) toupper(x[-(1:2)])), vapply(ln, `[`, "", 1))
}
sets <- lapply(read_gmt("data/genesets/hallmark.gmt"), intersect, rownames(ex_all))
sets <- sets[lengths(sets) >= 5]
R <- apply(ex_all, 2, function(z) rank(z, ties.method = "average") / length(z))
H <- t(vapply(sets, function(g) colMeans(R[g, , drop = FALSE]), numeric(ncol(ex_all))))
rownames(H) <- sub("^HALLMARK_", "", names(sets))
cat(sprintf("Hallmark sets scored: %d\n", nrow(H)))

h_un <- run_pair(H, FALSE, "hallmark"); h_ad <- run_pair(H, TRUE, "hallmark")
fwrite(merge(h_un$stats, h_ad$stats, by = "feature", suffixes = c("_unadj", "_adj")),
       "out/08_hallmark_stats.csv")

summ <- rbindlist(list(g_un$row, g_ad$row, h_un$row, h_ad$row))
fwrite(summ, "out/08_crossdrug_summary.csv")
print(summ[, .(level, model, features, spearman_rho = round(spearman_rho, 3),
               null_lo = round(null_lo, 3), null_hi = round(null_hi, 3),
               p = round(p_permutation, 4), sig_IFX_q10, sig_VDZ_q10)])

cat("\ntop Hallmark pathways by each drug, adjusted:\n")
hs <- h_ad$stats[order(q_IFX)][1:5]; print(hs[, .(feature, t_IFX, q_IFX, t_VDZ, q_VDZ)])
cat("\nwrote out/08_crossdrug_summary.csv, gene stats, hallmark stats\n")
