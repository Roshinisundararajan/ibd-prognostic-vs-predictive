# =====================================================================
# 06b_hallmark_interaction.R
#
# The same prognostic-versus-predictive test as script 06, run on all 50
# MSigDB Hallmark pathways instead of the published signatures.
#
# Why this is legitimate, and how to report it
# --------------------------------------------
# Hallmark pathways are already the pre-specified primary representation
# in script 04b, so testing the same 50 pathways for a treatment
# interaction is coherent. ALL 50 are tested and BH-corrected together,
# so nothing is cherry-picked.
#
# But it is a PATHWAY-LEVEL analysis, not a test of the published
# signatures (West, Martin, Czarnewski, Linggi). Report it as such.
# When those gene lists are curated, script 06 runs on them separately.
#
# Outputs are named 06b_* on purpose, so script 04 does NOT pick these
# scores up as its "published signature" rung.
#
# Endpoint: clinical remission only. Script 06 PART A showed the arms
# separate on remission but not on mucosal healing, so the predictive
# question is only well posed on remission.
#
# IN   out/02_cohort_GSE206285.rds, data/genesets/hallmark.gmt
# OUT  out/06b_hallmark_scores.csv, out/06b_hallmark_interaction.csv
# Runtime: a few minutes (100 Firth fits).
# =====================================================================

suppressPackageStartupMessages({ library(data.table) })
source("_helpers.R")

GSE_ID <- "GSE206285"
GMT    <- "data/genesets/hallmark.gmt"
Q      <- 0.10                      # BH threshold, fixed in advance

co  <- readRDS(file.path("out", paste0("02_cohort_", GSE_ID, ".rds")))
y   <- co$y
ex  <- if (!is.null(co$expr_all)) co$expr_all else co$expr   # unfiltered
arm <- factor(co$arm, levels = c("placebo", "treated"))

cat(sprintf("%d patients, %d events | placebo %d (%d events) | treated %d (%d events)\n",
            length(y), sum(y),
            sum(arm == "placebo"), sum(y[arm == "placebo"]),
            sum(arm == "treated"), sum(y[arm == "treated"])))

## ---- gene sets ---------------------------------------------------------
if (!file.exists(GMT)) stop(GMT, " not found. Run the msigdbr step first.")
read_gmt <- function(path) {
  ln <- strsplit(readLines(path, warn = FALSE), "\t")
  ln <- ln[lengths(ln) > 3]
  setNames(lapply(ln, function(z) toupper(trimws(z[-c(1, 2)]))),
           vapply(ln, function(z) z[1], character(1)))
}
sets <- read_gmt(GMT)
sets <- lapply(sets, function(g) intersect(g, rownames(ex)))
sets <- sets[lengths(sets) >= 5]
cat(sprintf("Hallmark sets usable (>= 5 genes on the array): %d\n", length(sets)))

## ---- score each pathway per patient (rank-based, outcome-blind) --------
if (!requireNamespace("singscore", quietly = TRUE))
  stop("singscore is not installed: BiocManager::install('singscore')")
rk <- singscore::rankGenes(ex)                      # rank once, reuse
scores <- vapply(sets, function(g)
  as.numeric(singscore::simpleScore(rk, upSet = g)$TotalScore),
  numeric(ncol(ex)))
colnames(scores) <- sub("^HALLMARK_", "", names(sets))
scores <- scale(scores)                              # ORs are per SD
fwrite(data.table(sample = co$ids, y = y, arm = as.character(arm), scores),
       "out/06b_hallmark_scores.csv")

## ---- interaction model ---------------------------------------------------
# Firth-penalised: ~60 events split across two arms plus an interaction
# term pushes ordinary maximum likelihood into separation.
has_firth <- requireNamespace("logistf", quietly = TRUE)
if (!has_firth)
  cat("*** logistf not installed, falling back to glm. Do not report that.\n")

fit_interaction <- function(dat, extra = character(0)) {
  f <- as.formula(paste("y ~ arm * sc", if (length(extra))
    paste("+", paste(extra, collapse = " + ")) else ""))
  if (has_firth) {
    m <- logistf::logistf(f, data = dat)
    k <- grep(":", names(coef(m)))[1]; j <- which(names(coef(m)) == "sc")
    list(i_or = exp(coef(m)[k]), i_lo = exp(m$ci.lower[k]), i_hi = exp(m$ci.upper[k]),
         i_p = m$prob[k], m_or = exp(coef(m)[j]), m_p = m$prob[j])
  } else {
    m <- suppressWarnings(glm(f, data = dat, family = binomial))
    s <- summary(m)$coefficients; ci <- suppressWarnings(confint.default(m))
    k <- grep(":", rownames(s))[1]; j <- which(rownames(s) == "sc")
    list(i_or = exp(s[k, 1]), i_lo = exp(ci[k, 1]), i_hi = exp(ci[k, 2]),
         i_p = s[k, 4], m_or = exp(s[j, 1]), m_p = s[j, 4])
  }
}

adj <- names(co$clin)[vapply(co$clin, function(z) mean(is.na(z)) < 0.2, logical(1))]
adj <- intersect(adj, c("CRP_log2", "Calprotectin_log2", "Mayo_score",
                        "Disease_duration", "Anti_TNF_history"))
cat("adjusted-model covariates:", paste(adj, collapse = ", "), "\n\n")

t0 <- Sys.time()
res <- rbindlist(lapply(seq_len(ncol(scores)), function(i) {
  s  <- colnames(scores)[i]
  sc <- scores[, i]
  dat <- data.frame(y = y, arm = arm, sc = sc, co$clin[, adj, drop = FALSE])
  u <- fit_interaction(dat)
  a <- tryCatch(fit_interaction(dat[stats::complete.cases(dat), ], adj),
                error = function(e) list(i_or = NA, i_lo = NA, i_hi = NA, i_p = NA))
  if (i %% 10 == 0) cat(sprintf("  %d of %d pathways (%s)\n", i, ncol(scores),
                                format(round(difftime(Sys.time(), t0, units = "mins"), 1))))
  # placebo is the reference level of arm, so `sc` = effect IN THE PLACEBO ARM
  data.table(pathway = s,
             OR_placebo_arm = u$m_or, p_placebo_arm = u$m_p,
             inter_OR = u$i_or, inter_lo = u$i_lo, inter_hi = u$i_hi, inter_p = u$i_p,
             inter_OR_adj = a$i_or, inter_lo_adj = a$i_lo,
             inter_hi_adj = a$i_hi, inter_p_adj = a$i_p)
}))

## ---- classification, fixed in advance: BH q < 0.10 on each axis ---------
res[, prog_q  := p.adjust(p_placebo_arm, method = "BH")]
res[, inter_q := p.adjust(inter_p,       method = "BH")]
res[, call := fifelse(inter_q < Q & prog_q < Q, "prognostic + predictive",
                      fifelse(inter_q < Q,               "predictive only",
                              fifelse(prog_q  < Q,               "prognostic only", "null")))]
setorder(res, inter_q, prog_q)
fwrite(res, "out/06b_hallmark_interaction.csv")

cat("\nHow many pathways fall in each cell:\n")
print(table(res$call))

cat("\nTop 15 by interaction q:\n")
print(as.data.frame(res[1:min(15, .N),
                        .(pathway, OR_placebo_arm, prog_q, inter_OR, inter_lo, inter_hi, inter_q, call)]),
      row.names = FALSE, digits = 3)

se  <- (log(res$inter_hi) - log(res$inter_lo)) / (2 * 1.96)
med <- median(se, na.rm = TRUE)
cat(sprintf("\nmedian SE of the log interaction OR: %.3f\n", med))
cat(sprintf("smallest interaction OR detectable (80%% power, two-sided 0.05): %.2f per SD\n",
            exp(2.802 * med)))

cat("\nNote: the placebo arm has only", sum(y[arm == "placebo"]),
    "remission events, so the\nprognostic (placebo-arm) axis has little power. A 'null' call on that axis\n",
    "is weak evidence; report the odds ratio and interval, not just the label.\n")
cat("\nwrote out/06b_hallmark_scores.csv and out/06b_hallmark_interaction.csv\n")
