# =====================================================================
# 06_signatures_interaction.R
#
# The part that makes this project different from every other
# response-prediction paper.
#
# A main effect says a signature tracks who does well whatever they are
# given — prognostic. Only a treatment-by-signature interaction says the
# signature tells you WHICH drug to give — predictive. A biomarker can be
# strongly significant in the treated arm and still be useless for drug
# choice, because the placebo arm does the same thing. That mistake has
# its own paper: Altman & Bland, BMJ 2003 — a difference in significance
# is not the significance of a difference.
#
# Signature input: data/signatures/
#     <name>.txt                  one gene symbol per line (up-regulated)
#     <name>_up.txt + <name>_dn.txt   if the published signature has a
#                                     direction, use the pair — a signature
#                                     scored as if every gene pointed the
#                                     same way is a different signature
#
# IN   out/02_cohort_<GSE>.rds, data/signatures/*.txt
# OUT  out/06_endpoint_audit.csv, out/06_signature_scores.csv,
#      out/06_interaction_results.csv
# =====================================================================

suppressPackageStartupMessages({ library(data.table) })
source("_helpers.R")

GSE_ID <- "GSE206285"
co  <- readRDS(file.path("out", paste0("02_cohort_", GSE_ID, ".rds")))
y   <- co$y

# The UNFILTERED matrix, deliberately. Rank-based scoring ranks a
# patient's genes against each other, so scoring on the variance-filtered
# matrix silently changes every rank. Same reason xCell must see the full
# matrix if you add deconvolution later.
if (!is.null(co$expr_all)) {
  ex <- co$expr_all
} else {
  ex <- co$expr
  cat("*** this cohort object has no expr_all: scoring on the variance-filtered\n",
      "    matrix, which distorts within-sample ranks. Rerun 02_build_cohort.R.\n", sep = "")
}
arm <- factor(co$arm, levels = c("placebo", "treated"))
if (all(is.na(arm))) stop("no treatment arm in this cohort: no interaction to test")

cat(sprintf("%d patients, %d events, %d placebo / %d treated\n\n",
            length(y), sum(y), sum(arm == "placebo"), sum(arm == "treated")))

## =====================================================================
## PART A — endpoint audit. Run this FIRST.
## =====================================================================
# An interaction can only exist on an endpoint where the drug has an
# average effect in the first place. If the arms do not separate on an
# endpoint, a biomarker cannot tell you who benefits on it, because
# nobody measurably benefits. Checking this before testing interactions
# is the step that keeps the negative results interpretable — and no
# published analysis of this dataset reports it.
audit <- rbindlist(lapply(list(remission = y, healing = co$y2), function(yy) {
  if (all(is.na(yy))) return(NULL)
  ok <- !is.na(yy) & !is.na(arm)
  m  <- suppressWarnings(glm(yy[ok] ~ arm[ok], family = binomial))
  ci <- suppressWarnings(confint.default(m))
  data.table(n = sum(ok), events = sum(yy[ok]),
             OR = exp(coef(m)[[2]]), lo = exp(ci[2, 1]), hi = exp(ci[2, 2]),
             p = summary(m)$coefficients[2, 4])
}), idcol = "endpoint")
if (nrow(audit)) {
  audit[, arms_separate := lo > 1 | hi < 1]
  cat("PART A — does the drug beat placebo at all, endpoint by endpoint?\n")
  # as.data.frame() on purpose: data.table suppresses the print that
  # immediately follows a `:=`, and the table silently vanishes.
  print(as.data.frame(audit), row.names = FALSE)
  cat("\nThe predictive question is only well posed on an endpoint where\n")
  cat("arms_separate is TRUE. Say so explicitly for any endpoint where it\n")
  cat("is FALSE, rather than reporting a null interaction without context.\n\n")
  fwrite(audit, "out/06_endpoint_audit.csv")
}

## =====================================================================
## PART B — score each published signature, per patient
## =====================================================================
# Scores are rank-based: each patient's genes are ranked within that
# patient, so the score does not depend on platform, batch or
# normalisation. That is what lets the same score be computed on the
# replication cohort without merging matrices or running ComBat.
use_singscore <- requireNamespace("singscore", quietly = TRUE)
cat("scoring engine:", if (use_singscore) "singscore" else
    "built-in rank score (singscore not installed)", "\n\n")

rank_score <- function(ex, up, dn = character(0)) {
  R  <- apply(ex, 2, function(z) rank(z, ties.method = "average") / length(z))
  su <- if (length(up)) colMeans(R[up, , drop = FALSE]) else rep(0.5, ncol(ex))
  if (length(dn)) su - colMeans(R[dn, , drop = FALSE]) else su - 0.5
}

sig_files <- list.files("data/signatures", pattern = "\\.txt$", full.names = TRUE)
if (!length(sig_files))
  stop("No signature files in data/signatures/. This is the long pole of the\n",
       "  project: pull the gene lists from the supplementary tables of the\n",
       "  papers you cite, one symbol per line. Nothing in PART B or C runs\n",
       "  without them.")

base <- sub("_(up|dn|down)$", "", sub("\\.txt$", "", basename(sig_files)))
genes_avail <- rownames(ex)
scores <- sapply(unique(base), function(nm) {
  f_up <- sig_files[basename(sig_files) %in% c(paste0(nm, ".txt"), paste0(nm, "_up.txt"))]
  f_dn <- sig_files[basename(sig_files) %in% c(paste0(nm, "_dn.txt"), paste0(nm, "_down.txt"))]
  up <- if (length(f_up)) toupper(trimws(readLines(f_up[1]))) else character(0)
  dn <- if (length(f_dn)) toupper(trimws(readLines(f_dn[1]))) else character(0)
  up <- up[nzchar(up)]; dn <- dn[nzchar(dn)]
  u <- intersect(up, genes_avail); d <- intersect(dn, genes_avail)
  cat(sprintf("  %-26s up %3d/%3d   down %3d/%3d\n", nm, length(u), length(up),
              length(d), length(dn)))
  if (length(u) + length(d) < 3) { cat("     too few genes found — skipped\n")
                                   return(rep(NA_real_, ncol(ex))) }
  if (use_singscore) {
    r <- singscore::rankGenes(ex)
    s <- if (length(d))
      singscore::simpleScore(r, upSet = u, downSet = d)$TotalScore
    else singscore::simpleScore(r, upSet = u)$TotalScore
    as.numeric(s)
  } else rank_score(ex, u, d)
})
scores <- scale(scores)          # so each OR is "per standard deviation"
fwrite(data.table(sample = co$ids, y = y, arm = as.character(arm), scores),
       "out/06_signature_scores.csv")

## =====================================================================
## PART C — interaction test
## =====================================================================
# Firth penalisation, not plain logistic regression. At ~60 events split
# across two arms and an interaction term, ordinary maximum likelihood
# drifts into separation and returns odds ratios in the hundreds with
# meaningless intervals. The dry run produced exactly that.
has_firth <- requireNamespace("logistf", quietly = TRUE)
if (!has_firth)
  cat("\n*** logistf is not installed, so this falls back to plain glm.\n",
      "    Install it before you report anything: install.packages('logistf')\n", sep = "")

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

# Routine clinical variables that are complete enough to adjust for.
adj <- names(co$clin)[vapply(co$clin, function(z) mean(is.na(z)) < 0.2, logical(1))]
adj <- intersect(adj, c("CRP_log2", "Calprotectin_log2", "Mayo_score",
                        "Disease_duration", "Anti_TNF_history"))
cat("\nadjusted model covariates:", if (length(adj)) paste(adj, collapse = ", ") else "none", "\n")

res <- rbindlist(lapply(colnames(scores), function(s) {
  sc <- scores[, s]
  if (all(is.na(sc))) return(NULL)
  dat <- data.frame(y = y, arm = arm, sc = sc, co$clin[, adj, drop = FALSE])
  dat <- dat[stats::complete.cases(dat[, c("y", "arm", "sc")]), ]
  u <- fit_interaction(dat)
  a <- tryCatch(fit_interaction(dat[stats::complete.cases(dat), ], adj),
                error = function(e) list(i_or = NA, i_lo = NA, i_hi = NA, i_p = NA))
  # arm is a factor with placebo as the reference level, so the `sc`
  # coefficient is the signature's effect IN THE PLACEBO ARM. That is
  # exactly what a clean prognostic estimate should be, and it is not the
  # pooled effect across both arms. Named accordingly so no reader can
  # mistake it.
  data.table(signature = s, n = nrow(dat),
             OR_placebo_arm = u$m_or, p_placebo_arm = u$m_p,
             inter_OR = u$i_or, inter_lo = u$i_lo, inter_hi = u$i_hi, inter_p = u$i_p,
             inter_OR_adj = a$i_or, inter_lo_adj = a$i_lo,
             inter_hi_adj = a$i_hi, inter_p_adj = a$i_p)
}))
if (!nrow(res)) stop("no signature produced a usable score")

## ---- classification rule, fixed in advance ---------------------------
# PRIMARY DECISION: BH-adjusted q at 0.10, on both axes. The confidence
# interval is reported beside it but does not make the call — naming one
# of them primary in advance is what stops anyone asking why you used
# whichever was kinder.
#
# Two independent yes/no questions give FOUR cells, not three. A marker
# can be prognostic in the placebo arm AND modified by treatment; forcing
# the labels to be exclusive hides that case.
Q <- 0.10
res[, inter_q := p.adjust(inter_p,      method = "BH")]
res[, prog_q  := p.adjust(p_placebo_arm, method = "BH")]
res[, prognostic := prog_q  < Q]
res[, predictive := inter_q < Q]
res[, call := fifelse( predictive &  prognostic, "prognostic + predictive",
              fifelse( predictive & !prognostic, "predictive only",
              fifelse(!predictive &  prognostic, "prognostic only", "null")))]

fwrite(res, "out/06_interaction_results.csv")
cat("\nPART C — interaction results\n")
cat("OR_placebo_arm = the signature's effect in the PLACEBO arm (prognosis).\n")
cat("inter_OR       = how that effect changes under treatment (prediction).\n")
cat(sprintf("call           = BH q < %.2f on each axis, decided in advance.\n\n", Q))
print(as.data.frame(res[, .(signature, OR_placebo_arm, prog_q, inter_OR,
                            inter_lo, inter_hi, inter_q, call)]), row.names = FALSE)

## ---- what a null result here actually rules out ----------------------
# A negative interaction is only usable if you say how large an effect the
# design could have detected. Without this number "no interaction" is
# indistinguishable from "underpowered", which is the criticism every
# negative biomarker paper gets.
se  <- (log(res$inter_hi) - log(res$inter_lo)) / (2 * 1.96)
med <- median(se, na.rm = TRUE)
cat(sprintf("\nmedian SE of the log interaction OR: %.3f\n", med))
cat(sprintf("smallest interaction OR this design could detect (80%% power, two-sided 0.05): %.2f\n",
            exp(2.802 * med)))
cat("Quote that figure next to every null call. It converts 'we found nothing'\n")
cat("into 'we ruled out effects larger than this', which is a result.\n")
cat("\nwrote out/06_endpoint_audit.csv, out/06_signature_scores.csv,\n")
cat("      out/06_interaction_results.csv\n")
