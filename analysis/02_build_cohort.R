# =====================================================================
# 02_build_cohort.R
#
# Metadata -> outcome (0/1), arm (placebo/treated), clinical covariates.
# Produces the single object every later script reads.
#
# IN   out/01_expr_<GSE>.rds, out/01_pheno_<GSE>.csv, out/01_clinical_<GSE>.csv
# OUT  out/02_cohort_<GSE>.rds
#      out/02_exclusions_<GSE>.csv   CONSORT-style record of what was dropped
#
# The GSE206285 mapping below is already filled in from the parsed field
# names, so in the normal case nothing here needs editing. If a field name
# has changed, the script stops and prints the fields that do exist.
# =====================================================================

suppressPackageStartupMessages({ library(data.table) })
source("_helpers.R")

GSE_ID <- "GSE206285"

## ---- cohort specification --------------------------------------------
SPEC <- list(
  GSE206285 = list(
    outcome_col   = "clinical.remission.at.week.8",
    outcome_yes   = c("y", "yes", "1", "remission", "responder"),
    outcome_no    = c("n", "no", "0", "non-responder", "nonresponder"),
    second_col    = "mucosal.healing.at.week.8",   # for the endpoint audit in 06
    arm_col       = "treatment",
    placebo_regex = "^placebo",
    treated_regex = "^ustekinumab",
    drop_col      = "diagnosis",
    drop_regex    = "healthy|control",              # healthy controls are not randomised
    clin_num      = c("CRP", "Calprotectin", "Mayo_score", "Disease_duration",
                      "Neutrophils_LP", "Neutrophils_Epithelium"),
    clin_yn       = c("Anti_TNF_history"),
    clin_log2     = c("CRP", "Calprotectin")        # both are strongly right-skewed
  ),
  # Replication cohort. CONFIRM these against out/01_pheno_GSE207022.csv
  # before trusting anything script 07 prints.
  GSE207022 = list(
    outcome_col   = "response",
    outcome_yes   = c("r", "responder", "yes", "y", "1"),
    outcome_no    = c("nr", "non-responder", "no", "n", "0"),
    second_col    = NA,
    arm_col       = NA,                             # single-arm: no interaction test
    placebo_regex = NA, treated_regex = NA,
    drop_col      = NA,  drop_regex = NA,
    clin_num      = character(0), clin_yn = character(0), clin_log2 = character(0)
  )
)
sp <- SPEC[[GSE_ID]]
if (is.null(sp)) stop("No specification for ", GSE_ID, ". Add one to SPEC.")

## ---- load -------------------------------------------------------------
ex    <- readRDS(file.path("out", paste0("01_expr_", GSE_ID, ".rds")))
pheno <- fread(file.path("out", paste0("01_pheno_", GSE_ID, ".csv")),
               data.table = FALSE, colClasses = "character")
rownames(pheno) <- pheno$gsm
pheno <- pheno[colnames(ex), , drop = FALSE]
stopifnot(identical(rownames(pheno), colnames(ex)))

need <- function(col) {
  if (is.na(col)) return(invisible(NULL))
  if (!col %in% colnames(pheno)) {
    cat("\nfields that DO exist:\n"); print(colnames(pheno))
    stop("field '", col, "' is not in the parsed metadata. ",
         "Pick the right name from the list above and edit SPEC.")
  }
}
invisible(lapply(c(sp$outcome_col, sp$arm_col, sp$drop_col), need))

log_ex <- excl_new()
n0 <- ncol(ex)

## ---- outcome ----------------------------------------------------------
norm <- function(v) {
  v <- tolower(trimws(as.character(v)))
  sub("^.*:\\s*", "", v)                 # survives "remission at week 8: y"
}
ov <- norm(pheno[[sp$outcome_col]])
y  <- ifelse(ov %in% sp$outcome_yes, 1L, ifelse(ov %in% sp$outcome_no, 0L, NA_integer_))
cat("outcome field '", sp$outcome_col, "' raw values:\n", sep = "")
print(table(ov, useNA = "ifany"))
if (all(is.na(y)))
  stop("every outcome value became NA. The values above do not match ",
       "outcome_yes/outcome_no in SPEC. Edit those two vectors.")

y2 <- if (!is.na(sp$second_col) && sp$second_col %in% colnames(pheno)) {
  v <- norm(pheno[[sp$second_col]])
  ifelse(v %in% sp$outcome_yes, 1L, ifelse(v %in% sp$outcome_no, 0L, NA_integer_))
} else rep(NA_integer_, nrow(pheno))

## ---- treatment arm ----------------------------------------------------
# Anything that is neither clearly placebo nor clearly treated becomes NA
# and is dropped. The common bug here is `ifelse(is_placebo, "placebo",
# "treated")`, which quietly relabels healthy controls as treated patients.
if (is.na(sp$arm_col)) {
  arm <- rep(NA_character_, nrow(pheno))
  cat("\nno arm field for this cohort: single-arm, no interaction test\n")
} else {
  av  <- norm(pheno[[sp$arm_col]])
  arm <- ifelse(grepl(sp$placebo_regex, av), "placebo",
         ifelse(grepl(sp$treated_regex, av), "treated", NA_character_))
  cat("\narm field '", sp$arm_col, "' mapping:\n", sep = "")
  print(table(raw = av, arm = arm, useNA = "ifany"))
}

## ---- explicit exclusions ----------------------------------------------
drop <- rep(FALSE, nrow(pheno))
if (!is.na(sp$drop_col) && sp$drop_col %in% colnames(pheno)) {
  drop <- grepl(sp$drop_regex, norm(pheno[[sp$drop_col]]))
  cat("\nexcluded by '", sp$drop_col, "' ~ /", sp$drop_regex, "/: ",
      sum(drop), " samples\n", sep = "")
}
keep <- !drop
log_ex <- excl_add(log_ex, "non-randomised samples", n0, sum(keep),
                   paste0(sp$drop_col, " matches /", sp$drop_regex, "/"))

n1 <- sum(keep); keep <- keep & !is.na(y)
log_ex <- excl_add(log_ex, "outcome missing", n1, sum(keep),
                   paste0(sp$outcome_col, " not Y/N"))

if (!is.na(sp$arm_col)) {
  n2 <- sum(keep); keep <- keep & !is.na(arm)
  log_ex <- excl_add(log_ex, "arm unassignable", n2, sum(keep),
                     "treatment field matched neither arm")
}

## ---- clinical covariates ----------------------------------------------
# For GSE206285 these live in a SEPARATE supplementary file, not in the
# series matrix. Without them the "clinical baseline" the transcriptome is
# supposed to beat collapses to age and sex, and the whole comparison
# stops meaning anything.
clin <- data.frame(row.names = pheno$gsm)
clin_path <- file.path("out", paste0("01_clinical_", GSE_ID, ".csv"))
if (file.exists(clin_path)) {
  add <- fread(clin_path, data.table = FALSE, colClasses = "character")
  if (!"gsm" %in% names(add)) stop("clinical table has no 'gsm' column")
  add <- add[!duplicated(add$gsm), ]
  i   <- match(pheno$gsm, add$gsm)

  for (v in intersect(sp$clin_num, names(add)))
    clin[[v]] <- suppressWarnings(as.numeric(add[[v]][i]))
  for (v in intersect(sp$clin_yn, names(add))) {
    z <- tolower(trimws(add[[v]][i]))
    clin[[v]] <- ifelse(z %in% c("y", "yes", "1"), 1L,
                 ifelse(z %in% c("n", "no", "0"), 0L, NA_integer_))
  }
  for (v in intersect(sp$clin_log2, names(clin))) {
    clin[[paste0(v, "_log2")]] <- log2(clin[[v]] + 1)
    clin[[v]] <- NULL
  }
  cat("\nclinical covariates assembled:\n")
  for (v in names(clin))
    cat(sprintf("  %-26s %3d missing (%.0f%%)\n", v, sum(is.na(clin[[v]])),
                100 * mean(is.na(clin[[v]]))))
  cat("\nMissing values are kept as NA on purpose. They are imputed inside\n")
  cat("each training fold in script 04, never on the full dataset.\n")
} else {
  cat("\n*** no clinical table at", clin_path, "\n")
  cat("    The clinical-only model cannot be built, so the central\n")
  cat("    comparison of the project is unavailable. Fix this before\n")
  cat("    running 04. See the note at the end of script 01.\n")
}

## ---- apply the filter and align everything ---------------------------
ex   <- ex[, keep, drop = FALSE]
y    <- y[keep]; y2 <- y2[keep]; arm <- arm[keep]
clin <- clin[keep, , drop = FALSE]
ids  <- colnames(ex)
stopifnot(length(y) == ncol(ex), nrow(clin) == ncol(ex),
          identical(rownames(clin), ids))

## ---- outcome-blind gene filter ---------------------------------------
# Low-variance genes only slow the models down. This uses no outcome
# information, so it is safe to do once on the whole cohort.
#
# BUT the unfiltered matrix is kept as well, and it is the one signature
# scoring and cell-type deconvolution must use. singscore, ssGSEA and
# xCell rank a patient's genes against each other, so dropping a quarter
# of the genome changes every rank. xCell in particular returns numbers
# that look plausible and mean nothing on a filtered matrix. None of
# these methods touches the outcome or other patients' labels, so
# computing them once on all samples leaks nothing.
ex_all <- ex
v  <- row_var(ex)
ex <- ex[v > quantile(v, 0.25), , drop = FALSE]

## ---- report -----------------------------------------------------------
cat(sprintf("\n%s cohort: %d patients\n", GSE_ID, ncol(ex)))
cat(sprintf("  expr     %5d genes (variance-filtered) — for modelling\n", nrow(ex)))
cat(sprintf("  expr_all %5d genes (unfiltered)        — for signature and deconvolution scoring\n",
            nrow(ex_all)))
cat(sprintf("events: %d of %d (%.1f%%)\n", sum(y), length(y), 100 * mean(y)))
if (!is.na(sp$arm_col)) { cat("\noutcome by arm:\n"); print(table(arm = arm, outcome = y)) }
if (!all(is.na(y2)))    { cat("\nsecondary endpoint by arm:\n"); print(table(arm = arm, y2 = y2, useNA = "ifany")) }
cat("\nexclusions:\n"); print(log_ex, row.names = FALSE)

if (sum(y) < 20)
  cat("\n*** fewer than 20 events. Elastic net and random forest will be\n",
      "    very unstable. Report calibration and interval width, not AUROC\n",
      "    point estimates.\n", sep = "")

fwrite(log_ex, file.path("out", paste0("02_exclusions_", GSE_ID, ".csv")))
saveRDS(list(expr = ex, expr_all = ex_all, y = y, y2 = y2, arm = arm, clin = clin,
             ids = ids, gse = GSE_ID, built = Sys.time()),
        file.path("out", paste0("02_cohort_", GSE_ID, ".rds")))
cat("\nwrote out/02_cohort_", GSE_ID, ".rds\n", sep = "")
cat("NEXT: 03_qc_batch.R\n")
co <- readRDS("out/02_cohort_GSE206285.rds")
colSums(is.na(co$clin))

