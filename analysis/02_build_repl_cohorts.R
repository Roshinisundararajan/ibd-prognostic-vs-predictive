# Build out/02_cohort_GSE207022.rds and out/02_cohort_GSE73661.rds in the
# format Roshini's scripts expect: expr, expr_all, y, y2, arm, clin, ids.
#
# Source of truth for outcome and arm is the baseline annotation CSV that
# already exists in UC_project/data/<GSE>/, not a re-parse of GEO.
suppressPackageStartupMessages({ library(matrixStats) })
U <- "/mnt/user-data/uploads/UC_project/data"
dir.create("out", showWarnings = FALSE)

build <- function(gse, ann, arm_vec, clin, y2 = NULL) {
  ex <- readRDS(file.path(U, gse, paste0(gse, "_expr_gene.rds")))
  keep <- ann$response %in% c("0", "1")
  ann <- ann[keep, ]; arm_vec <- arm_vec[keep]; clin <- clin[keep, , drop = FALSE]
  if (!is.null(y2)) y2 <- y2[keep]
  ex <- ex[, ann$geo_accession, drop = FALSE]
  y  <- as.integer(ann$response)

  ex_all <- ex
  v  <- rowVars(ex)
  ex <- ex[v > quantile(v, 0.25), , drop = FALSE]
  rownames(clin) <- colnames(ex)

  co <- list(expr = ex, expr_all = ex_all, y = y,
             y2 = if (is.null(y2)) rep(NA_integer_, length(y)) else y2,
             arm = arm_vec, clin = clin, ids = colnames(ex), gse = gse,
             built = Sys.time())
  saveRDS(co, file.path("out", paste0("02_cohort_", gse, ".rds")))
  cat(sprintf("\n%s: %d patients | expr %d genes | expr_all %d genes | events %d (%.1f%%)\n",
              gse, ncol(ex), nrow(ex), nrow(ex_all), sum(y), 100 * mean(y)))
  if (!all(is.na(arm_vec))) print(table(arm = arm_vec, y = y))
  co
}

## ---- GSE207022: Crohn's, ustekinumab vs placebo, randomised -----------
a <- read.csv(file.path(U, "GSE207022/GSE207022_baseline_annotation.csv"),
              stringsAsFactors = FALSE)
arm <- ifelse(a$arm == "PBO", "placebo", ifelse(a$arm == "UST", "treated", NA))
clin <- data.frame(inflamed_area = as.integer(a$inflamed_area_w0 == "Ileum and colon"),
                   dose_6mgkg    = as.integer(a$dose == "6mg/kg"))
build("GSE207022", a, arm, clin)

## ---- GSE73661: UC, infliximab and vedolizumab, no placebo induction ---
b <- read.csv(file.path(U, "GSE73661/GSE73661_baseline_annotation.csv"),
              stringsAsFactors = FALSE)
sd73 <- read.csv(file.path(U, "GSE73661/GSE73661_scan_dates.csv"), stringsAsFactors = FALSE)
b$scan_day <- sd73$scan_day[match(b$geo_accession, sd73$geo_accession)]
# 3 plac_plac patients exist but cannot serve as a randomised control arm;
# arm stays NA so no script tries an interaction test here.
clin2 <- data.frame(mayo_endo_w0 = suppressWarnings(as.numeric(b$mayo_endo_w0)))
co <- build("GSE73661", b, rep(NA_character_, nrow(b)), clin2)
# drug and scan day travel with the cohort for objective (iii)
co$drug     <- b$drug[b$response %in% c("0", "1")]
co$scan_day <- b$scan_day[b$response %in% c("0", "1")]
saveRDS(co, "out/02_cohort_GSE73661.rds")
cat("\ndrug:\n"); print(table(co$drug, co$y))

## ---- gene-symbol vocabulary overlap with discovery --------------------
d <- readRDS("out/02_cohort_GSE206285.rds")
for (g in c("GSE207022", "GSE73661")) {
  r <- readRDS(file.path("out", paste0("02_cohort_", g, ".rds")))
  cat(sprintf("\n%s shared with discovery: expr %d | expr_all %d\n", g,
              length(intersect(rownames(d$expr), rownames(r$expr))),
              length(intersect(rownames(d$expr_all), rownames(r$expr_all)))))
}
