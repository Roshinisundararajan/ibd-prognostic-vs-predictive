# =====================================================================
# 01_download_annotate.R
#
# Fetch a GEO series, parse the metadata into named columns, collapse
# probes to genes. Nothing here depends on the outcome.
#
# OUT  out/01_expr_<GSE>.rds      gene-level matrix, genes x samples
#      out/01_pheno_<GSE>.csv     tidy metadata, one row per sample
#      out/01_clinical_<GSE>.csv  supplementary clinical table, if present
#      out/01_processing_<GSE>.txt how the submitters processed it (methods section)
#
# Re-running is cheap: the series matrix is cached under data/.
# Memory: GSE206285 is ~550 samples. Give R at least 8 GB.
# =====================================================================

suppressPackageStartupMessages({
  library(GEOquery); library(Biobase); library(matrixStats); library(data.table)
})
source("_helpers.R")

GSE_ID <- "GSE206285"          # change to "GSE207022" for replication
dir.create("out",  showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)
options(timeout = 1800)

## ---- 1. series matrix -------------------------------------------------
raw_path <- file.path("data", paste0(GSE_ID, "_eset.rds"))
if (file.exists(raw_path)) {
  eset <- readRDS(raw_path)
  cat("loaded cached", raw_path, "\n")
} else {
  gse <- getGEO(GSE_ID, GSEMatrix = TRUE, AnnotGPL = FALSE, destdir = "data")
  cat("platforms in series:", length(gse), "\n")
  for (i in seq_along(gse))
    cat(sprintf("  [%d] %s  %d samples  %d features\n", i,
                annotation(gse[[i]]), ncol(gse[[i]]), nrow(gse[[i]])))
  eset <- gse[[which.max(sapply(gse, ncol))]]   # keep the largest platform
  saveRDS(eset, raw_path)
}
GPL <- annotation(eset)
pd  <- pData(eset)
cat(sprintf("\n%s on %s: %d samples\n", GSE_ID, GPL, ncol(eset)))

## ---- 2. metadata into named columns -----------------------------------
# This is the step that removes the "open the CSV and guess which column
# holds the outcome" risk. After this, fields are named for what they are.
pheno <- parse_characteristics(pd)
pheno$title <- pd$title[match(pheno$gsm, rownames(pd))]

# Technical processing-group codes the sponsor left in the CEL file names.
# Derived as metadata only; whether they matter is tested in script 03.
if ("supplementary_file" %in% colnames(pd)) {
  fn <- basename(as.character(pd$supplementary_file[match(pheno$gsm, rownames(pd))]))
  pheno$cel_file <- fn
  pheno$b_code <- ifelse(grepl("B[0-9]+",  fn), sub(".*?(B[0-9]+).*",  "\\1", fn), NA)
  pheno$p_code <- ifelse(grepl("_P[0-9]+", fn), sub(".*?(P[0-9]+).*", "\\1", fn), NA)
}

cat("\nparsed metadata fields:\n")
print(setdiff(colnames(pheno), c("gsm", "title", "cel_file")))
cat("\nvalue counts for each short field (candidates for outcome and arm):\n")
for (k in setdiff(colnames(pheno), c("gsm", "title", "cel_file"))) {
  u <- unique(na.omit(pheno[[k]]))
  if (length(u) <= 8) cat(sprintf("  %-34s %s\n", k, paste(u, collapse = " | ")))
}
fwrite(pheno, file.path("out", paste0("01_pheno_", GSE_ID, ".csv")))

## ---- 3. supplementary clinical table ----------------------------------
# For GSE206285 the routine clinical variables (CRP, calprotectin, Mayo,
# disease duration, prior anti-TNF, neutrophil scores) are NOT in the
# series matrix. They are in a separate supplementary file. Without it
# there is no clinical baseline to beat, so this is not optional.
sup_dir <- file.path("data", GSE_ID)
if (!any(grepl("additional_sample_metadata", list.files(sup_dir), ignore.case = TRUE))) {
  try(getGEOSuppFiles(GSE_ID, baseDir = "data", makeDirectory = TRUE,
                      filter_regex = "additional_sample_metadata"), silent = TRUE)
}
add_path <- list.files("data", pattern = "additional_sample_metadata",
                       full.names = TRUE, recursive = TRUE)[1]

if (is.na(add_path)) {
  cat("\n*** supplementary clinical table NOT downloaded.\n")
  cat("    Go to the GEO page for", GSE_ID, ", download the file whose name\n")
  cat("    contains 'additional_sample_metadata', put it anywhere under data/,\n")
  cat("    and rerun this script. The clinical baseline model needs it.\n")
} else {
  add <- read.delim(if (grepl("\\.gz$", add_path)) gzfile(add_path) else add_path,
                    check.names = FALSE, stringsAsFactors = FALSE)
  names(add) <- sub("^Field\\[(.*)\\]$", "\\1", names(add))   # Field[CRP] -> CRP
  names(add)[names(add) == "Accession"] <- "gsm"
  cat(sprintf("\nclinical table: %d rows x %d columns\n", nrow(add), ncol(add)))
  cat("columns:", paste(names(add), collapse = ", "), "\n")
  cat(sprintf("matched to samples: %d of %d\n",
              sum(pheno$gsm %in% add$gsm), nrow(pheno)))
  fwrite(add, file.path("out", paste0("01_clinical_", GSE_ID, ".csv")))
}

## ---- 4. expression matrix --------------------------------------------
ex <- exprs(eset)
qs <- quantile(ex, c(0, .5, .99, 1), na.rm = TRUE)
cat("\nvalue quantiles:", paste(round(qs, 2), collapse = "  "), "\n")
if (max(ex, na.rm = TRUE) > 50) {
  ex <- log2(ex + 1); cat("applied log2 transform\n")
} else {
  cat("already log2 / normalised — not re-normalising\n")
}
ex <- ex[rowSums(is.na(ex)) == 0, , drop = FALSE]

writeLines(c(paste("Series   :", GSE_ID), paste("Platform :", GPL),
             paste("Samples  :", ncol(ex)), paste("Probes   :", nrow(ex)),
             "", "Data processing, verbatim from GEO:",
             unique(as.character(pd$data_processing))),
           file.path("out", paste0("01_processing_", GSE_ID, ".txt")))

## ---- 5. probes to genes ----------------------------------------------
fd <- fData(eset)
sym_col <- grep("gene.?symbol|^symbol$", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
if (is.na(sym_col)) {
  cat("\nfData columns:", paste(colnames(fd), collapse = ", "), "\n")
  stop("No gene symbol column found. Pick one from the list above by hand.")
}
cat("\nsymbol column:", sym_col, "\n")

sym <- as.character(fd[rownames(ex), sym_col])
cat("probes with no symbol     :", sum(is.na(sym) | sym %in% c("", "---")), "\n")
cat("probes mapping to >1 gene :", sum(grepl("///", sym)), "\n")

# Multi-mapping probes are DROPPED, not truncated to their first symbol.
# A probe that hybridises to several genes cannot be attributed to one of
# them, and keeping it silently invents an expression value for whichever
# gene happens to be listed first.
keep <- !is.na(sym) & !sym %in% c("", "---") & !grepl("///", sym)
ex <- ex[keep, , drop = FALSE]; sym <- toupper(trimws(sym[keep]))

# Several probes per gene: keep the most variable one (IQR), which is the
# probe carrying the most information about between-patient differences.
iqrs  <- rowIQRs(ex)
ord   <- order(sym, -iqrs)
ex    <- ex[ord, , drop = FALSE]; sym <- sym[ord]
first <- !duplicated(sym)
ex    <- ex[first, , drop = FALSE]
rownames(ex) <- sym[first]

cat(sprintf("\n%s: %d genes x %d samples\n", GSE_ID, nrow(ex), ncol(ex)))
cat("XIST present (sanity check):", "XIST" %in% rownames(ex), "\n")
stopifnot(identical(colnames(ex), pheno$gsm))

saveRDS(ex, file.path("out", paste0("01_expr_", GSE_ID, ".rds")))
cat("\nwrote out/01_expr_", GSE_ID, ".rds\n", sep = "")
cat("NEXT: 02_build_cohort.R\n")
